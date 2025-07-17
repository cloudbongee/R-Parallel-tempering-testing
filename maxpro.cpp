
#include <Rcpp.h>
#include <random>

using namespace std;
using namespace Rcpp;

// maxpro parallel tempering
// Jaime Meyer Beilis Michel
// 7.4.2025 at Bloomington, Indiana

// Initialize the numeric matrix object.
NumericMatrix initialize_numMat(const vector<double>& sched, int n, int k, int index) {
    NumericMatrix result(n, k);

    const double* source_block = sched.data() + std::size_t(index) * n * k;
    std::memcpy( REAL(result),source_block, sizeof(double) * n * k );

    return result;
}

// Calculate the measure (Energy) of a design D,
// by the average interpoint distance reciprocal
// based on a uniform bayesian prior of the
// importance weight of each sub-dimension
// O(n^2 k) time.

double psi(const double* D, int n, int k){
    double summation = 0.0;
    for(int i = 0 ; i < n - 1 ; i++ ){
        for(int j = i + 1 ; j < n ; j++){
            double prod = 1.0;
            for(int l = 0; l < k; l++)
            {
                double diff=  (D[l * n + i] - D[l * n + j]);
                prod *= diff * diff;
            }
            summation += 1.0/prod;
        }
    }
    // average and the power of a large enough p
    return pow(2.0/(n * (n - 1)) * summation, (double)1.0/k);
}



// [[Rcpp::export]]
List pt_maxpro_lhd(int n, int k, int M, int Nmax, int Nswap, int tolerance,  NumericVector temp_sched){
    
    // start the RNG
    mt19937_64 RNG_engine{random_device{}()};
    uniform_int_distribution<int> distSwap(0, M-2); // for choosing random swap
    uniform_int_distribution<int> distCol(0, k-1); // for choosing random col
    uniform_int_distribution<int> distRow(0, n-1); // for choosing random row
    uniform_real_distribution<> dist(0.0, 1.0); // for choosing uniformly distributed prob

    ////////////////////// Initialize random design ///////////////////////

    vector<int> klist(n);
    // fill 1:k
    iota(klist.begin(),klist.end(),1);

    // Column major matrix
    // element i,j (row/column) are accessed as D0[j * n + i]
    double D0[n * k]; 

    for(int j = 0; j < k; j++){
        // shuffle klist
        for(int i = n - 1; i >= 0; --i){
            // shuffle the integers
            uniform_int_distribution<int> shuffler(0, i);
            swap(klist[i], klist[shuffler(RNG_engine)]);
            D0[j * n + i] = (klist[i] - 0.5)/n;
        }
    }
    
    
    /////////////////////// End random design initialization //////////////

    
    ///////////////////// Calculate initial psi //////////////////////////

    double psi0 = psi(D0, n, k);

    //////////////////// Make full replica schedule //////////////////////
    
    vector<double> design_sched;
    // allocate M nk doubles
    design_sched.resize(size_t(M) * n * k);

    vector<double> psi_sched;

    // Fill in lists 1 through M of design configurations
    // psi is kept here to avoid unnecessary recomputation
    for(int i = 0; i < M; i++){

        double * predestination = design_sched.data() + size_t(i) * n * k;

        memcpy(predestination, D0, sizeof(double) * n * k);

        psi_sched.push_back(psi0);
    }

    //////////////////// End of replica schedule generation ////////////

    // convergence testing
    int unsuccessful_global_min = 0;
    int global_min_index = 0;

    // since here, psi0 will be utilized for global minimum. 
    // To ensure (somewhat) good use of space.

    ///////////////////// Run the PT algorithm ////////////////////
    for(int iter = 0; iter < Nmax; iter++){

        ///////////// Metropolis algorithm ///////////////////////
        for(int ii = 0; ii < M; ii++){
            
            // double* D = design_sched[ii * n * k + i * k + j];

            // choose perturbation rows and column
            int h = distCol(RNG_engine);
            int w = distRow(RNG_engine);
            int v;
            do { // ensure distinction between indices
                v = distRow(RNG_engine);
            } while (v == w);

            // recalculate psi
            double psi_original = psi_sched[ii];



            //////////////////// psi RECALCULATION ///////////////////////////////

            // keep track of the memory block of the matrix
            int block = (ii * n * k);

            double summation = 0.0;
            for(int i = 0; i < n; i++){
                if(!(i == w || i == v)){
                    double prod1 = 1.0;
                    double prod2 = 1.0;
                    for(int l = 0; l < k ; l++)
                    {
                        double variable = design_sched[block + l * n + i];
                        double diff1 = (design_sched[block + l * n + w] - variable);
                        double diff2 = (design_sched[block + l * n + v] - variable);

                        prod1 *= diff1 * diff1;
                        prod2 *= diff2 * diff2;
                    }

                    double affected_col = design_sched[block + h * n + i];
                    double delta_v = design_sched[block + h * n + v] - affected_col;
                    double delta_w = design_sched[block + h * n + w] - affected_col;
                    double replacement = (delta_v * delta_v)/(delta_w * delta_w);
                    summation += 
                    (
                        // w term
                        1.0/(prod1 * replacement) - 1.0/prod1 
                        +
                        // v term
                        1.0/(prod2 * (1.0/replacement))  - 1.0/prod2

                    );
                }
            }

            double psi_try = pow(2.0/(n *(n-1)) * ((n * (n-1) * 0.5) * pow(psi_original, (double)k) + summation), (double)1.0/k);
            
            
            /////////////////// psi RECALCULATION END ////////////////////////////////


            // if the recalculated psi try is better than the original
            // or alternatively the probability of acceptance is high enough 
            // compared to a uniformly distributed random probability
            if(
                psi_try < psi_original 
            )
            {
                // accept the perturbation of the design
                swap(design_sched[block +  h * n + w], design_sched[block + h * n + v]);
                // the psi designated becomes that of the perturbed design
                psi_sched[ii] = psi_try;

                // check if psi_try is the global minimum
                if(psi_try < psi0){
                    psi0 = psi_try;
                    unsuccessful_global_min = 0;
                    global_min_index = ii;
                }

            }else if(dist(RNG_engine) <= exp(-(psi_try - psi_original)/temp_sched[ii])){
                swap(design_sched[block +  h * n + w], design_sched[block + h * n + v]);
                psi_sched[ii] = psi_try;
                if(ii == global_min_index){
                    psi0 = psi_try;
                }
            }

        }

        ////////////// Swap configuration ///////////////////////
        if((iter + 1) % Nswap == 0){

            int q = distSwap(RNG_engine);
            double psi1 = psi_sched[q];
            double psi2 = psi_sched[q+1];
            double t1 = temp_sched[q];
            double t2 = temp_sched[q+1];

            // Since we are keeping three lists instead of 1, Decide on swapping
            // the temperatures. If I am not wrong, it should be mathematically equivalent
            // to swapping the designs. Regardless of having broken the original order
            // of T1 < T2 < T3 < ... < TM, just an irrelevant -1 swap to worry about.
            if(dist(RNG_engine) <= exp((psi2 - psi1) * (1/t2 - 1/t1))){
                swap(temp_sched[q], temp_sched[q+1]);
            }

        }
        /////////////// End of configuration swap ///////////////


        // unsuccesful_global_minimum check
        if(++unsuccessful_global_min > tolerance){
            Nmax = iter;
            break;
        }
    }

    //////////////////// End of the PT algorithm //////////////////
    return List::create(
            Named("design") = initialize_numMat(design_sched,n,k,global_min_index),
            Named("measure") = psi0,
            Named("t0") = temp_sched,
            Named("ntotal") = Nmax
    );
}
// End of code, why are you still here?