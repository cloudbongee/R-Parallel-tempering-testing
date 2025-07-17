#include <random>
#include <Rcpp.h>
#include <vector>

using namespace std;
using namespace Rcpp;

struct Pair{
    int index;
    double Φ;
};

NumericMatrix initialize_numMat(const vector<double>& sched, int n, int k, int index) {
    NumericMatrix result(n, k);

    const double* source_block = sched.data() + size_t(index) * n * k;
    memcpy( REAL(result),source_block, sizeof(double) * n * k );

    return result;
}

inline __attribute__((always_inline))
double phi_recalculation(const vector<double>& des, double old_phi, int n, int k, int block, int w, int v, int h, double p, double n_choose_2, double inv_n_choose_2, double minus_p_halves){
    double summation = 0.0;

    int min_i = w < v ? w : v;
    int max_i = w < v ? v : w;

    int hn = block + h * n;;
    double v_element = des[hn + v];
    double w_element = des[hn + w];
    int jn;
    double curr_at_h;
    double delta1;
    double delta2;
    double s;
    double dist_w;
    double dist_v;
    double curr_at_col;
    double diff_w;
    double diff_v ;
    double new_dist_w;
    double new_dist_v;

    for(int i = 0; i < min_i ; i++){ // all three for loops are the same.

        curr_at_h = des[hn + i ];
        delta1 =  v_element - curr_at_h;
        delta2 =  w_element - curr_at_h;
        s = (delta1 * delta1) - (delta2 * delta2);

        dist_w = 0.0;
        dist_v = 0.0;
        for(int j = 0 ; j < k ; j ++ ){
            jn = block + j * n;
            curr_at_col = des[jn + i];
            diff_w = des[jn + w] - curr_at_col;
            diff_v = des[jn + v] - curr_at_col;
            dist_w += diff_w*diff_w;
            dist_v += diff_v*diff_v;
        }

        new_dist_w = dist_w + s;
        new_dist_v = dist_v - s;

        double new_dist_w_pow = std::pow(new_dist_w, minus_p_halves);
        double dist_w_pow = std::pow(dist_w, minus_p_halves);
        double new_dist_v_pow = std::pow(new_dist_v, minus_p_halves);
        double dist_v_pow = std::pow(dist_v, minus_p_halves);

        summation += (new_dist_w_pow - dist_w_pow) + (new_dist_v_pow - dist_v_pow);


    }

    for(int i = min_i + 1; i < max_i ; i++){
        curr_at_h = des[hn + i ];
        delta1 =  v_element - curr_at_h;
        delta2 =  w_element - curr_at_h;
        s = (delta1 * delta1) - (delta2 * delta2);

        dist_w = 0.0;
        dist_v = 0.0;
        for(int j = 0 ; j < k ; j ++ ){
            jn = block + j * n;
            curr_at_col = des[jn + i];
            diff_w = des[jn + w] - curr_at_col;
            diff_v = des[jn + v] - curr_at_col;
            dist_w += diff_w*diff_w;
            dist_v += diff_v*diff_v;
        }

        new_dist_w = dist_w + s;
        new_dist_v = dist_v - s;

        double new_dist_w_pow = std::pow(new_dist_w, minus_p_halves);
        double dist_w_pow = std::pow(dist_w, minus_p_halves);
        double new_dist_v_pow = std::pow(new_dist_v, minus_p_halves);
        double dist_v_pow = std::pow(dist_v, minus_p_halves);

        summation += (new_dist_w_pow - dist_w_pow) + (new_dist_v_pow - dist_v_pow);

    }
    for(int i = max_i + 1; i < n ; i++){

        curr_at_h = des[hn + i ];
        delta1 =  v_element - curr_at_h;
        delta2 =  w_element - curr_at_h;
        s = (delta1 * delta1) - (delta2 * delta2);

        dist_w = 0.0;
        dist_v = 0.0;
        for(int j = 0 ; j < k ; j ++ ){
            jn = block + j * n;
            curr_at_col = des[jn + i];
            diff_w = des[jn + w] - curr_at_col;
            diff_v = des[jn + v] - curr_at_col;
            dist_w += diff_w*diff_w;
            dist_v += diff_v*diff_v;
        }

        new_dist_w = dist_w + s;
        new_dist_v = dist_v - s;

        double new_dist_w_pow = std::pow(new_dist_w, minus_p_halves);
        double dist_w_pow = std::pow(dist_w, minus_p_halves);
        double new_dist_v_pow = std::pow(new_dist_v, minus_p_halves);
        double dist_v_pow = std::pow(dist_v, minus_p_halves);

        summation += (new_dist_w_pow - dist_w_pow) + (new_dist_v_pow - dist_v_pow);

    }
    return pow(inv_n_choose_2 * (n_choose_2 * pow(old_phi, (double)p) + summation) , 1.0/p);

}

inline __attribute__((always_inline))
double phi(const double* D, int n, int k, double p){
    double summation = 0.0;
    for(int i = 0 ; i < n - 1 ; i++ ){
        for(int j = i + 1 ; j < n ; j++){
            double dist_sqr = 0.0;
            for(int l = 0; l < k ; l++){
                double diff = D[l * n + i] - D[l * n + j];
                dist_sqr += diff * diff;
            }
            // take the square root and power of -p at the same time
            summation += pow(dist_sqr, -p * 0.5);
        }
    }
    // average and the power of a large enough p
    return pow((2.0 / (n * (n-1))) * summation, 1.0/p);
}


// [[Rcpp::export]]
List maximinLHD_geom(int n, int k, int M, double p, int Nswap, int Nmax, int tolerance, double alpha){

    mt19937_64 RNG_engine{random_device{}()};
    uniform_int_distribution<int> distSwap(0, M-2); // for choosing random swap
    uniform_int_distribution<int> distCol(0, k-1); // for choosing random col
    uniform_int_distribution<int> distRow(0, n-1); // for choosing random row
    uniform_real_distribution<> dist(0.0, 1.0); // for choosing uniformly distributed prob

    ////////////////////// Initialize random design and schedule ///////////////////////

    vector<int> klist(n);

    iota(klist.begin(),klist.end(),1);

    double D0[n * k];

    for(int j = 0; j < k; j++){
        for(int i = n - 1; i >= 0; --i){
            uniform_int_distribution<int> shuffler(0, i);
            swap(klist[i], klist[shuffler(RNG_engine)]);
            D0[j * n + i] = klist[i];
        }
    }

    double d2 = k/6.0 * n * (n+1);
    double δ0 = 1.0/sqrt(d2 - k) - 1/sqrt(d2);
    double T0 = -δ0 * (1.0/log(0.99));

    double Φ0 = phi(D0, n, k, p);

    vector<double> design_sched;
    design_sched.resize(size_t(M) * n * k);

    vector<double> temp_sched;
    vector<double> Φ_sched;

    for(int i = 0; i < M; i++){

        double * predestination = design_sched.data() + size_t(i) * n * k;
        memcpy(predestination, D0, sizeof(double) * n * k);
        // beta decay of the temperatres
        // double new_temp = T0 * exp(-0.95 * i);

        double new_temp = T0 * pow(alpha, M - i);

        temp_sched.push_back(new_temp);
        Φ_sched.push_back(Φ0);
    }

    int iter = 1;
    int no_improv = 0;

    Pair prev_best = Pair();
    prev_best.Φ = Φ0;
    prev_best.index = 0;

    Pair best = Pair();
    Pair curr = Pair();

    double minus_p_halves = (double)-p*0.5;
    double n_choose_2 = 0.5 * n * (n-1);
    double inv_n_choose_2 = 1.0/n_choose_2;


    while( iter < Nmax && no_improv < tolerance )
    {

        best.Φ = Φ_sched[0];
        best.index = 0;

        for( int ii = 0 ; ii < M ; ii++) // metropolis loop
        {
            // random rows and cols.
            int w = distRow(RNG_engine);
            int v;
            do{ v = distRow(RNG_engine); }while(w == v);
            int h = distCol(RNG_engine);

            int block = ii * n * k;
            int hn =  block + h * n;

            double Φ_old = Φ_sched[ii];

            // this is always inlined.
            double Φ_try = phi_recalculation(design_sched, Φ_old, n, k, block, w, v, h, p, n_choose_2,inv_n_choose_2,minus_p_halves);

            if(Φ_try < Φ_old || dist(RNG_engine) <  exp(- (Φ_try - Φ_old)/temp_sched[ii])){
                Φ_old = Φ_try;
                Φ_sched[ii] = Φ_try;
                swap(design_sched[hn + v] , design_sched[hn + w]);

            }

            curr.Φ     = Φ_old;
            curr.index = ii;
            best = Φ_old < best.Φ ? curr : best;

        }

        // is the minimum of the swaps an improvement from
        // the previous best?
        if(best.Φ < prev_best.Φ){
            no_improv = 0;
            prev_best = best;
        };

        if((iter % Nswap) == 0){
            int q = distSwap(RNG_engine);
            double Φ1 = Φ_sched[q];
            double Φ2 = Φ_sched[q+1];
            double t1 = temp_sched[q];
            double t2 = temp_sched[q+1];

            if(dist(RNG_engine) <= exp((Φ2 - Φ1) * (1/t2 - 1/t1))){
                swap(temp_sched[q], temp_sched[q+1]);
            }

        }
        iter++;
        no_improv++;

    }
    
    return List::create(
            Named("design") = initialize_numMat(design_sched,n,k,best.index),
            Named("measure") = best.Φ,
            Named("t0") = temp_sched,
            Named("ntotal") = iter
    );
}