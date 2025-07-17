
#include <Rcpp.h>
#include <vector>
#include <random>
using namespace std;
using namespace Rcpp;



// [[Rcpp::export]]
NumericMatrix randomLHD(int n, int k){
    mt19937_64 RNG_engine{random_device{}()};


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
            D0[j * n + i] = klist[i];
        }
    }

    NumericMatrix result(n, k);

    memcpy( REAL(result), D0, sizeof(double) * n * k );

    return result;

    
}