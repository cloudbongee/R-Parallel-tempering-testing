
#include<Rcpp.h>
using namespace Rcpp;

// Criterion utils

// [[Rcpp::export]]
double phi(NumericMatrix D, int n, int k, double p){
    double summation = 0.0;
    for(int i = 0 ; i < n - 1 ; i++ ){
        for(int j = i + 1 ; j < n ; j++){
            double dist = 0.0;
            double diff;
            for(int l = 0; l < k ; l++){
                diff = (D(i,l) - D(j,l));
                dist += diff * diff;
            }
            summation += pow(dist , -p * 0.5);
        }
    }
    // average and the power of a large enough p
    return pow((2.0 / (n * (n-1))) * summation, 1.0/p);
}

// [[Rcpp::export]]
double psi(NumericMatrix D, int n, int k){
    double summation = 0.0;
    for(int i = 0 ; i < n - 1 ; i++ ){
        for(int j = i + 1 ; j < n ; j++){
            double prod = 1.0;
            for(int l = 0; l < k; l++)
            {
                double diff=  (D(i,l) - D(j,l));
                prod *= diff * diff;
            }
            summation += 1.0/prod;
        }
    }
    // average and the power of a large enough p
    return pow(2.0/(n * (n - 1)) * summation, (double)1.0/k);
}