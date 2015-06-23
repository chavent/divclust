#include <Rcpp.h>
#include <cstdio>
using namespace Rcpp;

NumericVector int2part(int n, int v) { 
    NumericVector x;
    for(int i = n-1; i >= 0; i--) 
        if (('0'+ ((v >> i) & 1)) == '1') 
            x.push_front(n-i);
    return x;  
}

// [[Rcpp::export]]
std::vector<NumericVector> bipart(int n) {
    std::vector<NumericVector> parts;
    for(int i = 1; i < (1 << (n-1)); ++i) {
        parts.push_back(int2part(n, ~i + (2 << (n-1)) ));
    }
    return parts;
}
