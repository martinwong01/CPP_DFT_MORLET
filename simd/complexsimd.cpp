#include <iostream>
#include <immintrin.h>
#include <cstdio>
#include "../allocate.h"
#include "../complex.h"

using namespace std;

int main() {
    complex<double> *a;
    complex<double> *b;
    complex<double> *c;
    double cc = -3.;

    a = allocate1D<complex<double>>(4,32);
    b = allocate1D<complex<double>>(4,32);
    c = allocate1D<complex<double>>(4,32);


    a[0].setrealimga(3.,2.);
    a[1].setrealimga(4.,-1.);
    b[0].setrealimga(7.,5.);
    b[1].setrealimga(6.,-2.);


    __m256d a_vals = _mm256_load_pd((double *)&a[0]);                    //  3   2   4  -1
    __m256d b_vals = _mm256_load_pd((double *)&b[0]);                    //  7   5   6  -2
    __m256d c_vals = _mm256_mul_pd(a_vals,b_vals);                        // 21  10  24   2
    c_vals = _mm256_xor_pd(c_vals,_mm256_setr_pd(0.0,-0.0,0.0,-0.0));     // 21  -10  24  -2

  
  
    b_vals = _mm256_permute_pd(b_vals,0b0101);                           //  5   7   -2    6
    a_vals = _mm256_mul_pd(a_vals,b_vals);                               //  15  14  -8   -6
    b_vals = _mm256_hadd_pd(c_vals,a_vals);                              //  11  29   22   -14
    _mm256_store_pd((double *)&c[0],b_vals);


    c[0].print();
    c[1].print();

}

  
  
