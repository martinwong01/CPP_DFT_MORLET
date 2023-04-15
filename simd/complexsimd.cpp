/////    double complex, 256bit register

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
    complex<float> *af;
    complex<float> *bf;
    complex<float> *cf;

    double d[4];
    double e[4];

    float df[8];
    float ef[8];

    a = allocate1D<complex<double>>(4,32);
    b = allocate1D<complex<double>>(4,32);
    c = allocate1D<complex<double>>(4,32);
    af = allocate1D<complex<float>>(4,32);
    bf = allocate1D<complex<float>>(4,32);
    cf = allocate1D<complex<float>>(4,32);

    a[0].setrealimga(3.,2.);
    a[1].setrealimga(4.,-1.);
    b[0].setrealimga(7.,5.);
    b[1].setrealimga(6.,-2.);
    af[0].setrealimga(3.,2.);
    af[1].setrealimga(4.,-1.);
    af[2].setrealimga(-5.,5.);
    af[3].setrealimga(-1.,4.);
    bf[0].setrealimga(7.,5.);
    bf[1].setrealimga(6.,-2.);
    bf[2].setrealimga(2.,2.);
    bf[3].setrealimga(8.,-7.);

/*
//    __m256d a_vals = complex_mul_256register((double *)&a[0],(double *)&b[0]);
    __m256d a_vals = complex_mul_256register(a[0].getreal(),a[0].getimga(),a[1].getreal(),a[1].getimga(),(double *)&b[0]);
    _mm256_store_pd((double *)&c[0],a_vals);
    c[0].print();
    c[1].print();
*/

    __m256 e_vals = complex_mul_256register(af[0].getreal(),af[0].getimga(),af[1].getreal(),af[1].getimga(),af[2].getreal(),af[2].getimga(),af[3].getreal(),af[3].getimga(),(float *)&bf[0]);
//    __m256 e_vals = complex_mul_256register((float *)&af[0],(float *)&bf[0]);
    _mm256_store_ps((float *)&cf[0],e_vals);
    cf[0].print();
    cf[1].print();
    cf[2].print();
    cf[3].print();
}
