#include <cstring>
#include <immintrin.h>
using namespace std;

template <class Type>
class complex {
    Type real;
    Type imga;
    public:
    complex<Type>() { real = 0.; imga = 0.; }                                  // constructors
    complex<Type>(Type r,Type i) { real = r; imga = i; }
    complex<Type>(Type angle) { real = cos(angle); imga = sin(angle); } 
    void operator=(complex<Type> c) {
        real = c.real;
	imga = c.imga;
    }
    complex<Type> operator+(complex<Type> c) {
        complex<Type> cc;
	cc.real = real + c.real;
	cc.imga = imga + c.imga;
	return cc;
    }
    complex<Type> operator-(complex<Type> c) {
        complex<Type> cc;
	cc.real = real - c.real;
	cc.imga = imga - c.imga;
	return cc;
    }
    complex<Type> operator*(complex<Type> c) {
        complex<Type> cc;
	cc.real = real*c.real - imga*c.imga;
	cc.imga = real*c.imga + imga*c.real;
	return cc;
    }
    complex<Type> operator*(Type c) {
        complex<Type> cc;
	cc.real = real*c;
	cc.imga = imga*c;
	return cc;
    }
    complex<Type> operator/(complex<Type> c) {
        complex<Type> cc;
	Type a = c.real*c.real+c.imga*c.imga;
	cc.real = (real*c.real+imga*c.imga)/a;
	cc.imga = (imga*c.real-real*c.imga)/a;
	return cc;
    }
    complex<Type> swap() {
        complex<Type> cc;
	cc.real = imga;
	cc.imga = real;
        return cc;
    }
    complex<Type> turnleft() {
        complex<Type> cc;
	cc.real = -imga;
	cc.imga = real;
	return cc;
    }
    complex<Type> turnright() {
        complex<Type> cc;
	cc.real = imga;
	cc.imga = -real;
	return cc;
    }
    complex<Type> reverse() {
        complex<Type> cc;
	cc.real = -real;
	cc.imga = -imga;
	return cc;
    }
    complex<Type> conjugate() {
        complex<Type> cc;
        cc.real = real;
	cc.imga = -1.*imga;
	return cc;
    }
    void operator+=(complex<Type> c) {
	real += c.real;
	imga += c.imga;
    }
    void operator-=(complex<Type> c) {
	real -= c.real;
	imga -= c.imga;
    }
    void operator*=(complex<Type> c) {
        Type a,b;
	a = real*c.real - imga*c.imga;
	b = real*c.imga + imga*c.real;
	real = a;
	imga = b;
    }
    void operator*=(Type c) {
        real *= c;
	imga *= c;
    }
    void operator/=(complex<Type> c) {
        Type a,u,v;
	a = c.real*c.real+c.imga*c.imga;
        u = (real*c.real+imga*c.imga)/a;
	v = (imga*c.real-real*c.imga)/a;
	real = u;
	imga = v;
    }
    void setzero() {
        real = 0.;
	imga = 0.;
    }
    void setrealimga(Type r,Type i) {
        real = r;
	imga = i;
    }
    void setreal(Type r) { real = r; }
    void setimga(Type i) { imga = i; }
    void setangle(Type angle) { real = cos(angle); imga = sin(angle); }
    Type getreal() {
        return real;
    }
    Type getimga() {
        return imga;
    }
    void print() {
        cout << real << " " << imga << endl;
    }
};


inline __m256d complex_mul_256register(double a0r,double a0i,double a1r,double a1i,double *b) {
    __m256d a_vals = _mm256_setr_pd(a0r,a0i,a1r,a1i);
    __m256d b_vals = _mm256_load_pd(b);
    __m256d c_vals = _mm256_mul_pd(a_vals,b_vals);
    c_vals = _mm256_xor_pd(c_vals,_mm256_setr_pd(0.0,-0.0,0.0,-0.0));
    b_vals = _mm256_permute_pd(b_vals,0b0101);
    a_vals = _mm256_mul_pd(a_vals,b_vals);
    b_vals = _mm256_hadd_pd(c_vals,a_vals);                      // complex product
    return b_vals;
}

inline __m256d complex_mul_256register(double *a,double *b) {
    __m256d a_vals = _mm256_load_pd(a);                                  //  3   2   4  -1
    __m256d b_vals = _mm256_load_pd(b);                                  //  7   5   6  -2
    __m256d c_vals = _mm256_mul_pd(a_vals,b_vals);                       // 21  10  24   2
    c_vals = _mm256_xor_pd(c_vals,_mm256_setr_pd(0.0,-0.0,0.0,-0.0));    // 21 -10  24  -2
    b_vals = _mm256_permute_pd(b_vals,0b0101);                           //  5   7  -2   6
    a_vals = _mm256_mul_pd(a_vals,b_vals);                               // 15  14  -8  -6
    b_vals = _mm256_hadd_pd(c_vals,a_vals);                              // 11  29  22 -14
    return b_vals;
}

inline __m256 complex_mul_256register(float a0r,float a0i,float a1r,float a1i,float a2r,float a2i,float a3r,float a3i,float *b) {
    __m256 a_vals = _mm256_setr_ps(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i);
    __m256 b_vals = _mm256_load_ps(b);
    __m256 c_vals = _mm256_mul_ps(a_vals,b_vals);
    c_vals = _mm256_xor_ps(c_vals,_mm256_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0));
    b_vals = _mm256_permute_ps(b_vals,0b10110001);
    a_vals = _mm256_mul_ps(a_vals,b_vals);
    b_vals = _mm256_hadd_ps(c_vals,a_vals);
    b_vals = _mm256_permute_ps(b_vals,0b11011000);
    return b_vals;
}

inline __m256 complex_mul_256register(float *a,float *b) {
    __m256 a_vals = _mm256_load_ps(a);
    __m256 b_vals = _mm256_load_ps(b);
    __m256 c_vals = _mm256_mul_ps(a_vals,b_vals);
    c_vals = _mm256_xor_ps(c_vals,_mm256_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0));
    b_vals = _mm256_permute_ps(b_vals,0b10110001);
    a_vals = _mm256_mul_ps(a_vals,b_vals);
    b_vals = _mm256_hadd_ps(c_vals,a_vals);
    b_vals = _mm256_permute_ps(b_vals,0b11011000);
    return b_vals;
}











