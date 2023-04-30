#include <iostream>
#include <cmath>
#if AVX > 0 || AVX512 > 0 || FMA > 0
    #include <immintrin.h>
#endif
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




#if AVX512 > 0
// 512 bit double multiply
inline __m512d complex_mul_512register(register __m512d mw01_mul512d_a,register __m512d mw01_mul512d_b) {
    return _mm512_fmaddsub_pd(mw01_mul512d_a,_mm512_permute_pd(mw01_mul512d_b,0b00000000),_mm512_mul_pd(_mm512_permute_pd(mw01_mul512d_a,0b01010101),_mm512_permute_pd(mw01_mul512d_b,0b11111111)));
}

inline __m512d complex_mul_512register(double *a,double *b,int aligna,int alignb) {
    register __m512d mw01_mul512d_a; //asm("ymm13");
    register __m512d mw01_mul512d_b; //asm("ymm12");

    if(aligna == 0)
        mw01_mul512d_a = _mm512_load_pd(a);                                      
    else
        mw01_mul512d_a = _mm512_loadu_pd(a);
    if(alignb == 0)
        mw01_mul512d_b = _mm512_load_pd(b);                                      
    else
        mw01_mul512d_b = _mm512_loadu_pd(b);

    return complex_mul_512register(mw01_mul512d_a,mw01_mul512d_b);
}

inline __m512d complex_mul_512register(double a0r,double a0i,double a1r,double a1i,double a2r,double a2i,double a3r,double a3i,double *b,int alignb) {
    register __m512d mw01_mul512d_a; //asm("ymm13");
    register __m512d mw01_mul512d_b; //asm("ymm12");

    mw01_mul512d_a = _mm512_setr_pd(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i);
    if(alignb == 0)
        mw01_mul512d_b = _mm512_load_pd(b);                                      
    else
        mw01_mul512d_b = _mm512_loadu_pd(b);

    return complex_mul_512register(mw01_mul512d_a,mw01_mul512d_b);
}



// 512 bit float multiply
inline __m512 complex_mul_512register(register __m512 mw01_mul512f_a,register __m512 mw01_mul512f_b) {
    return _mm512_fmaddsub_ps(mw01_mul512f_a,_mm512_permute_ps(mw01_mul512f_b,0b10100000),_mm512_mul_ps(_mm512_permute_ps(mw01_mul512f_a,0b10110001),_mm512_permute_ps(mw01_mul512f_b,0b11110101)));
}

inline __m512 complex_mul_512register(float *a,float *b,int aligna,int alignb) {
    register __m512 mw01_mul512f_a; //asm("ymm13");
    register __m512 mw01_mul512f_b; //asm("ymm12");

    if(aligna == 0)
        mw01_mul512f_a = _mm512_load_ps(a);                                      
    else
        mw01_mul512f_a = _mm512_loadu_ps(a);
    if(alignb == 0)
        mw01_mul512f_b = _mm512_load_ps(b);                                      
    else
        mw01_mul512f_b = _mm512_loadu_ps(b);
    return complex_mul_512register(mw01_mul512f_a,mw01_mul512f_b);
}

inline __m512 complex_mul_512register(float a0r,float a0i,float a1r,float a1i,float a2r,float a2i,float a3r,float a3i,float a4r,float a4i,float a5r,float a5i,float a6r,float a6i,float a7r,float a7i,float *b,int alignb) {
    register __m512 mw01_mul512f_a; //asm("ymm13");
    register __m512 mw01_mul512f_b; //asm("ymm12");

    mw01_mul512f_a = _mm512_setr_ps(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i,a4r,a4i,a5r,a5i,a6r,a6i,a7r,a7i);
    if(alignb == 0)
        mw01_mul512f_b = _mm512_load_ps(b);                                      
    else
        mw01_mul512f_b = _mm512_loadu_ps(b);
    return complex_mul_512register(mw01_mul512f_a,mw01_mul512f_b);
}
#endif




#if AVX > 0
// 256 bit double multiply
inline __m256d complex_mul_256register(register __m256d mw01_mul256d_a,register __m256d mw01_mul256d_b) {
/*
                                                                          a:    3   2   4  -1
                                                                          b:    7   5   6  -2
    c_256_2 = _mm256_mul_pd(a_256_2,b_256_2);                                  21  10  24   2
    c_256_2 = _mm256_xor_pd(c_256_2,_mm256_setr_pd(0.0,-0.0,0.0,-0.0));        21 -10  24  -2
    b_256_2 = _mm256_permute_pd(b_256_2,0b0101);                                5   7  -2   6
    a_256_2 = _mm256_mul_pd(a_256_2,b_256_2);                                  15  14  -8  -6
    d_256_2 = _mm256_hadd_pd(c_256_2,a_256_2);                                 11  29  22 -14
    return d_256_2;
    return _mm256_hadd_pd(_mm256_xor_pd(_mm256_mul_pd(mw01_mul256d_a,mw01_mul256d_b),_mm256_setr_pd(0.0,-0.0,0.0,-0.0)),_mm256_mul_pd(mw01_mul256d_a,_mm256_permute_pd(mw01_mul256d_b,0b0101)));
*/
	
    return _mm256_fmaddsub_pd(mw01_mul256d_a,_mm256_permute_pd(mw01_mul256d_b,0b0000),_mm256_mul_pd(_mm256_permute_pd(mw01_mul256d_a,0b0101),_mm256_permute_pd(mw01_mul256d_b,0b1111)));
}

inline __m256d complex_mul_256register(double a0r,double a0i,double a1r,double a1i,double *b,int alignb) {
    register __m256d mw01_mul256d_a asm("ymm15");
    register __m256d mw01_mul256d_b asm("ymm14");

    mw01_mul256d_a = _mm256_setr_pd(a0r,a0i,a1r,a1i);
    if(alignb == 0)
        mw01_mul256d_b = _mm256_load_pd(b);
    else
        mw01_mul256d_b = _mm256_loadu_pd(b);

    return complex_mul_256register(mw01_mul256d_a,mw01_mul256d_b);
}

inline __m256d complex_mul_256register(double *a,double *b,int aligna,int alignb) {
    register __m256d mw01_mul256d_a asm("ymm13");
    register __m256d mw01_mul256d_b asm("ymm12");

    if(aligna == 0)
        mw01_mul256d_a = _mm256_load_pd(a);                                      
    else
        mw01_mul256d_a = _mm256_loadu_pd(a);
    if(alignb == 0)
        mw01_mul256d_b = _mm256_load_pd(b);                                      
    else
        mw01_mul256d_b = _mm256_loadu_pd(b);

    return complex_mul_256register(mw01_mul256d_a,mw01_mul256d_b);
}




// 256 bit float multiply
inline __m256 complex_mul_256register(register __m256 mw01_mul256f_a,register __m256 mw01_mul256f_b) {
    return _mm256_fmaddsub_ps(mw01_mul256f_a,_mm256_permute_ps(mw01_mul256f_b,0b10100000),_mm256_mul_ps(_mm256_permute_ps(mw01_mul256f_a,0b10110001),_mm256_permute_ps(mw01_mul256f_b,0b11110101)));
}

inline __m256 complex_mul_256register(float *a,float *b,int aligna,int alignb) {
    register __m256 mw01_mul256f_a; //asm("ymm13");
    register __m256 mw01_mul256f_b; //asm("ymm12");

    if(aligna == 0)
        mw01_mul256f_a = _mm256_load_ps(a);
    else
        mw01_mul256f_a = _mm256_loadu_ps(a);
    if(alignb == 0)
        mw01_mul256f_b = _mm256_load_ps(b);
    else
        mw01_mul256f_b = _mm256_loadu_ps(b);
    return complex_mul_256register(mw01_mul256f_a,mw01_mul256f_b);
}

inline __m256 complex_mul_256register(float a0r,float a0i,float a1r,float a1i,float a2r,float a2i,float a3r,float a3i,float *b,int alignb) {
    register __m256 mw01_mul256f_a asm("ymm15");
    register __m256 mw01_mul256f_b asm("ymm14");

    mw01_mul256f_a = _mm256_setr_ps(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i);
    if(alignb == 0)
        mw01_mul256f_b = _mm256_load_ps(b);
    else
        mw01_mul256f_b = _mm256_loadu_ps(b);
    return complex_mul_256register(mw01_mul256f_a,mw01_mul256f_b);
}
#endif


