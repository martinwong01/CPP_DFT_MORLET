#include <iostream>
#include <cmath>
#if AVX > 0 || AVX512F > 0 || AVX512VL > 0 || FMA > 0
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
    complex<Type>(Type r) { real = r; imga = 0.; } 
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
    complex<Type> realconjugate() {
        complex<Type> cc;
        cc.real = -real;
	cc.imga = imga;
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
    friend ostream &operator<< (ostream &s,complex<Type> &c) {
        if(c.imga >= 0)
            s << c.real << "+" << c.imga << "i";
        else
            s << c.real << c.imga << "i";
        return s;
    }
};




#if AVX512F > 0
// 512 bit double multiply
inline __m512d complex_mul_512register(__m512d mw01_mul512d_a,__m512d mw01_mul512d_b) {
    return _mm512_fmaddsub_pd(mw01_mul512d_a,_mm512_permute_pd(mw01_mul512d_b,0b00000000),_mm512_mul_pd(_mm512_permute_pd(mw01_mul512d_a,0b01010101),_mm512_permute_pd(mw01_mul512d_b,0b11111111)));
}

inline __m512d complex_mul_512register(double *a,double *b,int aligna,int alignb) {
    if(aligna == 0 && alignb == 0)	
        return complex_mul_512register(_mm512_load_pd(a),_mm512_load_pd(b));
    else if(aligna == 0 && alignb != 0)
        return complex_mul_512register(_mm512_load_pd(a),_mm512_loadu_pd(b));
    else if(aligna != 0 && alignb == 0)
        return complex_mul_512register(_mm512_loadu_pd(a),_mm512_load_pd(b));
    else
        return complex_mul_512register(_mm512_loadu_pd(a),_mm512_loadu_pd(b));
}

inline __m512d complex_mul_512register(double a0r,double a0i,double a1r,double a1i,double a2r,double a2i,double a3r,double a3i,double *b,int alignb) {
    if(alignb == 0)
        return complex_mul_512register(_mm512_setr_pd(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i),_mm512_load_pd(b));	    
    else
        return complex_mul_512register(_mm512_setr_pd(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i),_mm512_loadu_pd(b));
}



// 512 bit float multiply
inline __m512 complex_mul_512register(__m512 mw01_mul512f_a,__m512 mw01_mul512f_b) {
    return _mm512_fmaddsub_ps(mw01_mul512f_a,_mm512_permute_ps(mw01_mul512f_b,0b10100000),_mm512_mul_ps(_mm512_permute_ps(mw01_mul512f_a,0b10110001),_mm512_permute_ps(mw01_mul512f_b,0b11110101)));
}

inline __m512 complex_mul_512register(float *a,float *b,int aligna,int alignb) {
    if(aligna == 0 && alignb == 0)
        return complex_mul_512register(_mm512_load_ps(a),_mm512_load_ps(b));
    else if(aligna == 0 && alignb != 0)
        return complex_mul_512register(_mm512_load_ps(a),_mm512_loadu_ps(b));
    else if(aligna != 0 && alignb == 0)
        return complex_mul_512register(_mm512_loadu_ps(a),_mm512_load_ps(b));
    else
        return complex_mul_512register(_mm512_loadu_ps(a),_mm512_loadu_ps(b));
}

inline __m512 complex_mul_512register(float a0r,float a0i,float a1r,float a1i,float a2r,float a2i,float a3r,float a3i,float a4r,float a4i,float a5r,float a5i,float a6r,float a6i,float a7r,float a7i,float *b,int alignb) {
    if(alignb == 0)
        return complex_mul_512register(_mm512_setr_ps(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i,a4r,a4i,a5r,a5i,a6r,a6i,a7r,a7i),_mm512_load_ps(b));
    else
        return complex_mul_512register(_mm512_setr_ps(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i,a4r,a4i,a5r,a5i,a6r,a6i,a7r,a7i),_mm512_loadu_ps(b));
}
#endif




#if AVX > 0
// 256 bit double multiply
inline __m256d complex_mul_256register(__m256d mw01_mul256d_a,__m256d mw01_mul256d_b) {
    return _mm256_fmaddsub_pd(mw01_mul256d_a,_mm256_permute_pd(mw01_mul256d_b,0b0000),_mm256_mul_pd(_mm256_permute_pd(mw01_mul256d_a,0b0101),_mm256_permute_pd(mw01_mul256d_b,0b1111)));
}

inline __m256d complex_mul_256register(double a0r,double a0i,double a1r,double a1i,double *b,int alignb) {
    if(alignb == 0)
        return complex_mul_256register(_mm256_setr_pd(a0r,a0i,a1r,a1i),_mm256_load_pd(b));
    else
	return complex_mul_256register(_mm256_setr_pd(a0r,a0i,a1r,a1i),_mm256_loadu_pd(b));
}

inline __m256d complex_mul_256register(double *a,double *b,int aligna,int alignb) {
    if(aligna == 0 && alignb == 0)
        return complex_mul_256register(_mm256_load_pd(a),_mm256_load_pd(b));
    else if(aligna == 0 && alignb != 0)
        return complex_mul_256register(_mm256_load_pd(a),_mm256_loadu_pd(b));
    else if(aligna != 0 && alignb == 0)
        return complex_mul_256register(_mm256_loadu_pd(a),_mm256_load_pd(b));
    else
        return complex_mul_256register(_mm256_loadu_pd(a),_mm256_loadu_pd(b));
}




// 256 bit float multiply
inline __m256 complex_mul_256register(__m256 mw01_mul256f_a,__m256 mw01_mul256f_b) {
    return _mm256_fmaddsub_ps(mw01_mul256f_a,_mm256_permute_ps(mw01_mul256f_b,0b10100000),_mm256_mul_ps(_mm256_permute_ps(mw01_mul256f_a,0b10110001),_mm256_permute_ps(mw01_mul256f_b,0b11110101)));
}

inline __m256 complex_mul_256register(float *a,float *b,int aligna,int alignb) {
    if(aligna == 0 && alignb == 0)
        return complex_mul_256register(_mm256_load_ps(a),_mm256_load_ps(b));
    else if(aligna == 0 && alignb != 0)
        return complex_mul_256register(_mm256_load_ps(a),_mm256_loadu_ps(b));
    else if(aligna != 0 && alignb == 0)
        return complex_mul_256register(_mm256_loadu_ps(a),_mm256_load_ps(b));
    else
        return complex_mul_256register(_mm256_loadu_ps(a),_mm256_loadu_ps(b));
}

inline __m256 complex_mul_256register(float a0r,float a0i,float a1r,float a1i,float a2r,float a2i,float a3r,float a3i,float *b,int alignb) {
    if(alignb == 0)
        return complex_mul_256register(_mm256_setr_ps(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i),_mm256_load_ps(b));
    else
        return complex_mul_256register(_mm256_setr_ps(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i),_mm256_loadu_ps(b));
}
#endif


