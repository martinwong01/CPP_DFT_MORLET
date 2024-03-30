#include <iostream>
#include <cmath>
#if AVX > 0 || AVX512F > 0 || AVX512VL > 0 || FMA > 0
    #include <immintrin.h>
#endif
using namespace std;

template <class Type>
class complex {
    Type real;
    Type imag;
  public:
    constexpr complex<Type>() { real = 0.; imag = 0.; }                                  
    constexpr complex<Type>(const Type& r,const Type& i) { real = r; imag = i; }
    constexpr complex<Type>(const Type& r) { real = r; imag = 0.; } 

    constexpr complex(const complex<Type>& c) {                            
        real = c.real;
        imag = c.imag;
    }

    constexpr complex(complex<Type>&& c) {                   
        if(this == &c) return;
        real = c.real;
        imag = c.imag;
    }

    constexpr complex<Type>& operator=(complex<Type>&& c) {      
        if(this == &c) return *this;
        real = c.real;                                          
        imag = c.imag;
        return *this;
    }

    constexpr complex<Type>& operator=(const complex<Type>& c) {      
        if(this == &c) return *this;
        real = c.real;
        imag = c.imag;
        return *this;
    }

    constexpr complex<Type> operator+(const complex<Type>& c) {
        complex<Type> temp;
        temp.real = real + c.real;
        temp.imag = imag + c.imag;
        return temp;
    }
    constexpr complex<Type> operator-(const complex<Type>& c) {
        complex<Type> temp;
        temp.real = real - c.real;
        temp.imag = imag - c.imag;
        return temp;
    }
    constexpr complex<Type> operator*(const complex<Type>& c) {
        complex<Type> temp;
        temp.real = real*c.real-imag*c.imag;
        temp.imag = real*c.imag+imag*c.real;
        return temp;
    }
    constexpr complex<Type> operator*(const Type& a) {
        complex<Type> temp;
        temp.real = real*a;
        temp.imag = imag*a;
        return temp;
    }
    template <typename Type2>
    friend constexpr complex<Type2> operator*(const Type2& a,const complex<Type2>&);
    constexpr complex<Type> operator/(const complex<Type>& c) {
        complex<Type> cc;
	Type a = c.real*c.real+c.imag*c.imag;
	cc.real = (real*c.real+imag*c.imag)/a;
	cc.imag = (imag*c.real-real*c.imag)/a;
	return cc;
    }
    constexpr complex<Type> operator/(const Type& a) {
        complex<Type> c;
        c.real = real/a;
        c.imag = imag/a;
        return c; 
    }
    template <typename Type2>
    friend constexpr complex<Type2> operator/(const Type2&,const complex<Type2>&);
    template <typename Type2>
    friend constexpr complex<Type2> swapcomplex(const complex<Type2>&);
    template <typename Type2>
    friend constexpr complex<Type2> turnleft(const complex<Type2>&);
    template <typename Type2>
    friend constexpr complex<Type2> turnright(const complex<Type2>&);
    template <typename Type2>
    friend constexpr complex<Type2> reverse(const complex<Type2>&);
    template <typename Type2>
    friend constexpr complex<Type2> conj(const complex<Type2>&);
    template <typename Type2>
    friend constexpr complex<Type2> realconj(const complex<Type2>&);
    void operator+=(const complex<Type>& c) {
	real += c.real;
	imag += c.imag;
    }
    void operator-=(const complex<Type>& c) {
	real -= c.real;
	imag -= c.imag;
    }
    void operator*=(const complex<Type>& c) {
        Type a,b;
	a = real*c.real - imag*c.imag;
	b = real*c.imag + imag*c.real;
	real = a;
	imag = b;
    }
    void operator*=(const Type& a) {
        real *= a;
	imag *= a;
    }
    void operator/=(const complex<Type>& c) {
        Type a,u,v;
	a = c.real*c.real+c.imag*c.imag;
        u = (real*c.real+imag*c.imag)/a;
	v = (imag*c.real-real*c.imag)/a;
	real = u;
	imag = v;
    }
    void operator/=(const Type& a) {
        real /= a;
        imag /= a;
    }
    bool operator==(const complex<Type>& c) {
        if(real == c.real && imag == c.imag) return true; else return false;
    }
    bool operator!=(const complex<Type>& c) {
        if(real != c.real || imag != c.imag) return true; else return false;
    }
    template <typename Type2>
    friend constexpr Type2 real(const complex<Type2>&);
    template <typename Type2>
    friend constexpr Type2 imag(const complex<Type2>&);
    template <typename Type2>
    friend constexpr Type2 fabs(const complex<Type2>&);
    template <typename Type2>
    friend constexpr Type2 norm(const complex<Type2>&);
    template <typename Type2>
    friend ostream &operator<< (ostream &,const complex<Type2>&);
};

template <typename Type2>
constexpr complex<Type2> operator*(const Type2& a,const complex<Type2>& c) {
    complex<Type2> cc;
    cc.real = c.real*a;
    cc.imag = c.imag*a;
    return cc;    
}

template <typename Type2>
constexpr complex<Type2> operator/(const Type2& a,const complex<Type2>& c) {
    complex<Type2> cc;
    Type2 b = c.real*c.real+c.imag*c.imag;
    cc.real = a/b*c.real;
    cc.imag = -a/b*c.imag;
    return cc;
}

template <typename Type2>
constexpr complex<Type2> swapcomplex(const complex<Type2>& c) {
    complex<Type2> cc;
    cc.real = c.imag;
    cc.imag = c.real;
    return cc;
}

template <typename Type2>
constexpr complex<Type2> turnleft(const complex<Type2>& c) {
    complex<Type2> cc;
    cc.real = -c.imag;
    cc.imag = c.real;
    return cc;
}

template <typename Type2>
constexpr complex<Type2> turnright(const complex<Type2>& c) {
    complex<Type2> cc;
    cc.real = c.imag;
    cc.imag = -c.real;
    return cc;
}

template <typename Type2>
constexpr complex<Type2> reverse(const complex<Type2>& c) {
    complex<Type2> cc;
    cc.real = -c.imag;
    cc.imag = -c.real;
    return cc;
}

template <typename Type2>
constexpr complex<Type2> conj(const complex<Type2>& c) {
    complex<Type2> cc;
    cc.real = c.real;
    cc.imag = -c.imag;
    return cc;
}

template <typename Type2>
constexpr complex<Type2> realconj(const complex<Type2>& c) {
    complex<Type2> cc;
    cc.real = -c.real;
    cc.imag = c.imag;
    return cc;
}

template <typename Type2>
constexpr Type2 real(const complex<Type2>& c) {
    return c.real;
}

template <typename Type2>
constexpr Type2 imag(const complex<Type2>& c) {
    return c.imag;
}

template<class Type2>
constexpr Type2& real(complex<Type2>& c) {
    return reinterpret_cast<Type2*>(&c)[0];
}

template<class Type2>
constexpr Type2& imag(complex<Type2>& c) {
  return reinterpret_cast<Type2*>(&c)[1];
}

template <typename Type2>
constexpr Type2 fabs(const complex<Type2>& c) {
    Type2 a;
    a = sqrt(c.real*c.real+c.imag*c.imag);
    return a;
}

template <typename Type2>
constexpr Type2 norm(const complex<Type2>& c) {
    Type2 a;
    a = c.real*c.real+c.imag*c.imag;
    return a;
}

template <typename Type2>
ostream &operator<< (ostream &s,const complex<Type2>& c) {
    if(c.imag >= 0)
        s << c.real << "+" << c.imag << "i";
    else
        s << c.real << c.imag << "i";
    return s;
}



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

inline __m512d complex_mul_512register(double& a0r,double& a0i,double& a1r,double& a1i,double& a2r,double& a2i,double& a3r,double& a3i,double *b,int alignb) {
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

inline __m512 complex_mul_512register(float& a0r,float& a0i,float& a1r,float& a1i,float& a2r,float& a2i,float& a3r,float& a3i,float& a4r,float& a4i,float& a5r,float& a5i,float& a6r,float& a6i,float& a7r,float& a7i,float *b,int alignb) {
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

inline __m256d complex_mul_256register(double& a0r,double& a0i,double& a1r,double& a1i,double *b,int alignb) {
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

inline __m256 complex_mul_256register(float& a0r,float& a0i,float& a1r,float& a1i,float& a2r,float& a2i,float& a3r,float& a3i,float *b,int alignb) {
    if(alignb == 0)
        return complex_mul_256register(_mm256_setr_ps(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i),_mm256_load_ps(b));
    else
        return complex_mul_256register(_mm256_setr_ps(a0r,a0i,a1r,a1i,a2r,a2i,a3r,a3i),_mm256_loadu_ps(b));
}
#endif


