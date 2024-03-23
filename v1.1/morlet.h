#include "dft.h"
#ifndef MAXS
    #define MAXS 128
#endif

template <class Type>
void morlet(Type *data,complex<Type> **transform,int N,int S,Type param,Type dx) {
    Type thread_local pi = acos(-1.);
    int thread_local oldN = 0;
    int thread_local oldS = 0;
    Type thread_local oldparam = 0.;
    Type thread_local olddx = 0.;

    int thread_local k,n,s;
    Type thread_local a,b;
    alignas(ALIGN) Type thread_local freq[MAXN];
    alignas(ALIGN) Type thread_local scale[MAXS];
    alignas(ALIGN) Type thread_local wavefunc[MAXS][MAXN];
    alignas(ALIGN) complex<Type> thread_local dft[MAXN];
    alignas(ALIGN) complex<Type> thread_local dft_product[MAXS][MAXN];
    alignas(ALIGN) complex<Type> thread_local datacomplex[MAXN];
#if AVX512F > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m512d,__m512>::type avxtype;
    alignas(ALIGN) avxtype mw01_morlet_a;
    alignas(ALIGN) avxtype mw01_morlet_b;
#if SVML > 0
    alignas(ALIGN) avxtype mw01_morlet_c;
    alignas(ALIGN) avxtype mw01_morlet_d;
    alignas(ALIGN) avxtype mw01_morlet_e;
#endif
#elif AVX > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m256d,__m256>::type avxtype;
    alignas(ALIGN) avxtype mw01_morlet_a;
    alignas(ALIGN) avxtype mw01_morlet_b;
#if SVML > 0
    alignas(ALIGN) avxtype mw01_morlet_c;
    alignas(ALIGN) avxtype mw01_morlet_d;
    alignas(ALIGN) avxtype mw01_morlet_e;
#endif
#endif


    if(N != oldN || S != oldS || param != oldparam || dx != olddx) {
        a = 2.*pi/N/dx;
        b = 1./pow(pi,0.25);
        for(s=0;s<S;s++) scale[s] = pow(2.,s/4.);                                // set scales with dj=0.25
	for(k=0;k<=N/2;k++) freq[k] = a*k;                                       
//	for(;k<N;k++) freq[k] = a*(k-N);
        for(;k<N;k++) freq[k] = 0.;


#if !defined(AVX) || AVX == 0 || SVML == 0
        for(s=0;s<S;s++) {
  	    for(k=0;k<=N/2;k++)
                wavefunc[s][k] = exp(-0.5*pow(scale[s]*freq[k]-param,2.))*b;     // k>N/2 not set, should be left as 0s
	    for(;k<N;k++)
	        wavefunc[s][k] = 0.;
	}
#elif AVX512F == 0
        if constexpr(std::is_same_v<double,Type>) {
            for(s=0;s<S;s++) {
                mw01_morlet_a = _mm256_set1_pd(-0.5);
                mw01_morlet_b = _mm256_set1_pd(scale[s]);
                mw01_morlet_c = _mm256_set1_pd(param);
                mw01_morlet_d = _mm256_set1_pd(b);
                for(k=0;k<=N/2;k+=4) {
                    mw01_morlet_e = _mm256_fmsub_pd(mw01_morlet_b,_mm256_load_pd((double *)&freq[k]),mw01_morlet_c);
                    _mm256_store_pd((double *)&wavefunc[s][k],_mm256_mul_pd(_mm256_exp_pd(_mm256_mul_pd(_mm256_mul_pd(mw01_morlet_e,mw01_morlet_e),mw01_morlet_a)),mw01_morlet_d));
                }
		mw01_morlet_a = _mm256_setzero_pd();
		for(;k<N;k+=4) _mm256_store_pd((double *)&wavefunc[s][k],mw01_morlet_a);
            }
        } else if constexpr(std::is_same_v<float,Type>) {
            for(s=0;s<S;s++) {
                mw01_morlet_a = _mm256_set1_ps(-0.5);
                mw01_morlet_b = _mm256_set1_ps(scale[s]);
                mw01_morlet_c = _mm256_set1_ps(param);
                mw01_morlet_d = _mm256_set1_ps(b);
                for(k=0;k<=N/2;k+=8) {
                    mw01_morlet_e = _mm256_fmsub_ps(mw01_morlet_b,_mm256_load_ps((float *)&freq[k]),mw01_morlet_c);
                    _mm256_store_ps((float *)&wavefunc[s][k],_mm256_mul_ps(_mm256_exp_ps(_mm256_mul_ps(_mm256_mul_ps(mw01_morlet_e,mw01_morlet_e),mw01_morlet_a)),mw01_morlet_d));
                }
		mw01_morlet_a = _mm256_setzero_ps();
		for(;k<N;k+=8) _mm256_store_ps((float *)&wavefunc[s][k],mw01_morlet_a);
            }
        }
#else
        if constexpr(std::is_same_v<double,Type>) {
            for(s=0;s<S;s++) {
                mw01_morlet_a = _mm512_set1_pd(-0.5);
                mw01_morlet_b = _mm512_set1_pd(scale[s]);
                mw01_morlet_c = _mm512_set1_pd(param);
                mw01_morlet_d = _mm512_set1_pd(b);
                for(k=0;k<=N/2;k+=8) {
                    mw01_morlet_e = _mm512_fmsub_pd(mw01_morlet_b,_mm512_load_pd((double *)&freq[k]),mw01_morlet_c);
                    _mm512_store_pd((double *)&wavefunc[s][k],_mm512_mul_pd(_mm512_exp_pd(_mm512_mul_pd(_mm512_mul_pd(mw01_morlet_e,mw01_morlet_e),mw01_morlet_a)),mw01_morlet_d));
                }
		mw01_morlet_a = _mm512_setzero_pd();
		for(;k<N;k+=8) _mm512_store_pd((double *)&wavefunc[s][k],mw01_morlet_a);
            }
        } else if constexpr(std::is_same_v<float,Type>) {
            for(s=0;s<S;s++) {
                mw01_morlet_a = _mm512_set1_ps(-0.5);
                mw01_morlet_b = _mm512_set1_ps(scale[s]);
                mw01_morlet_c = _mm512_set1_ps(param);
                mw01_morlet_d = _mm512_set1_ps(b);
                for(k=0;k<=N/2;k+=16) {
                    mw01_morlet_e = _mm512_fmsub_ps(mw01_morlet_b,_mm512_load_ps((float *)&freq[k]),mw01_morlet_c);
                    _mm512_store_ps((float *)&wavefunc[s][k],_mm512_mul_ps(_mm512_exp_ps(_mm512_mul_ps(_mm512_mul_ps(mw01_morlet_e,mw01_morlet_e),mw01_morlet_a)),mw01_morlet_d));
                }
		mw01_morlet_a = _mm512_setzero_ps();
		for(;k<N;k+=16) _mm512_store_ps((float *)&wavefunc[s][k],mw01_morlet_a);
            }
        }
#endif
    }

    
#if !defined(AVX) || AVX == 0 || AVX512VL == 0    
    for(n=0;n<N;n++) datacomplex[n].setrealimga(data[n],0);
/*
#elif AVX512F == 0
    if constexpr(std::is_same_v<double,Type>) {
        mw01_morlet_b = _mm256_setzero_pd();
        for(n=0;n<N;n+=4) {
	    mw01_morlet_a = _mm256_load_pd((double *)&data[n]);
	    _mm256_store_pd((double *)&datacomplex[n],_mm256_permutex2var_pd(mw01_morlet_a,_mm256_setr_epi64x(0,4,1,5),mw01_morlet_b));
	    _mm256_store_pd((double *)&datacomplex[n+2],_mm256_permutex2var_pd(mw01_morlet_a,_mm256_setr_epi64x(2,6,3,7),mw01_morlet_b));
	}
    } else if constexpr(std::is_same_v<float,Type>) {
        mw01_morlet_b = _mm256_setzero_ps();
        for(n=0;n<N;n+=8) {
	    mw01_morlet_a = _mm256_load_ps((float *)&data[n]);
	    _mm256_store_ps((float *)&datacomplex[n],_mm256_permutex2var_ps(mw01_morlet_a,_mm256_setr_epi32(0,8,1,9,2,10,3,11),mw01_morlet_b));
	    _mm256_store_ps((float *)&datacomplex[n+4],_mm256_permutex2var_ps(mw01_morlet_a,_mm256_setr_epi32(4,12,5,13,6,14,7,15),mw01_morlet_b));
	}
    }
*/
#else
    if constexpr(std::is_same_v<double,Type>) {
        mw01_morlet_b = _mm512_setzero_pd();
        for(n=0;n<N;n+=8) {
	    mw01_morlet_a = _mm512_load_pd((double *)&data[n]);
	    _mm512_store_pd((double *)&datacomplex[n],_mm512_permutex2var_pd(mw01_morlet_a,_mm512_setr_epi64(0,8,1,9,2,10,3,11),mw01_morlet_b));
	    _mm512_store_pd((double *)&datacomplex[n+4],_mm512_permutex2var_pd(mw01_morlet_a,_mm512_setr_epi64(4,12,5,13,6,14,7,15),mw01_morlet_b));
	}
    } else if constexpr(std::is_same_v<float,Type>) {
        mw01_morlet_b = _mm512_setzero_ps();
        for(n=0;n<N;n+=16) {
	    mw01_morlet_a = _mm512_load_ps((float *)&data[n]);
	    _mm512_store_ps((float *)&datacomplex[n],_mm512_permutex2var_ps(mw01_morlet_a,_mm512_setr_epi32(0,16,1,17,2,18,3,19,4,20,5,21,6,22,7,23),mw01_morlet_b));
	    _mm512_store_ps((float *)&datacomplex[n+8],_mm512_permutex2var_ps(mw01_morlet_a,_mm512_setr_epi32(8,24,9,25,10,26,11,27,12,28,13,29,14,30,15,31),mw01_morlet_b));
	}
    }
#endif
    
    
    dft_func<Type>(datacomplex,dft,N,1,1);

    for(s=0;s<S;s++)
#if !defined(AVX) || AVX == 0    
        for(k=0;k<N;k++) dft_product[s][k] = dft[k]*wavefunc[s][k];
#elif AVX512F == 0
        if constexpr(std::is_same_v<double,Type>) {
	    for(k=0;k<N;k+=2) _mm256_store_pd((double *)&dft_product[s][k],_mm256_mul_pd(_mm256_load_pd((double *)&dft[k]),_mm256_setr_pd(wavefunc[s][k],wavefunc[s][k],wavefunc[s][k+1],wavefunc[s][k+1])));
	} else if constexpr(std::is_same_v<float,Type>) {
	    for(k=0;k<N;k+=4) _mm256_store_ps((float *)&dft_product[s][k],_mm256_mul_ps(_mm256_load_ps((float *)&dft[k]),_mm256_setr_ps(wavefunc[s][k],wavefunc[s][k],wavefunc[s][k+1],wavefunc[s][k+1],wavefunc[s][k+2],wavefunc[s][k+2],wavefunc[s][k+3],wavefunc[s][k+3])));
        }
#else
        if constexpr(std::is_same_v<double,Type>) {
	    for(k=0;k<N;k+=4) _mm512_store_pd((double *)&dft_product[s][k],_mm512_mul_pd(_mm512_load_pd((double *)&dft[k]),_mm512_setr_pd(wavefunc[s][k],wavefunc[s][k],wavefunc[s][k+1],wavefunc[s][k+1],wavefunc[s][k+2],wavefunc[s][k+2],wavefunc[s][k+3],wavefunc[s][k+3])));
	} else if constexpr(std::is_same_v<float,Type>) {
	    for(k=0;k<N;k+=8) _mm512_store_ps((float *)&dft_product[s][k],_mm512_mul_ps(_mm512_load_ps((float *)&dft[k]),_mm512_setr_ps(wavefunc[s][k],wavefunc[s][k],wavefunc[s][k+1],wavefunc[s][k+1],wavefunc[s][k+2],wavefunc[s][k+2],wavefunc[s][k+3],wavefunc[s][k+3],wavefunc[s][k+4],wavefunc[s][k+4],wavefunc[s][k+5],wavefunc[s][k+5],wavefunc[s][k+6],wavefunc[s][k+6],wavefunc[s][k+7],wavefunc[s][k+7])));
        }
#endif
    
    for(s=0;s<S;s++) dftinv_func<Type>(dft_product[s],transform[s],N);

    for(s=0;s<S;s++) {
        a = sqrt(2.*pi*scale[s]/dx);
#if !defined(AVX) || AVX == 0
        for(n=0;n<N;n++) transform[s][n] *= a;
#elif AVX512F == 0
        if constexpr(std::is_same_v<double,Type>) {
    	    mw01_morlet_a = _mm256_set1_pd(a);
	    for(n=0;n<N;n+=2) _mm256_store_pd((double *)&transform[s][n],_mm256_mul_pd(_mm256_load_pd((double *)&transform[s][n]),mw01_morlet_a));
	} else if constexpr(std::is_same_v<float,Type>) {
       	    mw01_morlet_a = _mm256_set1_ps(a);
	    for(n=0;n<N;n+=4) _mm256_store_ps((float *)&transform[s][n],_mm256_mul_ps(_mm256_load_ps((float *)&transform[s][n]),mw01_morlet_a));
	}
#else
        if constexpr(std::is_same_v<double,Type>) {
    	    mw01_morlet_a = _mm512_set1_pd(a);
	    for(n=0;n<N;n+=4) _mm512_store_pd((double *)&transform[s][n],_mm512_mul_pd(_mm512_load_pd((double *)&transform[s][n]),mw01_morlet_a));
	} else if constexpr(std::is_same_v<float,Type>) {
       	    mw01_morlet_a = _mm512_set1_ps(a);
	    for(n=0;n<N;n+=8) _mm512_store_ps((float *)&transform[s][n],_mm512_mul_ps(_mm512_load_ps((float *)&transform[s][n]),mw01_morlet_a));
	}
#endif
    }
    oldN = N;
    oldS = S;
    oldparam = param;
    olddx = dx;
}

template <class Type>
Type integral_reconstruction(Type param) {
    Type thread_local sum;
    Type thread_local lowerbnd = 0.0001;
    Type thread_local upperbnd = 1000.;
    Type thread_local dx = lowerbnd;
    int thread_local intervals;
    int thread_local i;
    Type thread_local x;
    Type thread_local pi = acos(-1.);

    intervals = (int)(upperbnd/dx);
    sum = 0.;

    for(i=1;i<=intervals;i++) {
        x = i*dx;
        sum += exp(-0.5*(x-param)*(x-param))/x;
    }
    sum *= pow(pi,-0.25)*dx;
    return sum;
}

template <class Type>
Type integral_covariance(Type param) {
    Type thread_local sum;
    Type thread_local lowerbnd = 0.0001;
    Type thread_local upperbnd = 1000.;
    Type thread_local dx = lowerbnd;
    int thread_local intervals;
    int thread_local i;
    Type thread_local x;
    Type thread_local pi = acos(-1.);

    intervals = (int)(upperbnd/dx);
    sum = 0.;

    for(i=1;i<=intervals;i++) {
        x = i*dx;
        sum += exp(-(x-param)*(x-param))/x;
    }
    sum *= pow(pi,-0.5)*dx;
    return sum;
}

