#include "dft.h"
#ifndef MAXS
    #define MAXS 128
#endif

template <class Type>
void morlet(Type *data,complex<Type> **transform,int N,int S,Type param,Type dx,Type pi,int init) {
    int thread_local k,n,s;
    Type thread_local a,b;
    alignas(ALIGN) Type thread_local freq[MAXN];
    alignas(ALIGN) Type thread_local scale[MAXS];
    alignas(ALIGN) Type thread_local wavefunc[MAXS][MAXN];
    alignas(ALIGN) complex<Type> thread_local dft[MAXN];
    alignas(ALIGN) complex<Type> thread_local dft_product[MAXS][MAXN];
    alignas(ALIGN) complex<Type> thread_local datacomplex[MAXN];
#if AVX512 > 0
    alignas(ALIGN) __m512d thread_local mw01_morlet_a;
    alignas(ALIGN) __m512 thread_local mw01_morlet_af;
#elif AVX > 0
    alignas(ALIGN) __m256d thread_local mw01_morlet_a;
    alignas(ALIGN) __m256 thread_local mw01_morlet_af;
#endif


    if(init == 1) {
        a = 2.*pi/N/dx;
        b = 1./pow(pi,0.25);
        for(s=0;s<S;s++) scale[s] = pow(2.,s/4.);                                        // set scales with dj=0.25
        for(k=0;k<N;k++)
            if(k <= N/2)
                freq[k] = a*k;
	    else
                freq[k] = a*(k-N); 
/*
        for(s=0;s<S;s++)
        for(k=0;k<N;k++)
            if(freq[k] > 0.)
                wavefunc[s][k] = exp(-0.5*pow(scale[s]*freq[k]-param,2.))*b;
	    else
	        wavefunc[s][k] = exp(-0.5*pow(scale[s]*freq[k]-param,2.))*b;
	        //wavefunc[s][k] = 0.;
*/
        for(s=0;s<S;s++)
	for(k=0;k<=N/2;k++)
            wavefunc[s][k] = exp(-0.5*pow(scale[s]*freq[k]-param,2.))*b;
    }
    
//    _mm256_exp_pd
    
    for(n=0;n<N;n++) datacomplex[n].setrealimga(data[n],0);
    
    dft_func<Type>(datacomplex,dft,N,1,pi,1,init);

    for(s=0;s<S;s++)
#if !defined(AVX) || AVX == 0    
        for(k=0;k<N;k++) dft_product[s][k] = dft[k]*wavefunc[s][k];
#elif AVX512 == 0
        if constexpr(sizeof(Type) == 8) {
	    for(k=0;k<N;k+=2) _mm256_store_pd((double *)&dft_product[s][k],_mm256_mul_pd(_mm256_load_pd((double *)&dft[k]),_mm256_setr_pd(wavefunc[s][k],wavefunc[s][k],wavefunc[s][k+1],wavefunc[s][k+1])));
	} else if constexpr(sizeof(Type) == 4) {
	    for(k=0;k<N;k+=4) _mm256_store_ps((float *)&dft_product[s][k],_mm256_mul_ps(_mm256_load_ps((float *)&dft[k]),_mm256_setr_ps(wavefunc[s][k],wavefunc[s][k],wavefunc[s][k+1],wavefunc[s][k+1],wavefunc[s][k+2],wavefunc[s][k+2],wavefunc[s][k+3],wavefunc[s][k+3])));
        }
#else
        if constexpr(sizeof(Type) == 8) {
	    for(k=0;k<N;k+=4) _mm512_store_pd((double *)&dft_product[s][k],_mm512_mul_pd(_mm512_load_pd((double *)&dft[k]),_mm512_setr_pd(wavefunc[s][k],wavefunc[s][k],wavefunc[s][k+1],wavefunc[s][k+1],wavefunc[s][k+2],wavefunc[s][k+2],wavefunc[s][k+3],wavefunc[s][k+3])));
	} else if constexpr(sizeof(Type) == 4) {
	    for(k=0;k<N;k+=8) _mm512_store_ps((float *)&dft_product[s][k],_mm512_mul_ps(_mm512_load_ps((float *)&dft[k]),_mm512_setr_ps(wavefunc[s][k],wavefunc[s][k],wavefunc[s][k+1],wavefunc[s][k+1],wavefunc[s][k+2],wavefunc[s][k+2],wavefunc[s][k+3],wavefunc[s][k+3],wavefunc[s][k+4],wavefunc[s][k+4],wavefunc[s][k+5],wavefunc[s][k+5],wavefunc[s][k+6],wavefunc[s][k+6],wavefunc[s][k+7],wavefunc[s][k+7])));
        }
#endif
    
    for(s=0;s<S;s++) dftinv_func<Type>(dft_product[s],transform[s],N,pi,0);

    for(s=0;s<S;s++) {
        a = sqrt(2.*pi*scale[s]/dx);
#if !defined(AVX) || AVX == 0
        for(n=0;n<N;n++) transform[s][n] *= a;
#elif AVX512 == 0
        if constexpr(sizeof(Type) == 8) {
    	    mw01_morlet_a = _mm256_set1_pd(a);
	    for(n=0;n<N;n+=2) _mm256_store_pd((double *)&transform[s][n],_mm256_mul_pd(_mm256_load_pd((double *)&transform[s][n]),mw01_morlet_a));
	} else if constexpr(sizeof(Type) == 4) {
       	    mw01_morlet_af = _mm256_set1_ps(a);
	    for(n=0;n<N;n+=4) _mm256_store_ps((float *)&transform[s][n],_mm256_mul_ps(_mm256_load_ps((float *)&transform[s][n]),mw01_morlet_af));
	}
#else
        if constexpr(sizeof(Type) == 8) {
    	    mw01_morlet_a = _mm512_set1_pd(a);
	    for(n=0;n<N;n+=4) _mm512_store_pd((double *)&transform[s][n],_mm512_mul_pd(_mm512_load_pd((double *)&transform[s][n]),mw01_morlet_a));
	} else if constexpr(sizeof(Type) == 4) {
       	    mw01_morlet_af = _mm512_set1_ps(a);
	    for(n=0;n<N;n+=8) _mm512_store_ps((float *)&transform[s][n],_mm512_mul_ps(_mm512_load_ps((float *)&transform[s][n]),mw01_morlet_af));
	}
#endif
    }
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
    Type thread_local pi = atan(1.)*4.;

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
    Type thread_local pi = atan(1.)*4.;

    intervals = (int)(upperbnd/dx);
    sum = 0.;

    for(i=1;i<=intervals;i++) {
        x = i*dx;
        sum += exp(-(x-param)*(x-param))/x;
    }
    sum *= pow(pi,-0.5)*dx;
    return sum;
}

