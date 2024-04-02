#include "complex.h"
#ifndef fft_bit_reverse
    #define fft_bit_reverse 0                                // in-place (1) or out-of-place (0) fft.
#endif
#ifndef MAXN
    #define MAXN 131072                                      // defined also in table.cpp. use 2^. arrays aligned to ALIGN  
#endif
#ifndef ALIGN
    #define ALIGN 64                                         // Bytes
#endif 
#define RaderMin 45                                          // if prime factor larger than this, use Rader algorithm
#if fft_bit_reverse == 1
    #include "table.h"
#endif

int primeFactors(int,int *);
int aligned_int(int,int);

template <class Type>
void Rader(complex<Type> *,complex<Type> *,complex<Type> *,int);

template <class Type>
void fft_func(complex<Type> *,complex<Type> *,int,int,int);

template <class Type>
void dft_func(complex<Type> *data,complex<Type> *out,int N,int Product,int sign) {
    Type thread_local pi = acos(-1.);      
    int thread_local oldN = 0;
    int thread_local i,j,k,m,n,p,q,r,t,inc;
    alignas(ALIGN) complex<Type> thread_local datasub1[MAXN];
    alignas(ALIGN) complex<Type> thread_local datasub2[MAXN];
    alignas(ALIGN) complex<Type> thread_local datasub3[MAXN];
    alignas(ALIGN) complex<Type> thread_local datasub4[MAXN];
    Type thread_local a;
    alignas(ALIGN) complex<Type> thread_local c[2];
    int thread_local PF,NoverPF;
    int thread_local kleft,kright,mleft,nright,tail;
    alignas(ALIGN) complex<Type> thread_local roots[MAXN];
    alignas(ALIGN) int thread_local Factor[100];
    int thread_local NFactor;
    alignas(ALIGN) complex<Type> thread_local datatemp[MAXN];
    int thread_local kkleft,kkright,mmleft,nnright;
    complex<Type> thread_local *dataptr,*outptr;
#if AVX512F > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m512d,__m512>::type avxtype;
    alignas(ALIGN) avxtype mw01_dft_a;
    alignas(ALIGN) avxtype mw01_dft_b;
    alignas(ALIGN) avxtype mw01_dft_c;
    alignas(ALIGN) avxtype mw01_dft_d;
    alignas(ALIGN) avxtype mw01_dft_e;
#elif AVX > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m256d,__m256>::type avxtype;
    alignas(ALIGN) avxtype mw01_dft_a;
    alignas(ALIGN) avxtype mw01_dft_b;
    alignas(ALIGN) avxtype mw01_dft_c;
    alignas(ALIGN) avxtype mw01_dft_d;
    alignas(ALIGN) avxtype mw01_dft_e;
#endif
    
    i = N & (N - 1);                                                               //  if 2^, use fft
    if(i == 0 && N >=4) {
        fft_func(data,out,N,Product,sign);
	return;
    }


    if(N != oldN) {
        if(sign > 0) a = -2.*pi/N; else a = 2.*pi/N;
        if(N%4 == 0) {
            k = N>>3;
            j = N>>2;
//	    _mm256_cos_pd(__m256d v1);   SVML
            for(i=0;i<=k;i++) roots[i] = complex<Type>(cos(a*i),sin(a*i));                   //     1/8 values
            if(sign > 0) {
                for(i=k+1;i<j;i++) roots[i] = reverse(swapcomplex(roots[j-i]));     //     values remaining in quadrant
#if !defined(AVX) || AVX == 0
                for(i=j;i<N/2;i++) roots[i] = turnright(roots[i-j]);           //     copy to next quadrant
#elif AVX512F == 0
		inc = 32/sizeof(complex<Type>);
                for(i=j;i<aligned_int(j,inc);i++) roots[i] = turnright(roots[i-j]);
                if constexpr(std::is_same_v<double,Type>) {
		    for(i=aligned_int(j,inc);i<N/2;i+=inc) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_permute_pd(_mm256_load_pd((double *)&roots[i-j]),0b0101),_mm256_setr_pd(0.0,-0.0,0.0,-0.0))); 
                } else if constexpr(std::is_same_v<float,Type>) {
                    for(i=aligned_int(j,inc);i<N/2;i+=inc) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_permute_ps(_mm256_load_ps((float *)&roots[i-j]),0b10110001),_mm256_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
                }
#else
		inc = 64/sizeof(complex<Type>);
                for(i=j;i<aligned_int(j,inc);i++) roots[i] = turnright(roots[i-j]);
                if constexpr(std::is_same_v<double,Type>) {
		    for(i=aligned_int(j,inc);i<N/2;i+=inc) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_permute_pd(_mm512_load_pd((double *)&roots[i-j]),0b01010101),_mm512_setr_pd(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));  
                } else if constexpr(std::is_same_v<float,Type>) {
		    for(i=aligned_int(j,inc);i<N/2;i+=inc) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_permute_ps(_mm512_load_ps((float *)&roots[i-j]),0b10110001),_mm512_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0))); 
                }
#endif
            } else {
                for(i=k+1;i<j;i++) roots[i] = swapcomplex(roots[j-i]);
#if !defined(AVX) || AVX == 0
                for(i=j;i<N/2;i++) roots[i] = turnleft(roots[i-j]);
#elif AVX512F == 0
	        inc = 32/sizeof(complex<Type>);
                for(i=j;i<aligned_int(j,inc);i++) roots[i] = turnleft(roots[i-j]);
                if constexpr(std::is_same_v<double,Type>) {
		    for(i=aligned_int(j,inc);i<N/2;i+=inc) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_permute_pd(_mm256_load_pd((double *)&roots[i-j]),0b0101),_mm256_setr_pd(-0.0,0.0,-0.0,0.0))); 
                } else if constexpr(std::is_same_v<float,Type>) {
		    for(i=aligned_int(j,inc);i<N/2;i+=inc) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_permute_ps(_mm256_load_ps((float *)&roots[i-j]),0b10110001),_mm256_setr_ps(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0)));
                }
#else
    	        inc = 64/sizeof(complex<Type>);
                for(i=j;i<aligned_int(j,inc);i++) roots[i] = turnleft(roots[i-j]);
                if constexpr(std::is_same_v<double,Type>) {
		    for(i=aligned_int(j,inc);i<N/2;i+=inc) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_permute_pd(_mm512_load_pd((double *)&roots[i-j]),0b01010101),_mm512_setr_pd(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0)));
                } else if constexpr(std::is_same_v<float,Type>) {
		    for(i=aligned_int(j,inc);i<N/2;i+=inc) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_permute_ps(_mm512_load_ps((float *)&roots[i-j]),0b10110001),_mm512_setr_ps(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0)));
                }
#endif
            }
#if !defined(AVX) || AVX == 0
            for(i=N/2;i<N;i++) roots[i] = reverse(roots[i-N/2]);                   // copy to whole circle
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
                for(i=N/2;i<N;i+=2) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_load_pd((double *)&roots[i-N/2]),_mm256_setr_pd(-0.0,-0.0,-0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
		inc = 32/sizeof(complex<Type>);
	        for(i=N/2;i<aligned_int(N/2,inc);i++) roots[i] = reverse(roots[i-N/2]);
		for(i=aligned_int(N/2,inc);i<N;i+=inc) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_load_ps((float *)&roots[i-N/2]),_mm256_setr_ps(-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0)));
            }
#else
            inc = 64/sizeof(complex<Type>);
            for(i=N/2;i<aligned_int(N/2,inc);i++) roots[i] = reverse(roots[i-N/2]);
            if constexpr(std::is_same_v<double,Type>) {
		for(i=aligned_int(N/2,inc);i<N;i+=inc) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_load_pd((double *)&roots[i-N/2]),_mm512_setr_pd(-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
		for(i=aligned_int(N/2,inc);i<N;i+=inc) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_load_ps((float *)&roots[i-N/2]),_mm512_setr_ps(-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0))); 
            }
#endif
        } else if(N%2 == 0) {
	    j = N>>2;
            for(i=0;i<=j;i++) roots[i] = complex<Type>(cos(a*i),sin(a*i));        // quadrant values
	    for(i=j+1;i<N/2;i++) roots[i] = realconj(roots[N/2-i]);               // copy to next quadrant
#if !defined(AVX) || AVX == 0
            for(i=N/2;i<N;i++) roots[i] = reverse(roots[i-N/2]);                       
#elif AVX512F == 0
            inc = 32/sizeof(complex<Type>);
            for(i=N/2;i<aligned_int(N/2,inc);i++) roots[i] = reverse(roots[i-N/2]);
            if constexpr(std::is_same_v<double,Type>) {    // N/2 must be odd
                for(i=aligned_int(N/2,inc);i<N;i+=inc) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_load_pd((double *)&roots[i-N/2]),_mm256_setr_pd(-0.0,-0.0,-0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
                for(i=aligned_int(N/2,inc);i<N;i+=inc) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_load_ps((float *)&roots[i-N/2]),_mm256_setr_ps(-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0))); 
            }
#else
            inc = 64/sizeof(complex<Type>);
            for(i=N/2;i<aligned_int(N/2,inc);i++) roots[i] = reverse(roots[i-N/2]);
            if constexpr(std::is_same_v<double,Type>) {    // N/2 must be odd
                for(i=aligned_int(N/2,inc);i<N;i+=inc) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_load_pd((double *)&roots[i-N/2]),_mm512_setr_pd(-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
                for(i=aligned_int(N/2,inc);i<N;i+=inc) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_load_ps((float *)&roots[i-N/2]),_mm512_setr_ps(-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0,-0.0))); 
            }
#endif
        } else {
            for(i=0;i<N/2+1;i++) roots[i] = complex<Type>(cos(a*i),sin(a*i));
            for(i=N/2+1;i<N;i++) roots[i] = conj(roots[N-i]);
        }
    } else {
        if((sign > 0 && imag(roots[1]) > 0.) || (sign < 0 && imag(roots[1]) < 0.)) {
#if !defined(AVX) || AVX == 0
            for(i=0;i<N;i++) roots[i] = conj(roots[i]);
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
                for(i=0;i<N;i+=2) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_load_pd((double *)&roots[i]),_mm256_setr_pd(0.0,-0.0,0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
                for(i=0;i<N;i+=4) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_load_ps((float *)&roots[i]),_mm256_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
                for(i=0;i<N;i+=4) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_load_pd((double *)&roots[i]),_mm512_setr_pd(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
                for(i=0;i<N;i+=8) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_load_ps((float *)&roots[i]),_mm512_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            }
#endif
        }
    }
   
 
    NFactor = primeFactors(N,Factor);
    for(j=0;j<NFactor;j++) {
        if(j == 0) {
            dataptr = data;
	    if(NFactor%2 == 0) outptr = datatemp; else outptr = out;
        } else if(j == 1) {
            if(NFactor%2 == 0) {
                dataptr = datatemp;
                outptr = out;
            } else {
                dataptr = out;
                outptr = datatemp;
            }
        }
   
        PF = Product*Factor[j];
        NoverPF = N/PF;
        kleft = NoverPF;
        kright = N/Product; 
        mleft = N/Factor[j];                                            
        nright = NoverPF;
        tail = NoverPF;

        if(Factor[j] == 2) {
            p = -NoverPF;
	    kkleft = -kleft;
	    kkright = -kright;
            for(k=0;k<Product;k++) {
                p += NoverPF;
	        kkleft += kleft;
	        kkright += kright;
                c[0] = roots[p];
#if !defined(AVX) || AVX == 0
                for(t=0;t<tail;t++) {
                    datasub2[t] = c[0]*dataptr[kkright+nright+t];
                    outptr[kkleft+t] = dataptr[kkright+t] + datasub2[t];
		    outptr[kkleft+mleft+t] = dataptr[kkright+t] - datasub2[t]; 
                }
#elif AVX512F == 0
                if constexpr(std::is_same_v<double,Type>) {
                    i = (kkright+nright)%2;
                    n = kkleft%2;
                    q = (kkleft+mleft)%2;
		    mw01_dft_c = _mm256_setr_pd(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                    for(t=0;t<tail-1;t+=2) {     // if tail=6,t=0,2,4,exit at 6; if tail=7,t=0,2,4,exit at 6
		        if(i == 0)
			    mw01_dft_b = complex_mul_256register(mw01_dft_c,_mm256_load_pd((double *)&dataptr[kkright+nright+t]));
			else
			    mw01_dft_b = complex_mul_256register(mw01_dft_c,_mm256_loadu_pd((double *)&dataptr[kkright+nright+t]));
                        //mw01_dft_b = complex_mul_256register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(double *)&dataptr[kkright+nright+t],i);              // nright=NoverPF
                        mw01_dft_a = _mm256_load_pd((double *)&dataptr[kkright+t]);  // kkright = k*N/P = even
                        if(n == 0)
                            _mm256_store_pd((double *)&outptr[kkleft+t],_mm256_add_pd(mw01_dft_a,mw01_dft_b));
                        else
                            _mm256_storeu_pd((double *)&outptr[kkleft+t],_mm256_add_pd(mw01_dft_a,mw01_dft_b));
                        if(q == 0)
                            _mm256_store_pd((double *)&outptr[kkleft+mleft+t],_mm256_sub_pd(mw01_dft_a,mw01_dft_b));
                        else
                            _mm256_storeu_pd((double *)&outptr[kkleft+mleft+t],_mm256_sub_pd(mw01_dft_a,mw01_dft_b));
                    }
                } else if constexpr(std::is_same_v<float,Type>) {
                    i = (kkright+nright)%4;
                    m = kkright%4;                 
                    n = kkleft%4;
                    q = (kkleft+mleft)%4;
    		    mw01_dft_c = _mm256_setr_ps(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                    for(t=0;t<tail-3;t+=4) {       // tail=1,2,3,exit at 0;  tail=4,5,6,7,exit at 4
                        if(i == 0)
                            mw01_dft_b = complex_mul_256register(mw01_dft_c,_mm256_load_ps((float *)&dataptr[kkright+nright+t]));
			else
                            mw01_dft_b = complex_mul_256register(mw01_dft_c,_mm256_loadu_ps((float *)&dataptr[kkright+nright+t]));
                        //mw01_dft_b = complex_mul_256register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(float *)&dataptr[kkright+nright+t],i);
                        if(m == 0)
                            mw01_dft_a = _mm256_load_ps((float *)&dataptr[kkright+t]);
                        else
                            mw01_dft_a = _mm256_loadu_ps((float *)&dataptr[kkright+t]); 
                        if(n == 0)
                            _mm256_store_ps((float *)&outptr[kkleft+t],_mm256_add_ps(mw01_dft_a,mw01_dft_b));
                        else
                            _mm256_storeu_ps((float *)&outptr[kkleft+t],_mm256_add_ps(mw01_dft_a,mw01_dft_b));
                        if(q == 0)
                            _mm256_store_ps((float *)&outptr[kkleft+mleft+t],_mm256_sub_ps(mw01_dft_a,mw01_dft_b));
                        else
                            _mm256_storeu_ps((float *)&outptr[kkleft+mleft+t],_mm256_sub_ps(mw01_dft_a,mw01_dft_b));
                    }
                }
#else
                if constexpr(std::is_same_v<double,Type>) {
                    i = (kkright+nright)%4;
                    m = kkright%4;
                    n = kkleft%4;
                    q = (kkleft+mleft)%4;
    		    mw01_dft_c = _mm512_setr_pd(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                    for(t=0;t<tail-3;t+=4) {     
		        if(i == 0)
			    mw01_dft_b = complex_mul_512register(mw01_dft_c,_mm512_load_pd((double *)&dataptr[kkright+nright+t]));
			else
			    mw01_dft_b = complex_mul_512register(mw01_dft_c,_mm512_loadu_pd((double *)&dataptr[kkright+nright+t]));
                        //mw01_dft_b = complex_mul_512register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(double *)&dataptr[kkright+nright+t],i);   
                        if(m == 0)
                            mw01_dft_a = _mm512_load_pd((double *)&dataptr[kkright+t]);
                        else
                            mw01_dft_a = _mm512_loadu_pd((double *)&dataptr[kkright+t]);
                        if(n == 0)
                            _mm512_store_pd((double *)&outptr[kkleft+t],_mm512_add_pd(mw01_dft_a,mw01_dft_b));
                        else
                            _mm512_storeu_pd((double *)&outptr[kkleft+t],_mm512_add_pd(mw01_dft_a,mw01_dft_b));
                        if(q == 0)
                            _mm512_store_pd((double *)&outptr[kkleft+mleft+t],_mm512_sub_pd(mw01_dft_a,mw01_dft_b));
                        else
                            _mm512_storeu_pd((double *)&outptr[kkleft+mleft+t],_mm512_sub_pd(mw01_dft_a,mw01_dft_b));
                    }
                } else if constexpr(std::is_same_v<float,Type>) {
                    i = (kkright+nright)%8;
                    m = kkright%8;                 
                    n = kkleft%8;
                    q = (kkleft+mleft)%8;
       		    mw01_dft_c = _mm512_setr_ps(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                    for(t=0;t<tail-7;t+=8) {   
		        if(i == 0)
			    mw01_dft_b = complex_mul_512register(mw01_dft_c,_mm512_load_ps((float *)&dataptr[kkright+nright+t]));
			else
 		            mw01_dft_b = complex_mul_512register(mw01_dft_c,_mm512_loadu_ps((float *)&dataptr[kkright+nright+t]));
                        //mw01_dft_b = complex_mul_512register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(float *)&dataptr[kkright+nright+t],i);
                        if(m == 0)
                            mw01_dft_a = _mm512_load_ps((float *)&dataptr[kkright+t]);
                        else
                            mw01_dft_a = _mm512_loadu_ps((float *)&dataptr[kkright+t]); 
                        if(n == 0)
                            _mm512_store_ps((float *)&outptr[kkleft+t],_mm512_add_ps(mw01_dft_a,mw01_dft_b));
                        else
                            _mm512_storeu_ps((float *)&outptr[kkleft+t],_mm512_add_ps(mw01_dft_a,mw01_dft_b));
                        if(q == 0)
                            _mm512_store_ps((float *)&outptr[kkleft+mleft+t],_mm512_sub_ps(mw01_dft_a,mw01_dft_b));
                        else
                            _mm512_storeu_ps((float *)&outptr[kkleft+mleft+t],_mm512_sub_ps(mw01_dft_a,mw01_dft_b));
                    }
                }
#endif
                for(;t<tail;t++) {
                    datasub2[t] = c[0]*dataptr[kkright+nright+t];
                    outptr[kkleft+t] = dataptr[kkright+t] + datasub2[t];
		    outptr[kkleft+mleft+t] = dataptr[kkright+t] - datasub2[t]; 
                }
	    }
        } else if(Factor[j] <= RaderMin) {
#if !defined(AVX) || AVX == 0
            memset(outptr,0,N*sizeof(complex<Type>));
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
                mw01_dft_a = _mm256_setzero_pd();
                for(i=0;i<N;i+=2)
                    _mm256_store_pd((double *)&outptr[i],mw01_dft_a);
            } else if constexpr(std::is_same_v<float,Type>) {
                mw01_dft_a = _mm256_setzero_ps();
                for(i=0;i<N;i+=4)
                    _mm256_store_ps((float *)&outptr[i],mw01_dft_a);
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
                mw01_dft_a = _mm512_setzero_pd();
                for(i=0;i<N;i+=4)
                    _mm512_store_pd((double *)&outptr[i],mw01_dft_a);
            } else if constexpr(std::is_same_v<float,Type>) {
                mw01_dft_a = _mm512_setzero_ps();
                for(i=0;i<N;i+=8)
                    _mm512_store_ps((float *)&outptr[i],mw01_dft_a);
            }
#endif
            kkleft = -kleft;
	    kkright = -kright;
	    for(k=0;k<Product;k++) {
	        kkleft += kleft;
	        kkright += kright;
	        nnright = -nright;
	        for(n=0;n<Factor[j];n++) {                                               // summation index
	            nnright += nright;
		    i = n*k*NoverPF;
                    c[0] = roots[i];
#if !defined(AVX) || AVX == 0
                    for(t=0;t<tail;t++) {
		        datasub2[t] = c[0]*dataptr[kkright+nnright+t];
		        outptr[kkleft+0*mleft+t] += datasub2[t];                            // m = 0
		        for(m=1;m<(Factor[j]+1)/2;m++) {
                            datasub1[m] = roots[n*m*mleft%N];
		            datasub3[m] = datasub2[t]*real(datasub1[m]);
			    datasub4[m] = turnleft(datasub2[t])*imag(datasub1[m]);
			    outptr[kkleft+m*mleft+t] += datasub3[m]+datasub4[m];
			    outptr[kkleft+(Factor[j]-m)*mleft+t] += datasub3[m]-datasub4[m];
		        }
		    }
#elif AVX512F == 0
                    if constexpr(std::is_same_v<double,Type>) {
                        q = (kkright+nnright)%2;   // kkright = k*N/P    nnright = n*NoverPF
                        p = kkleft%2;              // kkleft = k*NoverPF 
			mw01_dft_e = _mm256_setr_pd(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                        for(t=0;t<tail-1;t+=2) {   // tail = NoverPF, must be odd
			    if(q == 0)
			        mw01_dft_c = complex_mul_256register(mw01_dft_e,_mm256_load_pd((double *)&dataptr[kkright+nnright+t]));
			    else
    			        mw01_dft_c = complex_mul_256register(mw01_dft_e,_mm256_loadu_pd((double *)&dataptr[kkright+nnright+t]));
                            //mw01_dft_c = complex_mul_256register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(double *)&dataptr[kkright+nnright+t],q);
                            _mm256_store_pd((double *)&datasub2[t],mw01_dft_c);
                            if(p == 0)
                                _mm256_store_pd((double *)&outptr[kkleft+t],_mm256_add_pd(_mm256_load_pd((double *)&outptr[kkleft+t]),mw01_dft_c));  
                            else
                                _mm256_storeu_pd((double *)&outptr[kkleft+t],_mm256_add_pd(_mm256_loadu_pd((double *)&outptr[kkleft+t]),mw01_dft_c));

    		            //mw01_dft_c = _mm256_load_pd((double *)&datasub2[t]);
                            mw01_dft_d = _mm256_xor_pd(_mm256_permute_pd(mw01_dft_c,0b0101),_mm256_setr_pd(-0.0,0.0,-0.0,0.0));

				
                            for(m=1;m<(Factor[j]+1)/2;m++) {
                                datasub1[m] = roots[n*m*mleft%N];
				/*
                                for(r=0;r<2;r++) {
                                    datasub3[m] = datasub2[t+r]*real(datasub1[m]);                        
                                    datasub4[m] = turnleft(datasub2[t+r])*imag(datasub1[m]);             
                                    outptr[kkleft+m*mleft+t+r] += datasub3[m]+datasub4[m];                         
                                    outptr[kkleft+(Factor[j]-m)*mleft+t+r] += datasub3[m]-datasub4[m];        
                                }
				*/
				
                                mw01_dft_a = _mm256_set1_pd(real(datasub1[m]));
                                mw01_dft_b = _mm256_set1_pd(imag(datasub1[m]));
				
                                if((kkleft+m*mleft)%2 == 0)
				    //mw01_dft_e = _mm256_load_pd((double *)&outptr[kkleft+m*mleft+t]);
				    //mw01_dft_f = _mm_fmadd_pd(c,a,e);
				    //mw01_dft_g = _mm_fmadd_pd(d,b,f);
				    //_mm256_store_pd((double *)&outptr[kkleft+m*mleft+t],g);                
				    _mm256_store_pd((double *)&outptr[kkleft+m*mleft+t],_mm256_fmadd_pd(mw01_dft_d,mw01_dft_b,_mm256_fmadd_pd(mw01_dft_c,mw01_dft_a,_mm256_load_pd((double *)&outptr[kkleft+m*mleft+t]))));
				else    
    				    _mm256_storeu_pd((double *)&outptr[kkleft+m*mleft+t],_mm256_fmadd_pd(mw01_dft_d,mw01_dft_b,_mm256_fmadd_pd(mw01_dft_c,mw01_dft_a,_mm256_loadu_pd((double *)&outptr[kkleft+m*mleft+t]))));
                          
			        if((kkleft+(Factor[j]-m)*mleft)%2 == 0)
				    _mm256_store_pd((double *)&outptr[kkleft+(Factor[j]-m)*mleft+t],_mm256_fnmadd_pd(mw01_dft_d,mw01_dft_b,_mm256_fmadd_pd(mw01_dft_c,mw01_dft_a,_mm256_load_pd((double *)&outptr[kkleft+(Factor[j]-m)*mleft+t]))));
				else
				    _mm256_storeu_pd((double *)&outptr[kkleft+(Factor[j]-m)*mleft+t],_mm256_fnmadd_pd(mw01_dft_d,mw01_dft_b,_mm256_fmadd_pd(mw01_dft_c,mw01_dft_a,_mm256_loadu_pd((double *)&outptr[kkleft+(Factor[j]-m)*mleft+t]))));
			  
                            }
                        }
                    } else if constexpr(std::is_same_v<float,Type>) {
                        q = (kkright+nnright)%4;
                        p = kkleft%4;
			mw01_dft_e = _mm256_setr_ps(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                        for(t=0;t<tail-3;t+=4) { 
			    if(q == 0)
			        mw01_dft_c = complex_mul_256register(mw01_dft_e,_mm256_load_ps((float *)&dataptr[kkright+nnright+t]));
			    else
    			        mw01_dft_c = complex_mul_256register(mw01_dft_e,_mm256_loadu_ps((float *)&dataptr[kkright+nnright+t]));
                            //mw01_dft_c = complex_mul_256register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(float *)&dataptr[kkright+nnright+t],q);
                            _mm256_store_ps((float *)&datasub2[t],mw01_dft_c);
                            if(p == 0)
                                _mm256_store_ps((float *)&outptr[kkleft+t],_mm256_add_ps(_mm256_load_ps((float *)&outptr[kkleft+t]),mw01_dft_c));  
                            else
                                _mm256_storeu_ps((float *)&outptr[kkleft+t],_mm256_add_ps(_mm256_loadu_ps((float *)&outptr[kkleft+t]),mw01_dft_c));
				
    		            //mw01_dft_c = _mm256_load_ps((float *)&datasub2[t]);
                            mw01_dft_d = _mm256_xor_ps(_mm256_permute_ps(mw01_dft_c,0b10110001),_mm256_setr_ps(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0));
				
                            for(m=1;m<(Factor[j]+1)/2;m++) {
                                datasub1[m] = roots[n*m*mleft%N];
				/*
                                for(r=0;r<4;r++) {
                                    datasub3[m] = datasub2[t+r]*real(datasub1[m]);
                                    datasub4[m] = turnleft(datasub2[t+r])*imag(datasub1[m]);
                                    outptr[kkleft+m*mleft+t+r] += datasub3[m]+datasub4[m];
                                    outptr[kkleft+(Factor[j]-m)*mleft+t+r] += datasub3[m]-datasub4[m];
                                }
				*/
				
				mw01_dft_a = _mm256_set1_ps(real(datasub1[m]));
                                mw01_dft_b = _mm256_set1_ps(imag(datasub1[m]));
                                if((kkleft+m*mleft)%4 == 0)
				    _mm256_store_ps((float *)&outptr[kkleft+m*mleft+t],_mm256_fmadd_ps(mw01_dft_d,mw01_dft_b,_mm256_fmadd_ps(mw01_dft_c,mw01_dft_a,_mm256_load_ps((float *)&outptr[kkleft+m*mleft+t]))));
				else    
    				    _mm256_storeu_ps((float *)&outptr[kkleft+m*mleft+t],_mm256_fmadd_ps(mw01_dft_d,mw01_dft_b,_mm256_fmadd_ps(mw01_dft_c,mw01_dft_a,_mm256_loadu_ps((float *)&outptr[kkleft+m*mleft+t]))));
                          
			        if((kkleft+(Factor[j]-m)*mleft)%4 == 0)
				    _mm256_store_ps((float *)&outptr[kkleft+(Factor[j]-m)*mleft+t],_mm256_fnmadd_ps(mw01_dft_d,mw01_dft_b,_mm256_fmadd_ps(mw01_dft_c,mw01_dft_a,_mm256_load_ps((float *)&outptr[kkleft+(Factor[j]-m)*mleft+t]))));
				else
				    _mm256_storeu_ps((float *)&outptr[kkleft+(Factor[j]-m)*mleft+t],_mm256_fnmadd_ps(mw01_dft_d,mw01_dft_b,_mm256_fmadd_ps(mw01_dft_c,mw01_dft_a,_mm256_loadu_ps((float *)&outptr[kkleft+(Factor[j]-m)*mleft+t]))));
				
                            }
                        }
                    }
#else
                    if constexpr(std::is_same_v<double,Type>) {
                        q = (kkright+nnright)%4;   // kkright = k*N/P    nnright = n*NoverPF
                        p = kkleft%4;              // kkleft = k*NoverPF 
			mw01_dft_e = _mm512_setr_pd(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                        for(t=0;t<tail-3;t+=4) {   // tail = NoverPF, must be odd
			    if(q == 0)
			        mw01_dft_c = complex_mul_512register(mw01_dft_e,_mm512_load_pd((double *)&dataptr[kkright+nnright+t]));
			    else
    			        mw01_dft_c = complex_mul_512register(mw01_dft_e,_mm512_loadu_pd((double *)&dataptr[kkright+nnright+t]));
                            //mw01_dft_c = complex_mul_512register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(double *)&dataptr[kkright+nnright+t],q);
                            _mm512_store_pd((double *)&datasub2[t],mw01_dft_c);
                            if(p == 0)
                                _mm512_store_pd((double *)&outptr[kkleft+t],_mm512_add_pd(_mm512_load_pd((double *)&outptr[kkleft+t]),mw01_dft_c));  
                            else
                                _mm512_storeu_pd((double *)&outptr[kkleft+t],_mm512_add_pd(_mm512_loadu_pd((double *)&outptr[kkleft+t]),mw01_dft_c));
				
	  		    //mw01_dft_c = _mm512_load_pd((double *)&datasub2[t]);
                            mw01_dft_d = _mm512_xor_pd(_mm512_permute_pd(mw01_dft_c,0b01010101),_mm512_setr_pd(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0));
				
                            for(m=1;m<(Factor[j]+1)/2;m++) {
                                datasub1[m] = roots[n*m*mleft%N];
				/*
                                for(r=0;r<4;r++) {
                                    datasub3[m] = datasub2[t+r]*real(datasub1[m]);
                                    datasub4[m] = turnleft(datasub2[t+r])*imag(datasub1[m]);
                                    outptr[kkleft+m*mleft+t+r] += datasub3[m]+datasub4[m];
                                    outptr[kkleft+(Factor[j]-m)*mleft+t+r] += datasub3[m]-datasub4[m];
                                }
				*/
				
		                mw01_dft_a = _mm512_set1_pd(real(datasub1[m]));
                                mw01_dft_b = _mm512_set1_pd(imag(datasub1[m]));
                                if((kkleft+m*mleft)%4 == 0)
				    _mm512_store_pd((double *)&outptr[kkleft+m*mleft+t],_mm512_fmadd_pd(mw01_dft_d,mw01_dft_b,_mm512_fmadd_pd(mw01_dft_c,mw01_dft_a,_mm512_load_pd((double *)&outptr[kkleft+m*mleft+t]))));
				else    
    				    _mm512_storeu_pd((double *)&outptr[kkleft+m*mleft+t],_mm512_fmadd_pd(mw01_dft_d,mw01_dft_b,_mm512_fmadd_pd(mw01_dft_c,mw01_dft_a,_mm512_loadu_pd((double *)&outptr[kkleft+m*mleft+t]))));
                          
			        if((kkleft+(Factor[j]-m)*mleft)%4 == 0)
				    _mm512_store_pd((double *)&outptr[kkleft+(Factor[j]-m)*mleft+t],_mm512_fnmadd_pd(mw01_dft_d,mw01_dft_b,_mm512_fmadd_pd(mw01_dft_c,mw01_dft_a,_mm512_load_pd((double *)&outptr[kkleft+(Factor[j]-m)*mleft+t]))));
				else
				    _mm512_storeu_pd((double *)&outptr[kkleft+(Factor[j]-m)*mleft+t],_mm512_fnmadd_pd(mw01_dft_d,mw01_dft_b,_mm512_fmadd_pd(mw01_dft_c,mw01_dft_a,_mm512_loadu_pd((double *)&outptr[kkleft+(Factor[j]-m)*mleft+t]))));
				
                            }
                        }
                    } else if constexpr(std::is_same_v<float,Type>) {
                        q = (kkright+nnright)%8;
                        p = kkleft%8;
			mw01_dft_e = _mm512_setr_ps(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                        for(t=0;t<tail-7;t+=8) { 
			    if(q == 0)
			        mw01_dft_c = complex_mul_512register(mw01_dft_e,_mm512_load_ps((float *)&dataptr[kkright+nnright+t]));
			    else
    			        mw01_dft_c = complex_mul_512register(mw01_dft_e,_mm512_loadu_ps((float *)&dataptr[kkright+nnright+t]));
                            //mw01_dft_c = complex_mul_512register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(float *)&dataptr[kkright+nnright+t],q);
                            _mm512_store_ps((float *)&datasub2[t],mw01_dft_c);
                            if(p == 0)
                                _mm512_store_ps((float *)&outptr[kkleft+t],_mm512_add_ps(_mm512_load_ps((float *)&outptr[kkleft+t]),mw01_dft_c));  
                            else
                                _mm512_storeu_ps((float *)&outptr[kkleft+t],_mm512_add_ps(_mm512_loadu_ps((float *)&outptr[kkleft+t]),mw01_dft_c));
				
    		            //mw01_dft_c = _mm512_load_ps((float *)&datasub2[t]);
                            mw01_dft_d = _mm512_xor_ps(_mm512_permute_ps(mw01_dft_c,0b10110001),_mm512_setr_ps(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0));
				
                            for(m=1;m<(Factor[j]+1)/2;m++) {
                                datasub1[m] = roots[n*m*mleft%N];
				/*
                                for(r=0;r<8;r++) {
                                    datasub3[m] = datasub2[t+r]*real(datasub1[m]);
                                    datasub4[m] = turnleft(datasub2[t+r])*imag(datasub1[m]);
                                    outptr[kkleft+m*mleft+t+r] += datasub3[m]+datasub4[m];
                                    outptr[kkleft+(Factor[j]-m)*mleft+t+r] += datasub3[m]-datasub4[m];
                                }
				*/
				
				
    			        mw01_dft_a = _mm512_set1_ps(real(datasub1[m]));
                                mw01_dft_b = _mm512_set1_ps(imag(datasub1[m]));
                                if((kkleft+m*mleft)%8 == 0)
				    _mm512_store_ps((float *)&outptr[kkleft+m*mleft+t],_mm512_fmadd_ps(mw01_dft_d,mw01_dft_b,_mm512_fmadd_ps(mw01_dft_c,mw01_dft_a,_mm512_load_ps((float *)&outptr[kkleft+m*mleft+t]))));
				else    
    				    _mm512_storeu_ps((float *)&outptr[kkleft+m*mleft+t],_mm512_fmadd_ps(mw01_dft_d,mw01_dft_b,_mm512_fmadd_ps(mw01_dft_c,mw01_dft_a,_mm512_loadu_ps((float *)&outptr[kkleft+m*mleft+t]))));
                          
			        if((kkleft+(Factor[j]-m)*mleft)%8 == 0)
				    _mm512_store_ps((float *)&outptr[kkleft+(Factor[j]-m)*mleft+t],_mm512_fnmadd_ps(mw01_dft_d,mw01_dft_b,_mm512_fmadd_ps(mw01_dft_c,mw01_dft_a,_mm512_load_ps((float *)&outptr[kkleft+(Factor[j]-m)*mleft+t]))));
				else
				    _mm512_storeu_ps((float *)&outptr[kkleft+(Factor[j]-m)*mleft+t],_mm512_fnmadd_ps(mw01_dft_d,mw01_dft_b,_mm512_fmadd_ps(mw01_dft_c,mw01_dft_a,_mm512_loadu_ps((float *)&outptr[kkleft+(Factor[j]-m)*mleft+t]))));
				
                            }
                        }
                    }
#endif
                    for(;t<tail;t++) {
		        datasub2[t] = c[0]*dataptr[kkright+nnright+t];
		        outptr[kkleft+0*mleft+t] += datasub2[t];                            // m = 0
		        for(m=1;m<(Factor[j]+1)/2;m++) {
		            datasub1[m] = roots[n*m*mleft%N];
		            datasub3[m] = datasub2[t]*real(datasub1[m]);
		            datasub4[m] = turnleft(datasub2[t])*imag(datasub1[m]);
		            outptr[kkleft+m*mleft+t] += datasub3[m]+datasub4[m];
			    outptr[kkleft+(Factor[j]-m)*mleft+t] += datasub3[m]-datasub4[m];
		        }
		    }
	        } 
	    }
        } else {                                            // Rader
#if !defined(AVX) || AVX == 0
            memset(outptr,0,N*sizeof(complex<Type>));
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
                mw01_dft_a = _mm256_setzero_pd();
                for(i=0;i<N;i+=2)
                    _mm256_store_pd((double *)&outptr[i],mw01_dft_a);
            } else if constexpr(std::is_same_v<float,Type>) {
                mw01_dft_a = _mm256_setzero_ps();
                for(i=0;i<N;i+=4)
                    _mm256_store_ps((float *)&outptr[i],mw01_dft_a);
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
                mw01_dft_a = _mm512_setzero_pd();
                for(i=0;i<N;i+=4)
                    _mm512_store_pd((double *)&outptr[i],mw01_dft_a);
            } else if constexpr(std::is_same_v<float,Type>) {
                mw01_dft_a = _mm512_setzero_ps();
                for(i=0;i<N;i+=8)
                    _mm512_store_ps((float *)&outptr[i],mw01_dft_a);
            }
#endif
            kkleft = -kleft;
	    kkright = -kright;
            for(k=0;k<Product;k++) {
	        kkleft += kleft;
	        kkright += kright;
	        nnright = -nright;
	        for(n=0;n<Factor[j];n++) {                                                       //    m=0, summation index
		    p = n*k*NoverPF;
		    nnright += nright;
                    c[0] = roots[p];
#if !defined(AVX) || AVX == 0
                    for(t=0;t<tail;t++) 
		        outptr[kkleft+0*mleft+t] += c[0]*dataptr[kkright+nnright+t];
#elif AVX512F == 0
                    if constexpr(std::is_same_v<double,Type>) {
                        // tail is NoverPF, or product of the remaining factors. must be odd.
                        i = (kkright+nnright)%2;   // kkright = k*N/P   nnright = n*NoverPF
                        m = kkleft%2;              // kkleft = k*NoverPF
			mw01_dft_c = _mm256_setr_pd(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                        for(t=0;t<tail-1;t+=2) {
			    if(i == 0)
			        mw01_dft_b = complex_mul_256register(mw01_dft_c,_mm256_load_pd((double *)&dataptr[kkright+nnright+t]));
			    else
    			        mw01_dft_b = complex_mul_256register(mw01_dft_c,_mm256_loadu_pd((double *)&dataptr[kkright+nnright+t]));
                            //mw01_dft_b = complex_mul_256register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(double *)&dataptr[kkright+nnright+t],i);
                            if(m == 0)
                                mw01_dft_a = _mm256_load_pd((double *)&dataptr[kkleft+t]);
                            else
                                mw01_dft_a = _mm256_loadu_pd((double *)&dataptr[kkleft+t]);
                            if(m == 0)
                                _mm256_store_pd((double *)&outptr[kkleft+t],_mm256_add_pd(mw01_dft_a,mw01_dft_b));
                            else
                                _mm256_storeu_pd((double *)&outptr[kkleft+t],_mm256_add_pd(mw01_dft_a,mw01_dft_b));
                        }
                    } else if constexpr(std::is_same_v<float,Type>) {
                        i = (kkright+nnright)%4;
                        m = kkleft%4;
			mw01_dft_c = _mm256_setr_ps(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                        for(t=0;t<tail-3;t+=4) {
			    if(i == 0)
			        mw01_dft_b = complex_mul_256register(mw01_dft_c,_mm256_load_ps((float *)&dataptr[kkright+nnright+t]));
			    else
    			        mw01_dft_b = complex_mul_256register(mw01_dft_c,_mm256_loadu_ps((float *)&dataptr[kkright+nnright+t]));
                            //mw01_dft_b = complex_mul_256register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(float *)&dataptr[kkright+nnright+t],i);
                            if(m == 0)
                                mw01_dft_a = _mm256_load_ps((float *)&dataptr[kkleft+t]);
                            else
                                mw01_dft_a = _mm256_loadu_ps((float *)&dataptr[kkleft+t]);
                            if(m == 0)
                                _mm256_store_ps((float *)&outptr[kkleft+t],_mm256_add_ps(mw01_dft_a,mw01_dft_b));
                            else
                                _mm256_storeu_ps((float *)&outptr[kkleft+t],_mm256_add_ps(mw01_dft_a,mw01_dft_b));
                        }
                    }
#else
                    if constexpr(std::is_same_v<double,Type>) {
                        i = (kkright+nnright)%4;   
                        m = kkleft%4;   
			mw01_dft_c = _mm512_setr_pd(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                        for(t=0;t<tail-3;t+=4) {
			    if(i == 0)
			        mw01_dft_b = complex_mul_512register(mw01_dft_c,_mm512_load_pd((double *)&dataptr[kkright+nnright+t]));
			    else
   			        mw01_dft_b = complex_mul_512register(mw01_dft_c,_mm512_loadu_pd((double *)&dataptr[kkright+nnright+t]));
                            //mw01_dft_b = complex_mul_512register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(double *)&dataptr[kkright+nnright+t],i);
                            if(m == 0)
                                mw01_dft_a = _mm512_load_pd((double *)&dataptr[kkleft+t]);
                            else
                                mw01_dft_a = _mm512_loadu_pd((double *)&dataptr[kkleft+t]);
                            if(m == 0)
                                _mm512_store_pd((double *)&outptr[kkleft+t],_mm512_add_pd(mw01_dft_a,mw01_dft_b));
                            else
                                _mm512_storeu_pd((double *)&outptr[kkleft+t],_mm512_add_pd(mw01_dft_a,mw01_dft_b));
                        }
                    } else if constexpr(std::is_same_v<float,Type>) {
                        i = (kkright+nnright)%8;
                        m = kkleft%8;
			mw01_dft_c = _mm512_setr_ps(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                        for(t=0;t<tail-7;t+=8) {
			    if(i == 0)
			        mw01_dft_b = complex_mul_512register(mw01_dft_c,_mm512_load_ps((float *)&dataptr[kkright+nnright+t]));
			    else
  			        mw01_dft_b = complex_mul_512register(mw01_dft_c,_mm512_loadu_ps((float *)&dataptr[kkright+nnright+t]));
                            //mw01_dft_b = complex_mul_512register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(float *)&dataptr[kkright+nnright+t],i);
                            if(m == 0)
                                mw01_dft_a = _mm512_load_ps((float *)&dataptr[kkleft+t]);
                            else
                                mw01_dft_a = _mm512_loadu_ps((float *)&dataptr[kkleft+t]);
                            if(m == 0)
                                _mm512_store_ps((float *)&outptr[kkleft+t],_mm512_add_ps(mw01_dft_a,mw01_dft_b));
                            else
                                _mm512_storeu_ps((float *)&outptr[kkleft+t],_mm512_add_ps(mw01_dft_a,mw01_dft_b));
                        }
                    }
#endif
                    for(;t<tail;t++) 
		        outptr[kkleft+0*mleft+t] += c[0]*dataptr[kkright+nnright+t];
	        }
	    
	        for(q=1;q<Factor[j];q++)
                    datasub2[q] = roots[q*Product*NoverPF];
	        for(t=0;t<tail;t++) {                                                         //    m=1,.....
		    for(q=1;q<Factor[j];q++)
                        datasub1[q] = roots[q*k*NoverPF]*dataptr[kkright+q*nright+t];
		    Rader<Type>(datasub1,datasub2,datasub3,Factor[j]);
	            for(m=1;m<Factor[j];m++) 
		        outptr[kkleft+m*mleft+t] = dataptr[kkright+t] + datasub3[m];	    
	        }
	    }
        }

        if(PF != N) {
            if(j > 0) std::swap(outptr,dataptr);
        } else {
	    if(sign > 0) {
                a = 1./N;
#if !defined(AVX) || AVX == 0
                for(n=0;n<N;n++) out[n] *= a;
#elif AVX512F == 0
                if constexpr(std::is_same_v<double,Type>) {
		    mw01_dft_a = _mm256_set1_pd(a);
		    for(n=0;n<N;n+=2) _mm256_store_pd((double *)&out[n],_mm256_mul_pd(_mm256_load_pd((double *)&out[n]),mw01_dft_a));
		} else if constexpr(std::is_same_v<float,Type>) {
		    mw01_dft_a = _mm256_set1_ps(a);
    		    for(n=0;n<N;n+=4) _mm256_store_ps((float *)&out[n],_mm256_mul_ps(_mm256_load_ps((float *)&out[n]),mw01_dft_a));
		}
#else
                if constexpr(std::is_same_v<double,Type>) {
		    mw01_dft_a = _mm512_set1_pd(a);
    		    for(n=0;n<N;n+=4) _mm512_store_pd((double *)&out[n],_mm512_mul_pd(_mm512_load_pd((double *)&out[n]),mw01_dft_a));
		} else if constexpr(std::is_same_v<float,Type>) {
		    mw01_dft_a = _mm512_set1_ps(a);
       		    for(n=0;n<N;n+=8) _mm512_store_ps((float *)&out[n],_mm512_mul_ps(_mm512_load_ps((float *)&out[n]),mw01_dft_a));
		}
#endif
	    }
            break;	
        }

        Product = PF;
    }
    oldN = N;
}

#if fft_bit_reverse != 1
// out-of-place fft, N>=4
template <class Type>
void fft_func(complex<Type> *data,complex<Type> *out,int N,int Product,int sign) {    
    Type thread_local pi = acos(-1.);
    int thread_local oldN = 0;
    int thread_local i,j,k,m,n,p,q,t;
    alignas(ALIGN) complex<Type> thread_local c[8];
    Type thread_local a;
    int thread_local PF,NoverPF;
    int thread_local kleft,kright,mleft,nright,tail;
    alignas(ALIGN) complex<Type> thread_local roots[MAXN];
    alignas(ALIGN) complex<Type> thread_local datatemp[MAXN];
    int thread_local kkleft,kkright,mmleft,nnright;
    complex<Type> thread_local *dataptr,*outptr;
#if AVX512F > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m512d,__m512>::type avxtype;
    alignas(ALIGN) avxtype mw01_fft0_a;
    alignas(ALIGN) avxtype mw01_fft0_b;
    alignas(ALIGN) avxtype mw01_fft0_c;
#elif AVX > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m256d,__m256>::type avxtype;
    alignas(ALIGN) avxtype mw01_fft0_a;
    alignas(ALIGN) avxtype mw01_fft0_b;
    alignas(ALIGN) avxtype mw01_fft0_c;
#endif


    if(N != oldN) {
        if(sign > 0) a = -2.*pi/N; else a = 2.*pi/N;
        k = N>>3;
        j = N>>2;
        for(i=0;i<=k;i++) roots[i] = complex<Type>(cos(a*i),sin(a*i));          //     1/8 values
        if(sign > 0) {
            for(i=k+1;i<j;i++) roots[i] = reverse(swapcomplex(roots[j-i]));          //     values remaining in quadrant
#if !defined(AVX) || AVX == 0 
            for(i=j;i<N/2;i++) roots[i] = turnright(roots[i-j]);
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
	        for(i=j;i<aligned_int(j,2);i++) roots[i] = turnright(roots[i-j]);
		for(i=aligned_int(j,2);i<N/2;i+=2) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_permute_pd(_mm256_load_pd((double *)&roots[i-j]),0b0101),_mm256_setr_pd(0.0,-0.0,0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
	        for(i=j;i<aligned_int(j,4);i++) roots[i] = turnright(roots[i-j]);
		for(i=aligned_int(j,4);i<N/2;i+=4) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_permute_ps(_mm256_load_ps((float *)&roots[i-j]),0b10110001),_mm256_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0))); 
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
	        for(i=j;i<aligned_int(j,4);i++) roots[i] = turnright(roots[i-j]);
		for(i=aligned_int(j,4);i<N/2;i+=4) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_permute_pd(_mm512_load_pd((double *)&roots[i-j]),0b01010101),_mm512_setr_pd(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
  	        for(i=j;i<aligned_int(j,8);i++) roots[i] = turnright(roots[i-j]);
		for(i=aligned_int(j,8);i<N/2;i+=8) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_permute_ps(_mm512_load_ps((float *)&roots[i-j]),0b10110001),_mm512_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0))); 
            }
#endif
        } else {
            for(i=k+1;i<j;i++) roots[i] = swapcomplex(roots[j-i]);
#if !defined(AVX) || AVX == 0
            for(i=j;i<N/2;i++) roots[i] = turnleft(roots[i-j]);
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
	        for(i=j;i<aligned_int(j,2);i++) roots[i] = turnleft(roots[i-j]);
		for(i=aligned_int(j,2);i<N/2;i+=2) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_permute_pd(_mm256_load_pd((double *)&roots[i-j]),0b0101),_mm256_setr_pd(-0.0,0.0,-0.0,0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
	        for(i=j;i<aligned_int(j,4);i++) roots[i] = turnleft(roots[i-j]);
		for(i=aligned_int(j,4);i<N/2;i+=4) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_permute_ps(_mm256_load_ps((float *)&roots[i-j]),0b10110001),_mm256_setr_ps(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0))); 
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
	        for(i=j;i<aligned_int(j,4);i++) roots[i] = turnleft(roots[i-j]);
		for(i=aligned_int(j,4);i<N/2;i+=4) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_permute_pd(_mm512_load_pd((double *)&roots[i-j]),0b01010101),_mm512_setr_pd(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
	        for(i=j;i<aligned_int(j,8);i++) roots[i] = turnleft(roots[i-j]);
		for(i=aligned_int(j,8);i<N/2;i+=8) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_permute_ps(_mm512_load_ps((float *)&roots[i-j]),0b10110001),_mm512_setr_ps(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0))); 
            }
#endif
        }
    } else {
        if((sign > 0 && imag(roots[1]) > 0.) || (sign < 0 && imag(roots[1]) < 0.)) {
#if !defined(AVX) || AVX == 0
            for(i=0;i<N;i++) roots[i] = conj(roots[i]);
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
                for(i=0;i<N;i+=2) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_load_pd((double *)&roots[i]),_mm256_setr_pd(0.0,-0.0,0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
                for(i=0;i<N;i+=4) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_load_ps((float *)&roots[i]),_mm256_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
                for(i=0;i<N;i+=4) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_load_pd((double *)&roots[i]),_mm512_setr_pd(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
                for(i=0;i<N;i+=8) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_load_ps((float *)&roots[i]),_mm512_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            }
#endif
        }
    }
    
    i = N;                                                       
    j = 0;
    while(i != 1) { i>>=1; j++; }                                             // N = 2^j
    j = j%2;
   
 
    NoverPF = N;
    while(1) {
        if(Product == 1) {
            dataptr = data;
            if(j == 0) outptr = datatemp; else outptr = out;
        } else if(Product == 2) {
            if(j == 0) {
                dataptr = datatemp;
                outptr = out;
            } else {
                dataptr = out;
                outptr = datatemp;
            }
        }
  
        PF = Product<<1;
        NoverPF >>= 1;
        kleft = NoverPF;
        kright = NoverPF<<1; 
        mleft = N>>1;                                                 
        nright = NoverPF;
        tail = NoverPF;

        kkleft = -kleft;
        p = -NoverPF;
        for(k=0;k<Product;k++) {
            kkleft += kleft;
	    kkright = kkleft<<1;
            p += NoverPF;
            c[0] = roots[p];
#if !defined(AVX) || AVX == 0 
            for(t=0;t<tail;t++) {
                c[1] = c[0]*dataptr[kkright+nright+t];
                outptr[kkleft+t] = dataptr[kkright+t] + c[1];
	        outptr[kkleft+mleft+t] = dataptr[kkright+t] - c[1]; 
            }
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
	        mw01_fft0_c = _mm256_setr_pd(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                for(t=0;t<tail-1;t+=2) {  // tail = 1,2,4,.... = kleft = nright    mleft = N/2   kkright = 2*kkleft     
                    mw01_fft0_b = complex_mul_256register(mw01_fft0_c,_mm256_load_pd((double *)&dataptr[kkright+nright+t]));
		    //mw01_fft0_b = complex_mul_256register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(double *)&dataptr[kkright+nright+t],0);
                    mw01_fft0_a = _mm256_load_pd((double *)&dataptr[kkright+t]);
                    _mm256_store_pd((double *)&outptr[kkleft+t],_mm256_add_pd(mw01_fft0_a,mw01_fft0_b));
                    _mm256_store_pd((double *)&outptr[kkleft+mleft+t],_mm256_sub_pd(mw01_fft0_a,mw01_fft0_b));
                }
            } else if constexpr(std::is_same_v<float,Type>) {
	    	mw01_fft0_c = _mm256_setr_ps(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                for(t=0;t<tail-3;t+=4) {   
                    mw01_fft0_b = complex_mul_256register(mw01_fft0_c,_mm256_load_ps((float *)&dataptr[kkright+nright+t]));
                    //mw01_fft0_b = complex_mul_256register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(float *)&dataptr[kkright+nright+t],0);
                    mw01_fft0_a = _mm256_load_ps((float *)&dataptr[kkright+t]);
                    _mm256_store_ps((float *)&outptr[kkleft+t],_mm256_add_ps(mw01_fft0_a,mw01_fft0_b));
                    _mm256_store_ps((float *)&outptr[kkleft+mleft+t],_mm256_sub_ps(mw01_fft0_a,mw01_fft0_b));
                }
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
       	    	mw01_fft0_c = _mm512_setr_pd(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                for(t=0;t<tail-3;t+=4) {       
                    mw01_fft0_b = complex_mul_512register(mw01_fft0_c,_mm512_load_pd((double *)&dataptr[kkright+nright+t]));
                    //mw01_fft0_b = complex_mul_512register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(double *)&dataptr[kkright+nright+t],0);
                    mw01_fft0_a = _mm512_load_pd((double *)&dataptr[kkright+t]);
                    _mm512_store_pd((double *)&outptr[kkleft+t],_mm512_add_pd(mw01_fft0_a,mw01_fft0_b));
                    _mm512_store_pd((double *)&outptr[kkleft+mleft+t],_mm512_sub_pd(mw01_fft0_a,mw01_fft0_b));
                }
            } else if constexpr(std::is_same_v<float,Type>) {
       	    	mw01_fft0_c = _mm512_setr_ps(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]));
                for(t=0;t<tail-7;t+=8) {
		    mw01_fft0_b = complex_mul_512register(mw01_fft0_c,_mm512_load_ps((float *)&dataptr[kkright+nright+t]));
                    //mw01_fft0_b = complex_mul_512register(real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),real(c[0]),imag(c[0]),(float *)&dataptr[kkright+nright+t],0);
                    mw01_fft0_a = _mm512_load_ps((float *)&dataptr[kkright+t]);
                    _mm512_store_ps((float *)&outptr[kkleft+t],_mm512_add_ps(mw01_fft0_a,mw01_fft0_b));
                    _mm512_store_ps((float *)&outptr[kkleft+mleft+t],_mm512_sub_ps(mw01_fft0_a,mw01_fft0_b));
                }
            }
#endif
            for(;t<tail;t++) {
                c[1] = c[0]*dataptr[kkright+nright+t];
                outptr[kkleft+t] = dataptr[kkright+t] + c[1];
                outptr[kkleft+mleft+t] = dataptr[kkright+t] - c[1];
            }
	}

        if(PF != N) {
            if(Product > 1) std::swap(outptr,dataptr);
        } else {
	    if(sign > 0) {
                a = 1./N;
#if !defined(AVX) || AVX == 0
                for(n=0;n<N;n++) out[n] *= a;
#elif AVX512F == 0
                if constexpr(std::is_same_v<double,Type>) {
		    mw01_fft0_a = _mm256_set1_pd(a);
		    for(n=0;n<N;n+=2) _mm256_store_pd((double *)&out[n],_mm256_mul_pd(_mm256_load_pd((double *)&out[n]),mw01_fft0_a));
		} else if constexpr(std::is_same_v<float,Type>) {
		    mw01_fft0_a = _mm256_set1_ps(a);
    		    for(n=0;n<N;n+=4) _mm256_store_ps((float *)&out[n],_mm256_mul_ps(_mm256_load_ps((float *)&out[n]),mw01_fft0_a));
		}
#else
                if constexpr(std::is_same_v<double,Type>) {
		    mw01_fft0_a = _mm512_set1_pd(a);
    		    for(n=0;n<N;n+=4) _mm512_store_pd((double *)&out[n],_mm512_mul_pd(_mm512_load_pd((double *)&out[n]),mw01_fft0_a));
		} else if constexpr(std::is_same_v<float,Type>) {
		    mw01_fft0_a = _mm512_set1_ps(a);
       		    for(n=0;n<N;n+=8) _mm512_store_ps((float *)&out[n],_mm512_mul_ps(_mm512_load_ps((float *)&out[n]),mw01_fft0_a));
		}
#endif
	    }
            break;	
        }

        Product = PF;
    }
    oldN = N;
}
#endif

#if fft_bit_reverse == 1
// in-place fft, N>=4
template <class Type>
void fft_func(complex<Type> *data,complex<Type> *out,int N,int Product,int sign) {
    Type thread_local pi = acos(-1.);
    int thread_local oldN = 0;      
    int thread_local i,j,k,m,n,p,q,h;
    alignas(ALIGN) complex<Type> thread_local c[8];
    Type thread_local a,b;
    int thread_local PF,NoverPF;
    int thread_local head,hhead;
    alignas(ALIGN) complex<Type> thread_local roots[MAXN];
#if AVX512F > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m512d,__m512>::type avxtype;
    alignas(ALIGN) avxtype mw01_fft1_a;
    alignas(ALIGN) avxtype mw01_fft1_b;
#elif AVX > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m256d,__m256>::type avxtype;
    alignas(ALIGN) avxtype mw01_fft1_a;
    alignas(ALIGN) avxtype mw01_fft1_b;
#endif


    if(N != oldN) {
        if(sign > 0) a = -2.*pi/N; else a = 2.*pi/N;
        k = N>>3;
        j = N>>2;
        for(i=0;i<=k;i++) roots[i] = complex<Type>(cos(a*i),sin(a*i));           //     1/8 values
        if(sign > 0) {
            for(i=k+1;i<j;i++) roots[i] = reverse(swapcomplex(roots[j-i]));           //     values remaining in quadrant
#if !defined(AVX) || AVX == 0
            for(i=j;i<N/2;i++) roots[i] = turnright(roots[i-j]);
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
	        for(i=j;i<aligned_int(j,2);i++) roots[i] = turnright(roots[i-j]);
		for(i=aligned_int(j,2);i<N/2;i+=2) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_permute_pd(_mm256_load_pd((double *)&roots[i-j]),0b0101),_mm256_setr_pd(0.0,-0.0,0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
	        for(i=j;i<aligned_int(j,4);i++) roots[i] = turnright(roots[i-j]);
		for(i=aligned_int(j,4);i<N/2;i+=4) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_permute_ps(_mm256_load_ps((float *)&roots[i-j]),0b10110001),_mm256_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
	        for(i=j;i<aligned_int(j,4);i++) roots[i] = turnright(roots[i-j]);
		for(i=aligned_int(j,4);i<N/2;i+=4) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_permute_pd(_mm512_load_pd((double *)&roots[i-j]),0b01010101),_mm512_setr_pd(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0))); 
            } else if constexpr(std::is_same_v<float,Type>) {
	        for(i=j;i<aligned_int(j,8);i++) roots[i] = turnright(roots[i-j]);
		for(i=aligned_int(j,8);i<N/2;i+=8) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_permute_ps(_mm512_load_ps((float *)&roots[i-j]),0b10110001),_mm512_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            }
#endif
        } else {
            for(i=k+1;i<j;i++) roots[i] = swapcomplex(roots[j-i]);
#if !defined(AVX) || AVX == 0
            for(i=j;i<N/2;i++) roots[i] = turnleft(roots[i-j]);
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
	        for(i=j;i<aligned_int(j,2);i++) roots[i] = turnleft(roots[i-j]);
		for(i=aligned_int(j,2);i<N/2;i+=2) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_permute_pd(_mm256_load_pd((double *)&roots[i-j]),0b0101),_mm256_setr_pd(-0.0,0.0,-0.0,0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
	        for(i=j;i<aligned_int(j,4);i++) roots[i] = turnleft(roots[i-j]);
		for(i=aligned_int(j,4);i<N/2;i+=4) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_permute_ps(_mm256_load_ps((float *)&roots[i-j]),0b10110001),_mm256_setr_ps(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0)));
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
	        for(i=j;i<aligned_int(j,4);i++) roots[i] = turnleft(roots[i-j]);
		for(i=aligned_int(j,4);i<N/2;i+=4) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_permute_pd(_mm512_load_pd((double *)&roots[i-j]),0b01010101),_mm512_setr_pd(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
	        for(i=j;i<aligned_int(j,8);i++) roots[i] = turnleft(roots[i-j]);
		for(i=aligned_int(j,8);i<N/2;i+=8) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_permute_ps(_mm512_load_ps((float *)&roots[i-j]),0b10110001),_mm512_setr_ps(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0)));
            }
#endif
        }
    } else {
        if((sign > 0 && imag(roots[1]) > 0.) || (sign < 0 && imag(roots[1]) < 0.)) {
#if !defined(AVX) || AVX == 0
            for(i=0;i<N;i++) roots[i] = conj(roots[i]);
#elif AVX512F == 0
            if constexpr(std::is_same_v<double,Type>) {
                for(i=0;i<N;i+=2) _mm256_store_pd((double *)&roots[i],_mm256_xor_pd(_mm256_load_pd((double *)&roots[i]),_mm256_setr_pd(0.0,-0.0,0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
                for(i=0;i<N;i+=4) _mm256_store_ps((float *)&roots[i],_mm256_xor_ps(_mm256_load_ps((float *)&roots[i]),_mm256_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            }
#else
            if constexpr(std::is_same_v<double,Type>) {
                for(i=0;i<N;i+=4) _mm512_store_pd((double *)&roots[i],_mm512_xor_pd(_mm512_load_pd((double *)&roots[i]),_mm512_setr_pd(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            } else if constexpr(std::is_same_v<float,Type>) {
                for(i=0;i<N;i+=8) _mm512_store_ps((float *)&roots[i],_mm512_xor_ps(_mm512_load_ps((float *)&roots[i]),_mm512_setr_ps(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)));
            }
#endif
        }
    }


    i = N;
    j = 0;
    while(i != 1) { i>>=1; j++; }

    NoverPF = N;
    while(1) {
        PF = Product<<1;
        NoverPF >>= 1;
        head = NoverPF;
        hhead = -PF;

        for(h=0;h<head;h++) {
            hhead += PF;
            p = -NoverPF;
            if(Product == 1) {
                p += NoverPF;
                c[0] = data[table[j][hhead]];
                c[1] = roots[p]*data[table[j][hhead+Product]];
                out[hhead] = c[0] + c[1];
                out[hhead+Product] = c[0] - c[1];
            } else {
#if !defined(AVX) || AVX == 0
                for(k=0;k<Product;k++) {
                    p += NoverPF;
                    c[0] = out[hhead+k];
                    c[1] = roots[p]*out[hhead+Product+k];
                    out[hhead+k] = c[0] + c[1];
                    out[hhead+Product+k] = c[0] - c[1];
                }
#elif AVX512F == 0
                if constexpr(std::is_same_v<double,Type>) {
                    for(k=0;k<Product-1;k+=2) {
                        p += NoverPF; c[0] = roots[p];
                        p += NoverPF; c[1] = roots[p]; 
                        mw01_fft1_b = complex_mul_256register((double *)&c[0],(double *)&out[hhead+Product+k],0,0);
                        mw01_fft1_a = _mm256_load_pd((double *)&out[hhead+k]);  // hhead = h*PF 
                        _mm256_store_pd((double *)&out[hhead+k],_mm256_add_pd(mw01_fft1_a,mw01_fft1_b));
                        _mm256_store_pd((double *)&out[hhead+Product+k],_mm256_sub_pd(mw01_fft1_a,mw01_fft1_b));
                    }
                } else if constexpr(std::is_same_v<float,Type>) {
                    for(k=0;k<Product-3;k+=4) {
                        for(i=0;i<4;i++) { p += NoverPF; c[i] = roots[p]; }
                        mw01_fft1_b = complex_mul_256register((float *)&c[0],(float *)&out[hhead+Product+k],0,0);
                        mw01_fft1_a = _mm256_load_ps((float *)&out[hhead+k]);   // hhead = h*PF
                        _mm256_store_ps((float *)&out[hhead+k],_mm256_add_ps(mw01_fft1_a,mw01_fft1_b));
                        _mm256_store_ps((float *)&out[hhead+Product+k],_mm256_sub_ps(mw01_fft1_a,mw01_fft1_b));
                    }
                }
#else
                if constexpr(std::is_same_v<double,Type>) {
                    for(k=0;k<Product-3;k+=4) {
                        for(i=0;i<4;i++) { p += NoverPF; c[i] = roots[p]; }
                        mw01_fft1_b = complex_mul_512register((double *)&c[0],(double *)&out[hhead+Product+k],0,0);
                        mw01_fft1_a = _mm512_load_pd((double *)&out[hhead+k]);   
                        _mm512_store_pd((double *)&out[hhead+k],_mm512_add_pd(mw01_fft1_a,mw01_fft1_b));
                        _mm512_store_pd((double *)&out[hhead+Product+k],_mm512_sub_pd(mw01_fft1_a,mw01_fft1_b));
                    }
                } else if constexpr(std::is_same_v<float,Type>) {
                    for(k=0;k<Product-7;k+=8) {
                        for(i=0;i<8;i++) { p += NoverPF; c[i] = roots[p]; } 
                        mw01_fft1_b = complex_mul_512register((float *)&c[0],(float *)&out[hhead+Product+k],0,0);
                        mw01_fft1_a = _mm512_load_ps((float *)&out[hhead+k]);   
                        _mm512_store_ps((float *)&out[hhead+k],_mm512_add_ps(mw01_fft1_a,mw01_fft1_b));
                        _mm512_store_ps((float *)&out[hhead+Product+k],_mm512_sub_ps(mw01_fft1_a,mw01_fft1_b));
                    }
                }
#endif
                for(;k<Product;k++) {
                    p += NoverPF;
                    c[0] = out[hhead+k];
                    c[1] = roots[p]*out[hhead+Product+k];
                    out[hhead+k] = c[0] + c[1];
                    out[hhead+Product+k] = c[0] - c[1];
                }
            }
        }

        if(PF == N) {
            if(sign > 0) {
                a = 1./N;
#if !defined(AVX) || AVX == 0
                for(n=0;n<N;n++) out[n] *= a;
#elif AVX512F == 0
                if constexpr(std::is_same_v<double,Type>) {
		    mw01_fft1_a = _mm256_set1_pd(a);
		    for(n=0;n<N;n+=2) _mm256_store_pd((double *)&out[n],_mm256_mul_pd(_mm256_load_pd((double *)&out[n]),mw01_fft1_a));
		} else if constexpr(std::is_same_v<float,Type>) {
		    mw01_fft1_a = _mm256_set1_ps(a);
    		    for(n=0;n<N;n+=4) _mm256_store_ps((float *)&out[n],_mm256_mul_ps(_mm256_load_ps((float *)&out[n]),mw01_fft1_a));
		}
#else
                if constexpr(std::is_same_v<double,Type>) {
		    mw01_fft1_a = _mm512_set1_pd(a);
    		    for(n=0;n<N;n+=4) _mm512_store_pd((double *)&out[n],_mm512_mul_pd(_mm512_load_pd((double *)&out[n]),mw01_fft1_a));
		} else if constexpr(std::is_same_v<float,Type>) {
		    mw01_fft1_a = _mm512_set1_ps(a);
       		    for(n=0;n<N;n+=8) _mm512_store_ps((float *)&out[n],_mm512_mul_ps(_mm512_load_ps((float *)&out[n]),mw01_fft1_a));
		}
#endif
            }
            break;
        }

        Product = PF;
    }
    oldN = N;
}
#endif

template <class Type>
void dftinv_func(complex<Type> *data,complex<Type> *out,int N) {
    dft_func<Type>(data,out,N,1,-1);
}

template <class Type>
void fftinv_func(complex<Type> *data,complex<Type> *out,int N) {
    fft_func<Type>(data,out,N,1,-1);
}

int primeFactors(int N,int *f) {
    int nfactors = 0;
    while(N%2 == 0) {
        N = N>>1;
        f[nfactors] = 2;
        nfactors++;
    }
    for(int i=3;i<=sqrt(N);i+=2) {
        while(N%i == 0) {
            f[nfactors] = i;
            N = N/i;
            nfactors++;
        }
    }
    if(N > 2) {
        f[nfactors] = N;
        nfactors++;
    }
    return nfactors;
}

int powmod(int a,int b,int p) {
    int res = 1;
    while (b)
        if (b&1)
            res = (int)(res*1ll*a%p), --b;       // 1ll   is long*long
        else
            a = (int)(a*1ll*a%p), b >>= 1;
    return res;
}

int generator(int p) {
    alignas(ALIGN) int fact[100];
    int m=0;
    int phi=p-1, n=phi;
    for (int i=2;i*i<=n;++i)
        if (n % i == 0) {
            fact[m] = i;
            m++;
            while(n%i == 0) n/=i;
        }
    if (n > 1) { fact[m] = n; m++; }

    for (int res=2;res<=p;++res) {
        bool ok = true;
        for (int i=0;i<m && ok;++i) ok &= powmod(res,phi/fact[i],p) != 1;
        if (ok) return res;
    }
    return -1;
}

int gcd1(int a,int b,int &x,int &y) {
    if(b == 0) {
        x = 1;
	y = 0;
	return a;
    }
    int x1, y1, gcd = gcd1(b,a%b,x1,y1);
    x = y1;
    y = x1 - (a/b)*y1;
    return gcd;
}

// p is prime, g is generator
int modulo_inverse(int g,int p) {
    int x,y;
    int z = gcd1(g,p,x,y);
    x = (x%p + p)%p;
    return (int)x;
}

int aligned_int(int start,int increment) {
    return start+(increment-start%increment)%increment;
}

template <class Type>
void Rader(complex<Type> *datasub1,complex<Type> *datasub2,complex<Type> *out,int N) {
    Type thread_local pi = acos(-1.);
    int thread_local oldN = 0;
    int thread_local g,ginv;
    alignas(ALIGN) int thread_local mapg[MAXN];
    alignas(ALIGN) int thread_local mapginv[MAXN];
    int thread_local newN;
    Type thread_local newN2;
    int thread_local i,p,q;
    alignas(ALIGN) complex<Type> thread_local padded[MAXN];
    alignas(ALIGN) complex<Type> thread_local result1[MAXN];
    alignas(ALIGN) complex<Type> thread_local result2[MAXN];
#if AVX512F > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m512d,__m512>::type avxtype;
    alignas(ALIGN) avxtype mw01_rader_a;
#elif AVX > 0
    typedef typename std::conditional<std::is_same_v<double,Type>,__m256d,__m256>::type avxtype;
    alignas(ALIGN) avxtype mw01_rader_a;
#endif


    if(N != oldN) {
        g = generator(N);
        ginv = modulo_inverse(g,N);
        mapg[0] = 1;
        for(q=1;q<N-1;q++) mapg[q] = mapg[q-1]*g%N;                     // n = g^q    (mod N)  , q=[0,N-2]
        mapginv[0] = 1;
        for(p=1;p<N-1;p++) mapginv[p] = mapg[N-p-1];                    // k/m = ginv^p (mod N)  , p=[0,N-2]

        // padding to 2^
        newN = 2;
        for(i=0;i<1000000;i++) {
            if(newN >= 2*(N-1)-1) 
	        break;
            else 
                newN <<= 1;
        }
    }
    
#if !defined(AVX) || AVX == 0
    memset(padded,0,newN*sizeof(complex<Type>));
#elif AVX512F == 0
    if constexpr(std::is_same_v<double,Type>) {
        mw01_rader_a = _mm256_setzero_pd();
        for(i=0;i<newN;i+=2)
            _mm256_store_pd((double *)&padded[i],mw01_rader_a);
    } else if constexpr(std::is_same_v<float,Type>) {
        mw01_rader_a = _mm256_setzero_ps();
        for(i=0;i<newN;i+=4)
            _mm256_store_ps((float *)&padded[i],mw01_rader_a);
    }
#else
    if constexpr(std::is_same_v<double,Type>) {
        mw01_rader_a = _mm512_setzero_pd();
        for(i=0;i<newN;i+=4)
            _mm512_store_pd((double *)&padded[i],mw01_rader_a);
    } else if constexpr(std::is_same_v<float,Type>) {
        mw01_rader_a = _mm512_setzero_ps();
        for(i=0;i<newN;i+=8)
            _mm512_store_ps((float *)&padded[i],mw01_rader_a);
    }
#endif
    for(q=0;q<N-1;q++) padded[q] = datasub1[mapg[q]];
    fft_func<Type>(padded,result1,newN,1,1);
#if !defined(AVX) || AVX == 0
    memset(padded,0,newN*sizeof(complex<Type>));
#elif AVX512F == 0
    if constexpr(std::is_same_v<double,Type>) {
        mw01_rader_a = _mm256_setzero_pd();
        for(i=0;i<newN;i+=2)
            _mm256_store_pd((double *)&padded[i],mw01_rader_a);
    } else if constexpr(std::is_same_v<float,Type>) {
        mw01_rader_a = _mm256_setzero_ps();
        for(i=0;i<newN;i+=4)
            _mm256_store_ps((float *)&padded[i],mw01_rader_a);
    }
#else
    if constexpr(std::is_same_v<double,Type>) {
        mw01_rader_a = _mm512_setzero_pd();
        for(i=0;i<newN;i+=4)
            _mm512_store_pd((double *)&padded[i],mw01_rader_a);
    } else if constexpr(std::is_same_v<float,Type>) {
        mw01_rader_a = _mm512_setzero_ps();
        for(i=0;i<newN;i+=8)
            _mm512_store_ps((float *)&padded[i],mw01_rader_a);
    }
#endif
    for(int q=0;q<N-1;q++) padded[q] = datasub2[mapginv[q]];
#if !defined(AVX) || AVX == 0
    memcpy(&padded[newN-N+2],&padded[1],(N-2)*sizeof(complex<Type>));
    //for(q=N-2;q>0;q--) padded[q+newN-N+1] = padded[q];
#elif AVX512F == 0
    if constexpr(std::is_same_v<double,Type>) {
        if((N-3)%2 == 0)
            for(q=N-3;q>0;q-=2) _mm256_store_pd((double *)&padded[q+newN-N+1],_mm256_load_pd((double *)&padded[q]));
        else
            for(q=N-3;q>0;q-=2) _mm256_store_pd((double *)&padded[q+newN-N+1],_mm256_loadu_pd((double *)&padded[q]));
        q = q + 2 - 1;
    } else if constexpr(std::is_same_v<float,Type>) {
        if((N-5)%4 == 0)
            for(q=N-5;q>0;q-=4) _mm256_store_ps((float *)&padded[q+newN-N+1],_mm256_load_ps((float *)&padded[q]));
        else
            for(q=N-5;q>0;q-=4) _mm256_store_ps((float *)&padded[q+newN-N+1],_mm256_loadu_ps((float *)&padded[q]));
        q = q + 4 - 1;
    }
    for(;q>0;q--) padded[q+newN-N+1] = padded[q];
#else
    if constexpr(std::is_same_v<double,Type>) {
        if((N-5)%4 == 0)
            for(q=N-5;q>0;q-=4) _mm512_store_pd((double *)&padded[q+newN-N+1],_mm512_load_pd((double *)&padded[q]));
        else
            for(q=N-5;q>0;q-=4) _mm512_store_pd((double *)&padded[q+newN-N+1],_mm512_loadu_pd((double *)&padded[q]));
        q = q + 4 - 1;
    } else if constexpr(std::is_same_v<float,Type>) {
        if((N-9)%8 == 0)
            for(q=N-9;q>0;q-=8) _mm512_store_ps((float *)&padded[q+newN-N+1],_mm512_load_ps((float *)&padded[q]));
        else
            for(q=N-9;q>0;q-=8) _mm512_store_ps((float *)&padded[q+newN-N+1],_mm512_loadu_ps((float *)&padded[q]));
        q = q + 8 - 1;
    }
    for(;q>0;q--) padded[q+newN-N+1] = padded[q];
#endif


    
    fft_func<Type>(padded,result2,newN,1,1);
    newN2 = 1.*newN;
#if !defined(AVX) || AVX == 0 
    for(q=0;q<newN;q++) result1[q] *= result2[q]*newN2;
#elif AVX512F == 0
    if constexpr(std::is_same_v<double,Type>) {
        mw01_rader_a = _mm256_set1_pd(newN2);
        for(q=0;q<newN;q+=2)
            _mm256_store_pd((double *)&result1[q],_mm256_mul_pd(mw01_rader_a,complex_mul_256register((double *)&result1[q],(double *)&result2[q],0,0)));  
    } else if constexpr(std::is_same_v<float,Type>) {
        mw01_rader_a = _mm256_set1_ps(newN2); 
        for(q=0;q<newN;q+=4)
            _mm256_store_ps((float *)&result1[q],_mm256_mul_ps(mw01_rader_a,complex_mul_256register((float *)&result1[q],(float *)&result2[q],0,0)));  
    }     
#else
    if constexpr(std::is_same_v<double,Type>) {
        mw01_rader_a = _mm512_set1_pd(newN2);
        for(q=0;q<newN;q+=4)
            _mm512_store_pd((double *)&result1[q],_mm512_mul_pd(mw01_rader_a,complex_mul_512register((double *)&result1[q],(double *)&result2[q],0,0)));
    } else if constexpr(std::is_same_v<float,Type>) {
        mw01_rader_a = _mm512_set1_ps(newN2);
        for(q=0;q<newN;q+=8)
            _mm512_store_ps((float *)&result1[q],_mm512_mul_ps(mw01_rader_a,complex_mul_512register((float *)&result1[q],(float *)&result2[q],0,0)));
    }
#endif
    fftinv_func<Type>(result1,result2,newN); 
    for(p=0;p<N-1;p++) out[mapginv[p]] = result2[p];                                     // rearrange
    oldN = N;    
}


