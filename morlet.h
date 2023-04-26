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

    if(init == 1) {
        a = 2.*pi/N/dx;
        b = 1./pow(pi,0.25);
        for(s=0;s<S;s++) scale[s] = pow(2.,s/4.);                                        // set scales with dj=0.25
        for(k=0;k<N;k++)
            if(k <= N/2)
                freq[k] = a*k;
	    else
                freq[k] = a*(k-N); 
        for(s=0;s<S;s++)
        for(k=0;k<N;k++)
            if(freq[k] > 0.)
                wavefunc[s][k] = exp(-0.5*pow(scale[s]*freq[k]-param,2.))*b;
	    else
	        wavefunc[s][k] = exp(-0.5*pow(scale[s]*freq[k]-param,2.))*b;
	        //wavefunc[s][k] = 0.;
    }
    
    for(n=0;n<N;n++) datacomplex[n].setrealimga(data[n],0);
    
    dft_func<Type>(datacomplex,dft,N,1,pi,1,init);

    for(s=0;s<S;s++)
        for(k=0;k<N;k++) dft_product[s][k] = dft[k]*wavefunc[s][k];
    
    for(s=0;s<S;s++) dftinv_func<Type>(dft_product[s],transform[s],N,pi,0);

    for(s=0;s<S;s++) {
        a = sqrt(2.*pi*scale[s]/dx);
        for(n=0;n<N;n++) transform[s][n] *= a;
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

