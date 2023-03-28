#define maxN 65536
#define maxS 128

template <class Type>
void morlet(Type *data,complex<Type> **transform,int N,int S,Type param,Type dx,Type pi) {
    int static thread_local k,n,s;
    Type static thread_local a2[maxN][maxS],b2[maxN][maxS];
    Type static thread_local a,b;
    Type static thread_local freq[maxN];
    Type static thread_local scale[maxS];
    Type static thread_local wavefunc[maxS][maxN];
    complex<Type> static thread_local dft[maxN];
    complex<Type> static thread_local dft_product[maxS][maxN];
    complex<Type> static thread_local datacomplex[maxN];
//    #pragma omp threadprivate(k,n,s,a2,b2,a,b,freq,scale,wavefunc,dft,dft_product,datacomplex)

    for(s=0;s<S;s++) scale[s] = pow(2.,s/4.);                                        // set scales with dj=0.25

    for(k=0;k<N;k++)
        if(k <= N/2)
            freq[k] = 2.*pi*k/N/dx;
	else
            freq[k] = 2.*pi*(k-N)/N/dx; 

    for(s=0;s<S;s++)
    for(k=0;k<N;k++)
	if(freq[k] > 0.) {
            wavefunc[s][k] = exp(-0.5*pow(scale[s]*freq[k]-param,2.))/pow(pi,0.25);
	} else {
	    wavefunc[s][k] = exp(-0.5*pow(scale[s]*freq[k]-param,2.))/pow(pi,0.25);
	    //wavefunc[s][k] = 0.;
	}
    
    for(n=0;n<N;n++) datacomplex[n].setrealimga(data[n],0);
    
    dft_func(datacomplex,dft,N,1,pi,1);

    for(s=0;s<S;s++)
    for(k=0;k<N;k++) {
        dft_product[s][k] = dft[k]*wavefunc[s][k];
    }
    
    for(s=0;s<S;s++) {
	dftinv_func(dft_product[s],transform[s],N,pi);        // this also modifies dft_product[s]
    }

    for(s=0;s<S;s++) {
        a = sqrt(2.*pi*scale[s]/dx);
        for(n=0;n<N;n++) {
            transform[s][n] *= a;
        }
    }
}



