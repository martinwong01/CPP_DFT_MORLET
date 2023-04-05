#define maxN 65536
#define maxS 128

template <class Type>
void morlet(Type *data,complex<Type> **transform,int N,int S,Type param,Type dx,Type pi) {
    int thread_local k,n,s;
    Type thread_local a2[maxN][maxS],b2[maxN][maxS];
    Type thread_local a,b;
    Type thread_local freq[maxN];
    Type thread_local scale[maxS];
    Type thread_local wavefunc[maxS][maxN];
    complex<Type> thread_local dft[maxN];
    complex<Type> thread_local dft_product[maxS][maxN];
    complex<Type> thread_local datacomplex[maxN];

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
    
    dft_func<Type>(datacomplex,dft,N,1,pi,1);

    for(s=0;s<S;s++)
    for(k=0;k<N;k++) {
        dft_product[s][k] = dft[k]*wavefunc[s][k];
    }
    
    for(s=0;s<S;s++) {
	dftinv_func<Type>(dft_product[s],transform[s],N,pi);        // this also modifies dft_product[s]
    }

    for(s=0;s<S;s++) {
        a = sqrt(2.*pi*scale[s]/dx);
        for(n=0;n<N;n++) {
            transform[s][n] *= a;
        }
    }
}

template <class Type>
Type integral_reconstruction(Type param) {
    Type thread_local sum;
    Type thread_local localbnd = 0.0001;
    Type thread_local upperbnd = 10000.;
    Type thread_local dx = localbnd;
    unsigned int thread_local intervals;
    unsigned int i;
    Type thread_local x;
    Type thread_local pi = atan(1.)*4.;

    intervals = (unsigned int)(upperbnd/localbnd);
    sum = 0.;

    for(i=1;i<=intervals;i++) {
        x = i*localbnd;
        sum += exp(-0.5*(x-param)*(x-param))/x;
    }
    sum *= pow(pi,-0.25)*dx;
    return sum;
}

template <class Type>
Type integral_covariance(Type param) {
    Type thread_local sum;
    Type thread_local localbnd = 0.0001;
    Type thread_local upperbnd = 10000.;
    Type thread_local dx = localbnd;
    unsigned int thread_local intervals;
    unsigned int i;
    Type thread_local x;
    Type thread_local pi = atan(1.)*4.;

    intervals = (unsigned int)(upperbnd/localbnd);
    sum = 0.;

    for(i=1;i<=intervals;i++) {
        x = i*localbnd;
        sum += exp(-(x-param)*(x-param))/x;
    }
    sum *= pow(pi,-0.5)*dx;
    return sum;
}

