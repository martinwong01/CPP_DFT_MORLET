#define maxN 65536 
#define maxCooleyTukey 400                                       // if prime factor larger than this, use Rader algorithm

int smallfactor(int,int);

template <class Type>
void dft_func(complex<Type> *data,complex<Type> *out,int N,int Product,Type pi,int sign) {             // if called from main, set sign=1
    int static thread_local i,j,k,m,n,p,q,t;
    complex<Type> static thread_local datasub1[maxN];
    complex<Type> static thread_local datasub2[maxN];
    complex<Type> static thread_local datasub3[maxN];
    complex<Type> static thread_local c1[1];                                                                    // default constructor not called
    complex<Type> static thread_local c2[1];
    Type static thread_local a;
    int static thread_local PF,NoverPF;
    int static thread_local kleft,kright,mleft,nright,tail;
    complex<Type> static thread_local roots[maxN];
    int static thread_local Factor;
//    #pragma omp threadprivate(i,j,k,m,n,p,q,t,datasub1,datasub2,datasub3,c1,c2,a,PF,NoverPF,kleft,kright,mleft,nright,tail,roots,Factor) 

    roots[0].setrealimga(1.,0.);
    if(sign > 0)
        c1[0].setangle(-2.*pi/N);
    else
        c1[0].setangle(2.*pi/N);

    if(N%2 == 0) {
        for(i=1;i<N/2;i++) { roots[i] = roots[i-1]*c1[0]; roots[N-i].setrealimga(roots[i].realpart(),-roots[i].imgapart()); }
	roots[N/2].setrealimga(-1.,0.);
    } else {
        for(i=1;i<N/2+1;i++) { roots[i] = roots[i-1]*c1[0]; roots[N-i].setrealimga(roots[i].realpart(),-roots[i].imgapart()); }
    }


  while(1) {
    for(n=0;n<N;n++) out[n].setzero();
    Factor = smallfactor(N,Product);                                                  // find the smallest factor of N/Product

    PF = Product*Factor;
    NoverPF = N/PF;

    kleft = NoverPF;
    kright = N/Product; 
    mleft = N/Factor;                                                 // when m increases by 1, n increases by this
    nright = NoverPF;
    tail = NoverPF;

    if(Factor == 2) {
        for(k=0;k<Product;k++) {
            p = k*NoverPF; 
            for(t=0;t<tail;t++) {
                c2[0] = roots[p]*data[k*kright+nright+t];
                out[k*kleft+t] += data[k*kright+t] + c2[0];
		out[k*kleft+mleft+t] += data[k*kright+t] - c2[0]; 
            }
	}
    } else if(Factor <= maxCooleyTukey) {
	for(k=0;k<Product;k++) 
	    for(m=0;m<Factor;m++) {
	        j = k+m*Product;
		for(n=0;n<Factor;n++) {                                               // summation index
		    p = (n*j)%PF;
		    p = p*NoverPF;
                    for(t=0;t<tail;t++) {
                        out[k*kleft+m*mleft+t] += roots[p]*data[k*kright+n*nright+t];
		    }
		} 
            }
    } else {
        for(k=0;k<Product;k++) {
	    for(n=0;n<Factor;n++) {                                                       //    m=0, summation index
		p = n*k;
		p = p%PF;
		p = p*NoverPF;
                for(t=0;t<tail;t++) 
		    out[k*kleft+0*mleft+t] += roots[p]*data[k*kright+n*nright+t];
	    }

	    for(t=0;t<tail;t++) {                                                         //    m=1,.....
		for(q=1;q<Factor;q++) {
                    p = q*k;
                    p = p%PF;		
	            p = p*NoverPF;	    
                    datasub1[q] = roots[p]*data[k*kright+q*nright+t];
		    datasub2[q] = roots[((q*Product)%PF)*NoverPF];
		}
		Rader(datasub1,datasub2,datasub3,Factor,pi);
	        for(m=1;m<Factor;m++) 
		    out[k*kleft+m*mleft+t] = data[k*kright+t] + datasub3[m];	    
	    }
	}
    }

    if(PF != N) {
	memcpy(data,out,N*sizeof(complex<Type>));
    } else {
	a = 1./N;
        for(n=0;n<N;n++) out[n] *= a;
        break;	
    }

    Product = PF;
  }
}

template <class Type>
void dft_func2(complex<Type> *data,complex<Type> *out,int N,int Product,Type pi,int sign) {             // if called from main, set sign=1
    int static thread_local i,j,k,m,n,p,q,t;
    complex<Type> static thread_local datasub1[maxN];
    complex<Type> static thread_local datasub2[maxN];
    complex<Type> static thread_local datasub3[maxN];
    complex<Type> static thread_local c1[1];                                                                    // default constructor not called
    complex<Type> static thread_local c2[1];
    Type static thread_local a;
    int static thread_local PF,NoverPF;
    int static thread_local kleft,kright,mleft,nright,tail;
    complex<Type> static thread_local roots[maxN];
    int static thread_local Factor;
//    #pragma omp threadprivate(i,j,k,m,n,p,q,t,datasub1,datasub2,datasub3,c1,c2,a,PF,NoverPF,kleft,kright,mleft,nright,tail,roots,Factor) 

    roots[0].setrealimga(1.,0.);
    if(sign > 0)
        c1[0].setangle(-2.*pi/N);
    else
        c1[0].setangle(2.*pi/N);

    if(N%2 == 0) {
        for(i=1;i<N/2;i++) { roots[i] = roots[i-1]*c1[0]; roots[N-i].setrealimga(roots[i].realpart(),-roots[i].imgapart()); }
	roots[N/2].setrealimga(-1.,0.);
    } else {
        for(i=1;i<N/2+1;i++) { roots[i] = roots[i-1]*c1[0]; roots[N-i].setrealimga(roots[i].realpart(),-roots[i].imgapart()); }
    }

  while(1) {
    for(n=0;n<N;n++) out[n].setzero();
    Factor = smallfactor(N,Product);                                                  // find the smallest factor of N/Product

    PF = Product*Factor;
    NoverPF = N/PF;

    kleft = NoverPF;
    kright = N/Product; 
    mleft = N/Factor;                                                 // when m increases by 1, n increases by this
    nright = NoverPF;
    tail = NoverPF;

    if(Factor == 2) {
        for(k=0;k<Product;k++) {
            p = k*NoverPF; 
            for(t=0;t<tail;t++) {
                c2[0] = roots[p]*data[k*kright+nright+t];
                out[k*kleft+t] += data[k*kright+t] + c2[0];
	        out[k*kleft+mleft+t] += data[k*kright+t] - c2[0]; 
            }
	}
    } else {
        printf("ERROR Rader\n");
    }

    if(PF != N) {
	memcpy(data,out,N*sizeof(complex<Type>));
    } else {
	a = 1./N;
        for(n=0;n<N;n++) out[n] *= a;
        break;	
    }

    Product = PF;
  }
}

template <class Type>
void dftinv_func(complex<Type> *data,complex<Type> *out,int N,Type pi) {
    int i,j;
    complex<Type> NN(N,0);
    for(i=0;i<N;i++) data[i] *= NN;
    dft_func<Type>(data,out,N,1,pi,-1);
}

template <class Type>
void dftinv_func2(complex<Type> *data,complex<Type> *out,int N,Type pi) {
    int i,j;
    complex<Type> NN(N,0);
    for(i=0;i<N;i++) data[i] *= NN;
    dft_func2<Type>(data,out,N,1,pi,-1);
}

// find the smallest factor
int smallfactor(int n,int p) {
    n = n/p;
    if(n%2 == 0) return 2;
    for(int i=3;i*i<=n;i+=2)
        if(n%i == 0) return i;
    if(n > 2) return n;
    return -1;
}

int powmod(int a,int b,int p) {
    int res = 1;
    while (b)
        if (b&1)
            res = int(res*1ll*a%p), --b;       // 1ll   is long*long
        else
            a = int(a*1ll*a%p), b >>= 1;
    return res;
}

int generator(int p) {
    int fact[100];
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
    return x;
}

template <class Type>
void Rader(complex<Type> *datasub1,complex<Type> *datasub2,complex<Type> *out,int N,Type pi) {
    int static thread_local g,ginv;
    int static thread_local mapg[maxN];
    int static thread_local mapginv[maxN];
    int static thread_local newN;
    int static thread_local i;
    complex<Type> static thread_local padded1[maxN];
    complex<Type> static thread_local padded2[maxN];
    complex<Type> static thread_local result1[maxN];
    complex<Type> static thread_local result2[maxN];
//    #pragma omp threadprivate(g,ginv,mapg,mapginv,newN,padded1,padded2,result1,result2)

    g = generator(N);
    ginv = modulo_inverse(g,N);

    mapg[0] = 1;
    for(int q=1;q<N-1;q++) mapg[q] = mapg[q-1]*g%N;                                         // n = g^q    (mod N)  , q=[0,N-2]
    mapginv[0] = 1;
    for(int p=1;p<N-1;p++) mapginv[p] = mapg[N-p-1];                                        // k/m = ginv^p (mod N)  , p=[0,N-2]

    // padding to 2^
    newN = 2;
    for(i=0;i<1000000;i++) {
	if(newN >= 2*(N-1)-1) 
	    break;
        else 
            newN <<= 1;
    }

    for(int q=0;q<newN;q++) { padded1[q].setzero(); padded2[q].setzero(); }
    for(int q=0;q<N-1;q++) padded1[q] = datasub1[mapg[q]];
    for(int q=0;q<N-1;q++) padded2[q] = datasub2[mapginv[q]];
    for(int q=1;q<N-1;q++) padded2[newN-N+1+q] = padded2[q];
    dft_func2<Type>(padded1,result1,newN,1,pi,1);    
    dft_func2<Type>(padded2,result2,newN,1,pi,1);
    for(int q=0;q<newN;q++) result1[q] *= result2[q]*newN; 
    dftinv_func2<Type>(result1,result2,newN,pi); 
    for(int p=0;p<N-1;p++) out[mapginv[p]] = result2[p];                                     // rearrange    
}


