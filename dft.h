#define maxN 65536 
#define maxCooleyTukey 100                                       // if prime factor larger than this, use Rader algorithm

int smallfactor(int,int);

template <class Type>
void Rader(complex<Type> *,complex<Type> *,complex<Type> *,int,Type);

template <class Type>
void dft_func(complex<Type> *data,complex<Type> *out,int N,int Product,Type pi,int sign) {             // if called from main, set sign=1
    int thread_local i,j,k,m,n,p,q,t;
    complex<Type> thread_local datasub1[maxN];
    complex<Type> thread_local datasub2[maxN];
    complex<Type> thread_local datasub3[maxN];
    complex<Type> thread_local c1[1];                                                                    // default constructor not called
    complex<Type> thread_local c2[1];
    complex<Type> thread_local c3[1];
    complex<Type> thread_local c4[1];
    Type thread_local a;
    int thread_local PF,NoverPF;
    int thread_local kleft,kright,mleft,nright,tail;
    complex<Type> thread_local roots[maxN];
    int thread_local Factor;
    complex<Type> thread_local datatemp[maxN];
    complex<Type> thread_local outtemp[maxN];
    int thread_local kkleft,kkright,mmleft,nnright;
    complex<Type> thread_local *dataptr,*outptr,*swapptr;

    memcpy(datatemp,data,N*sizeof(complex<Type>));
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

  outptr = outtemp;
  dataptr = datatemp;
  while(1) {
    for(n=0;n<N;n++) outptr[n].setzero();
    Factor = smallfactor(N,Product);                                                  // find the smallest factor of N/Product

    PF = Product*Factor;
    NoverPF = N/PF;

    kleft = NoverPF;
    kright = N/Product; 
    mleft = N/Factor;                                                 // when m increases by 1, n increases by this
    nright = NoverPF;
    tail = NoverPF;

    if(Factor == 2) {
        p = -NoverPF;
	kkleft = -kleft;
	kkright = -kright;
        for(k=0;k<Product;k++) {
            p += NoverPF;
	    kkleft += kleft;
	    kkright += kright;
            for(t=0;t<tail;t++) {
                c2[0] = roots[p]*dataptr[kkright+nright+t];
                outptr[kkleft+t] += dataptr[kkright+t] + c2[0];
		outptr[kkleft+mleft+t] += dataptr[kkright+t] - c2[0]; 
            }
	}
    } else if(Factor <= maxCooleyTukey) {
        kkleft = -kleft;
	kkright = -kright;
	for(k=0;k<Product;k++) {
	    kkleft += kleft;
	    kkright += kright;
	    nnright = -nright;
	    for(n=0;n<Factor;n++) {                                               // summation index
	        nnright += nright;
		i = n*k*NoverPF;
                for(t=0;t<tail;t++) {
		    c2[0] = roots[i]*dataptr[kkright+nnright+t];
		    outptr[kkleft+0*mleft+t] += c2[0];                            // m = 0
		    for(m=1;m<(Factor+1)/2;m++) {
		        j = n*m*mleft%N;
		        c3[0] = c2[0]*roots[j].realpart();
			c4[0] = c2[0].turnleft()*roots[j].imgapart();
			outptr[kkleft+m*mleft+t] += c3[0]+c4[0];
			outptr[kkleft+(Factor-m)*mleft+t] += c3[0]-c4[0];
		    }
		}
	    } 
	}
    } else {
        kkleft = -kleft;
	kkright = -kright;
        for(k=0;k<Product;k++) {
	    kkleft += kleft;
	    kkright += kright;
	    nnright = -nright;
	    for(n=0;n<Factor;n++) {                                                       //    m=0, summation index
		p = n*k*NoverPF;
		nnright += nright;
                for(t=0;t<tail;t++) 
		    outptr[kkleft+0*mleft+t] += roots[p]*dataptr[kkright+nnright+t];
	    }

	    for(t=0;t<tail;t++) {                                                         //    m=1,.....
		for(q=1;q<Factor;q++) {
	            p = q*k*NoverPF;	    
                    datasub1[q] = roots[p]*dataptr[kkright+q*nright+t];
		    datasub2[q] = roots[q*Product*NoverPF];
		}
		Rader<Type>(datasub1,datasub2,datasub3,Factor,pi);
	        for(m=1;m<Factor;m++) 
		    outptr[kkleft+m*mleft+t] = dataptr[kkright+t] + datasub3[m];	    
	    }
	}
    }

    if(PF != N) {
        std::swap(outptr,dataptr);
    } else {
        memcpy(out,outptr,N*sizeof(complex<Type>));
	if(sign > 0) {
            a = 1./N;
            for(n=0;n<N;n++) out[n] *= a;
	}
        break;	
    }

    Product = PF;
  }
}

template <class Type>
void dft_func2(complex<Type> *data,complex<Type> *out,int N,int Product,Type pi,int sign) {             // if called from main, set sign=1
    int thread_local i,j,k,m,n,p,q,t;
    complex<Type> thread_local datasub1[maxN];
    complex<Type> thread_local datasub2[maxN];
    complex<Type> thread_local datasub3[maxN];
    complex<Type> thread_local c1[1];                                                                    // default constructor not called
    complex<Type> thread_local c2[1];
    Type thread_local a;
    int thread_local PF,NoverPF;
    int thread_local kleft,kright,mleft,nright,tail;
    complex<Type> thread_local roots[maxN];
    int thread_local Factor;
    complex<Type> thread_local datatemp[maxN];
    complex<Type> thread_local outtemp[maxN];
    int thread_local kkleft,kkright,mmleft,nnright;
    complex<Type> thread_local *dataptr,*outptr,*swapptr;

    memcpy(datatemp,data,N*sizeof(complex<Type>));
    roots[0].setrealimga(1.,0.);
    if(sign > 0)
        c1[0].setangle(-2.*pi/N);
    else
        c1[0].setangle(2.*pi/N);
	
    k = N>>3;
    j = N>>2;
    for(i=1;i<=k;i++) roots[i] = roots[i-1]*c1[0];                   //     1/8 quadrant values
    if(sign > 0) {
        for(i=1;i<k;i++) roots[k-i] = roots[i].swap().reverse();   
    } else {
        for(i=1;i<k;i++) roots[k-i] = roots[i].swap();    
    }
    if(sign > 0) {
        for(i=0;i<j;i++) {
            roots[i+j] = roots[i].turnright();
	    roots[i+(j<<1)] = roots[i].reverse();
            roots[i+(j<<2)-j] = roots[i].turnleft();
        }
    } else {
        for(i=0;i<j;i++) {
            roots[i+j] = roots[i].turnleft();
	    roots[i+(j<<1)] = roots[i].reverse();
            roots[i+(j<<2)-j] = roots[i].turnright();
        }
    }

  outptr = outtemp;
  dataptr = datatemp;
  NoverPF = N;
  while(1) {
    for(n=0;n<N;n++) outptr[n].setzero();
    Factor = 2;                                                       // find the smallest factor of N/Product

    PF = Product<<1;
    NoverPF >>= 1;

    kleft = NoverPF;
    kright = NoverPF<<1; 
    mleft = N>>1;                                                 // when m increases by 1, n increases by this
    nright = NoverPF;
    tail = NoverPF;

    if(Factor == 2) {
        kkleft = -kleft;
	kkright = -kright;
	p = -NoverPF;
        for(k=0;k<Product;k++) {
	    kkleft += kleft;
	    kkright += kright;
            p += NoverPF; 
            for(t=0;t<tail;t++) {
                c2[0] = roots[p]*dataptr[kkright+nright+t];
                outptr[kkleft+t] += dataptr[kkright+t] + c2[0];
	        outptr[kkleft+mleft+t] += dataptr[kkright+t] - c2[0]; 
            }
	}
    } else {
        printf("ERROR Rader\n");
    }

    if(PF != N) {
        std::swap(outptr,dataptr);
    } else {
        memcpy(out,outptr,N*sizeof(complex<Type>));
	if(sign > 0) {
            a = 1./N;
            for(n=0;n<N;n++) out[n] *= a;
	}
        break;	
    }

    Product = PF;
  }
}

template <class Type>
void dftinv_func(complex<Type> *data,complex<Type> *out,int N,Type pi) {
    dft_func<Type>(data,out,N,1,pi,-1);
}

template <class Type>
void dftinv_func2(complex<Type> *data,complex<Type> *out,int N,Type pi) {
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
    int thread_local g,ginv;
    int thread_local mapg[maxN];
    int thread_local mapginv[maxN];
    int thread_local newN;
    int thread_local i;
    complex<Type> thread_local padded1[maxN];
    complex<Type> thread_local padded2[maxN];
    complex<Type> thread_local result1[maxN];
    complex<Type> thread_local result2[maxN];

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


