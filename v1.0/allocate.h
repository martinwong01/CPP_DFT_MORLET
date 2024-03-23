// these subroutines return pointers to multidimensional C-style arrays. Last dimension
// is modified so that the vectors begin and end at the specified alignment (ALIGN)
// see wavelet2d.cpp for usage including file read/write 
#include <cstdlib>
#ifndef ALIGN
    #define ALIGN 64
#endif
#include <cstring>


template <typename Type>
int aligned_dim(int i) {
    int j;
    j = i+(ALIGN-i*sizeof(Type)%ALIGN)/sizeof(Type);    
    return j;
}

template <typename Type>
Type *allocate1D(int i) {
    Type *pointer;
 
    i = aligned_dim<Type>(i);
    pointer = (Type *)aligned_alloc(ALIGN,i*sizeof(Type));
    memset(pointer,0,i*sizeof(Type));
    return pointer;
}

template <typename Type>
Type **allocate2D(int i,int j) {
    Type **pointer;
    int m;
    
    j = aligned_dim<Type>(j);
    pointer = (Type **)aligned_alloc(ALIGN,i*sizeof(Type*));
    pointer[0] = (Type *)aligned_alloc(ALIGN,i*j*sizeof(Type));
    memset(pointer[0],0,i*j*sizeof(Type));
    for(m=1;m<i;m++) pointer[m] = pointer[m-1] + j;
    return pointer;
}

template <typename Type>
Type ***allocate3D(int i,int j,int k) {
    Type ***pointer;
    int m,n;
    int oldm,oldn;

    k = aligned_dim<Type>(k);
    pointer = (Type ***)aligned_alloc(ALIGN,i*sizeof(Type**));
    for(m=0;m<i;m++) pointer[m] = (Type **)aligned_alloc(ALIGN,j*sizeof(Type*));
    pointer[0][0] = (Type *)aligned_alloc(ALIGN,i*j*k*sizeof(Type));
    memset(pointer[0][0],0,i*j*k*sizeof(Type));

    oldm = 0;
    oldn = 0;
    for(m=0;m<i;m++)
    for(n=0;n<j;n++) {
        if(m == 0 && n == 0) continue;
        pointer[m][n] = pointer[oldm][oldn] + k;
        oldm = m;
        oldn = n;
    }
    return pointer;
}

template <typename Type>
Type ****allocate4D(int i,int j,int k,int l) {
    Type ****pointer;
    int m,n,o;
    int oldm,oldn,oldo;

    l = aligned_dim<Type>(l);
    pointer = (Type ****)aligned_alloc(ALIGN,i*sizeof(Type***));
    for(m=0;m<i;m++) pointer[m] = (Type ***)aligned_alloc(ALIGN,j*sizeof(Type**));
    for(m=0;m<i;m++) for(n=0;n<j;n++) pointer[m][n] = (Type **)aligned_alloc(ALIGN,k*sizeof(Type*));
    pointer[0][0][0] = (Type *)aligned_alloc(ALIGN,i*j*k*l*sizeof(Type));  
    memset(pointer[0][0][0],0,i*j*k*l*sizeof(Type));

    oldm = 0;
    oldn = 0;
    oldo = 0;
    for(m=0;m<i;m++)
    for(n=0;n<j;n++)
    for(o=0;o<k;o++) {
        if(m == 0 && n == 0 && o == 0) continue;
	pointer[m][n][o] = pointer[oldm][oldn][oldo] + l;
	oldm = m;
	oldn = n;
	oldo = o;
    }
    return pointer;
}

template <typename Type>
Type *****allocate5D(int i,int j,int k,int l,int m) {
    Type *****pointer;
    int n,o,p,q;
    int oldn,oldo,oldp,oldq;

    m = aligned_dim<Type>(m);
    pointer = (Type *****)aligned_alloc(ALIGN,i*sizeof(Type****)); 
    for(n=0;n<i;n++) pointer[n] = (Type ****)aligned_alloc(ALIGN,j*sizeof(Type***)); 
    for(n=0;n<i;n++) for(o=0;o<j;o++) pointer[n][o] = (Type ***)aligned_alloc(ALIGN,k*sizeof(Type**));
    for(n=0;n<i;n++) for(o=0;o<j;o++) for(p=0;p<k;p++) pointer[n][o][p] = (Type **)aligned_alloc(ALIGN,l*sizeof(Type*));
    pointer[0][0][0][0] = (Type *)aligned_alloc(ALIGN,i*j*k*l*m*sizeof(Type));
    memset(pointer[0][0][0][0],0,i*j*k*l*m*sizeof(Type));

    oldn = 0;
    oldo = 0;
    oldp = 0;
    oldq = 0;
    for(n=0;n<i;n++)
    for(o=0;o<j;o++)
    for(p=0;p<k;p++)
    for(q=0;q<l;q++) {
        if(n == 0 && o == 0 && p == 0 && q == 0) continue;
	pointer[n][o][p][q] = pointer[oldn][oldo][oldp][oldq] + m;
	oldn = n;
	oldo = o;
	oldp = p;
	oldq = q;
    }
    return pointer;
}
