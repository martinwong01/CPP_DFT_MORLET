#include <cstring>

template <typename Type>
Type *allocate1D(int i,int align) {
    Type *pointer;

    pointer = (Type *)aligned_alloc(align,i*sizeof(Type));
    return pointer;
}

template <typename Type>
Type **allocate2D(int i,int j,int align) {
    Type **pointer;
    int m;

    pointer = (Type **)aligned_alloc(align,i*sizeof(Type*));
    pointer[0] = (Type *)aligned_alloc(align,i*j*sizeof(Type));
    memset(pointer[0],0,i*j*sizeof(Type));
    for(m=1;m<i;m++) pointer[m] = pointer[m-1] + j;
    return pointer;
}

template <typename Type>
Type ***allocate3D(int i,int j,int k,int align) {
    Type ***pointer;
    int m,n;
    int oldm,oldn;

    pointer = (Type ***)aligned_alloc(align,i*sizeof(Type**));
    for(m=0;m<i;m++) pointer[m] = (Type **)aligned_alloc(align,j*sizeof(Type*));
    pointer[0][0] = (Type *)aligned_alloc(align,i*j*k*sizeof(Type));
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
Type ****allocate4D(int i,int j,int k,int l,int align) {
    Type ****pointer;
    int m,n,o;
    int oldm,oldn,oldo;

    pointer = (Type ****)aligned_alloc(align,i*sizeof(Type***));
    for(m=0;m<i;m++) pointer[m] = (Type ***)aligned_alloc(align,j*sizeof(Type**));
    for(m=0;m<i;m++) for(n=0;n<j;n++) pointer[m][n] = (Type **)aligned_alloc(align,k*sizeof(Type*));
    pointer[0][0][0] = (Type *)aligned_alloc(align,i*j*k*l*sizeof(Type));  
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
Type *****allocate5D(int i,int j,int k,int l,int m,int align) {
    Type *****pointer;
    int n,o,p,q;
    int oldn,oldo,oldp,oldq;

    pointer = (Type *****)aligned_alloc(align,i*sizeof(Type****)); 
    for(n=0;n<i;n++) pointer[n] = (Type ****)aligned_alloc(align,j*sizeof(Type***)); 
    for(n=0;n<i;n++) for(o=0;o<j;o++) pointer[n][o] = (Type ***)aligned_alloc(align,k*sizeof(Type**));
    for(n=0;n<i;n++) for(o=0;o<j;o++) for(p=0;p<k;p++) pointer[n][o][p] = (Type **)aligned_alloc(align,l*sizeof(Type*));
    pointer[0][0][0][0] = (Type *)aligned_alloc(align,i*j*k*l*m*sizeof(Type));
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
