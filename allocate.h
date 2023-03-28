template <typename T>
T **allocate2D(int i,int j) {
    T **pointer;
    int m;

    pointer = new T*[i];
    pointer[0] = new T[i*j];
    for(m=1;m<i;m++) pointer[m] = pointer[m-1] + j;

    return pointer;
}

template <typename T>
T ***allocate3D(int i,int j,int k) {
    T ***pointer;
    int m,n;
    int oldm,oldn;

    pointer = new T**[i];
    for(m=0;m<i;m++) pointer[m] = new T*[j];
    pointer[0][0] = new T[i*j*k];

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

template <typename T>
T ****allocate4D(int i,int j,int k,int l) {
    T ****pointer;
    int m,n,o;
    int oldm,oldn,oldo;

    pointer = new T***[i];
    for(m=0;m<i;m++) pointer[m] = new T**[j];
    for(m=0;m<i;m++) for(n=0;n<j;n++) pointer[m][n] = new T*[k];
    pointer[0][0][0] = new T[i*j*k*l];

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

template <typename T>
T *****allocate5D(int i,int j,int k,int l,int m) {
    T *****pointer;
    int n,o,p,q;
    int oldn,oldo,oldp,oldq;

    pointer = new T****[i];
    for(n=0;n<i;n++) pointer[n] = new T***[j];
    for(n=0;n<i;n++) for(o=0;o<j;o++) pointer[n][o] = new T**[k];
    for(n=0;n<i;n++) for(o=0;o<j;o++) for(p=0;p<k;p++) pointer[n][o][p] = new T*[l];
    pointer[0][0][0][0] = new T[i*j*k*l*m];

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
