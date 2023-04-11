#include <iostream>
#include <immintrin.h>
#include <cstdio>


int main() {
    size_t n = 1000000000;
    float a[n];
    float b[n];
    float c[n];

    size_t i = 0;

    for(i=0;i<n;i++) {
        c[i] = a[i] + b[i];
    }

    printf("%f\n",c[0]);


/*
    n = ((uintptr_t)(&a[0]))%64;
    printf("%d\n",n);
    n = ((uintptr_t)(&b[0]))%64;
    printf("%d\n",n);
    n = ((uintptr_t)(&c[0]))%64;
    printf("%d\n",n);
*/
}
