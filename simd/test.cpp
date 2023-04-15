#include <iostream>
#include <immintrin.h>
#include <cstdio>
#define align 32


int main() {
    size_t n = 100000;
    alignas(align) float a[n];
    alignas(align) float b[n];
    alignas(align) float c[n];

    size_t i = 0;
    size_t j = 0;

    for(j=0;j<10000;j++) {
        for(i=0;i<n;i++) {
            c[i] = a[i] + b[i];
        }
    }

    printf("%f\n",c[0]);
}
