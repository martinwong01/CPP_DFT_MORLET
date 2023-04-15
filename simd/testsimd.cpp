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
    size_t num_simd_elements = 8;

    for(j=0;j<10000;j++) {
        for(;n-i>=num_simd_elements;i+=num_simd_elements) {
            __m256 x_vals = _mm256_load_ps(&a[i]);                      // simd intrinsics
            __m256 y_vals = _mm256_load_ps(&b[i]);
            __m256 z_vals = _mm256_add_ps(x_vals,y_vals);
            _mm256_store_ps(&c[i],z_vals);
        }
    }

    printf("%f\n",c[0]);

}
