#include <immintrin.h>
#include <cstdio>


int main() {
    size_t n = 1000000000;
    float a[n];
    float b[n];
    float c[n];

    size_t i = 0;
    size_t num_simd_elements = 8;

    for(;n-i>=num_simd_elements;i+=num_simd_elements) {
        __m256 x_vals = _mm256_loadu_ps(&a[i]);                      // simd intrinsics
        __m256 y_vals = _mm256_loadu_ps(&b[i]);
        __m256 z_vals = _mm256_add_ps(x_vals,y_vals);
        _mm256_storeu_ps(&c[i],z_vals);
    }

    printf("%f\n",c[0]);

}
