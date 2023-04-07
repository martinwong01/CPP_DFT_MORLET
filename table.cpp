#include <iostream>
#include <complex>
#include <cmath>
#include <iterator>

using namespace std;

unsigned int bitReverse(unsigned int x, int log2n) {
  int n = 0;
  int mask = 0x1;
  for (int i=0; i < log2n; i++) {
    n <<= 1;
    n |= (x & 1);
    x >>= 1;
  }
  return n;
}
  
int main() {
  int maxn = 131072;
  int m,n,i;

  m = log2(maxn);

  for(int n=0;n<=m;n++) {
    i = 1;
    for(int j=1;j<=n;j++) i<<=1;
    printf("int table%d[] = {",n);
    for(int j=0;j<i-1;j++) printf("%d,",bitReverse(j,n));
    printf("%d};\n",bitReverse(i-1,n));
  }

  printf("\n");
  printf("int *table[] = {");

  for(int n=0;n<=m-1;n++) {
    printf("table%d,",n);
  }
  printf("table%d};\n",m);
}
