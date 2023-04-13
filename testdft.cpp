#include <iostream>
#include <cmath>
#include <cstring>
#include "allocate.h"
#include "complex.h"
#include "dft.h"
#define align 32
using namespace std;

int main(int argc, char *argv[]) {
    int i,j,k,n;
    FILE *fp;
   
    complex<double> *result1_d,*result2_d;
    double pi_d;
    double a_d,b_d;
    complex<double> *data_d;

    complex<float> *result1_f,*result2_f;
    float pi_f;
    float a_f,b_f;
    complex<float> *data_f;

    pi_d = atan(1.)*4.;
    pi_f = atan(1.)*4.;


  // test 1, double computations
  if(1 == 1) {
    data_d = allocate1D<complex<double>>(24,align);
    result1_d = allocate1D<complex<double>>(24,align);
    result2_d = allocate1D<complex<double>>(24,align);
    data_d[0].setrealimga(0.,1.);
    data_d[1].setrealimga(2.,-1.);
    data_d[2].setrealimga(1.,3.);
    data_d[3].setrealimga(-2.,-2.);
    data_d[4].setrealimga(-1.,2.5);
    data_d[5].setrealimga(-1.5,0.8);
    data_d[6].setrealimga(20.,10.);
    data_d[7].setrealimga(-7.,1.4);
    data_d[8].setrealimga(-0.02,3.);
    data_d[9].setrealimga(-4.3,-2.5);
    data_d[10].setrealimga(102.,0.7);
    data_d[11].setrealimga(-1.5555,0.35);
    data_d[12].setrealimga(0.33,1.33);
    data_d[13].setrealimga(2.66,-1.66);
    data_d[14].setrealimga(-1.27,5.09);
    data_d[15].setrealimga(-2.44,-2.88);
    data_d[16].setrealimga(1.369,2.28);
    data_d[17].setrealimga(1.5757,0.437);
    data_d[18].setrealimga(15.3,-4.333);
    data_d[19].setrealimga(-5.,1.1111);
    data_d[20].setrealimga(20.2,-4.5);
    data_d[21].setrealimga(-34.3,-15.5);
    data_d[22].setrealimga(21.777,0.39);
    data_d[23].setrealimga(201.55,1000.35);

    printf("test1\n");
    dft_func<double>(data_d,result1_d,24,1,pi_d,1);
    
    printf("TRANSFORM\n");  
    for(n=0;n<24;n++) {
        result1_d[n].print();
    }
    dftinv_func<double>(result1_d,result2_d,24,pi_d);

    printf("INV TRANSFORM\n");
    for(n=0;n<24;n++) {
        result2_d[n].print();
    }

    free(data_d);
    free(result1_d);
    free(result2_d);
  } 


  // test 2, float computations
  if(1 == 1) {
    data_f = allocate1D<complex<float>>(24,align);
    result1_f = allocate1D<complex<float>>(24,align);
    result2_f = allocate1D<complex<float>>(24,align);
    data_f[0].setrealimga(0.,1.);
    data_f[1].setrealimga(2.,-1.);
    data_f[2].setrealimga(1.,3.);
    data_f[3].setrealimga(-2.,-2.);
    data_f[4].setrealimga(-1.,2.5);
    data_f[5].setrealimga(-1.5,0.8);
    data_f[6].setrealimga(20.,10.);
    data_f[7].setrealimga(-7.,1.4);
    data_f[8].setrealimga(-0.02,3.);
    data_f[9].setrealimga(-4.3,-2.5);
    data_f[10].setrealimga(102.,0.7);
    data_f[11].setrealimga(-1.5555,0.35);
    data_f[12].setrealimga(0.33,1.33);
    data_f[13].setrealimga(2.66,-1.66);
    data_f[14].setrealimga(-1.27,5.09);
    data_f[15].setrealimga(-2.44,-2.88);
    data_f[16].setrealimga(1.369,2.28);
    data_f[17].setrealimga(1.5757,0.437);
    data_f[18].setrealimga(15.3,-4.333);
    data_f[19].setrealimga(-5.,1.1111);
    data_f[20].setrealimga(20.2,-4.5);
    data_f[21].setrealimga(-34.3,-15.5);
    data_f[22].setrealimga(21.777,0.39);
    data_f[23].setrealimga(201.55,1000.35);

    printf("test2\n");
    dft_func<float>(data_f,result1_f,24,1,pi_f,1);
    
    printf("TRANSFORM\n");  
    for(n=0;n<24;n++) {
        result1_f[n].print();
    }
    dftinv_func<float>(result1_f,result2_f,24,pi_f);
    printf("INV TRANSFORM\n");
    for(n=0;n<24;n++) {
        result2_f[n].print();
    }
    free(data_f);
    free(result1_f);
    free(result2_f);
  } 



  // test 3, large prime. read double data, use double to compute
  if(1 == 1) {
    i = 10007;
    data_d = allocate1D<complex<double>>(i,align); 
    result1_d = allocate1D<complex<double>>(i,align); 
    result2_d = allocate1D<complex<double>>(i,align); 
    fp = fopen("data10007","r");
    for(int y=0;y<i;y++) {
        n = fscanf(fp,"%lf %lf",&a_d,&b_d);
	data_d[y].setrealimga(a_d,b_d);
    }
    fclose(fp);

    printf("test3\n");
    dft_func<double>(data_d,result1_d,i,1,pi_d,1);
   
    printf("TRANSFORM\n");
    for(n=0;n<i;n++) {
        result1_d[n].print();
    }
    
    dftinv_func<double>(result1_d,result2_d,i,pi_d);
    printf("INV TRANSFORM\n");
    for(n=0;n<i;n++) {
        result2_d[n].print();
    }
    free(data_d);
    free(result1_d);
    free(result2_d);
    
  }

  // test 4, large prime. read double data, use float to compute
  if(1 == 1) {
    i = 10007;
    data_f = allocate1D<complex<float>>(i,align); 
    result1_f = allocate1D<complex<float>>(i,align); 
    result2_f = allocate1D<complex<float>>(i,align); 
    fp = fopen("data10007","r");
    for(int y=0;y<i;y++) {
        n = fscanf(fp,"%lf %lf",&a_d,&b_d);
	data_f[y].setrealimga((float)a_d,(float)b_d);
    }
    fclose(fp);

    printf("test4\n");
    dft_func<float>(data_f,result1_f,i,1,pi_f,1);
   
    printf("TRANSFORM\n");
    for(n=0;n<i;n++) {
        result1_f[n].print();
    }
    
    dftinv_func<float>(result1_f,result2_f,i,pi_f);
    printf("INV TRANSFORM\n");
    for(n=0;n<i;n++) {
        result2_f[n].print();
    }
    free(data_f);
    free(result1_f);
    free(result2_f);
    
  }


  // test 5, composite number, read double data, double compute  
  if(1 == 1) {
    i = 16071;
    data_d = allocate1D<complex<double>>(i,align); 
    result1_d = allocate1D<complex<double>>(i,align); 
    result2_d = allocate1D<complex<double>>(i,align); 
    fp = fopen("data16071","r");
    for(int y=0;y<i;y++) {
        n = fscanf(fp,"%lf %lf",&a_d,&b_d);
	data_d[y].setrealimga(a_d,b_d);
    }
    fclose(fp);

    printf("test5\n");
    dft_func<double>(data_d,result1_d,i,1,pi_d,1);
    
    printf("TRANSFORM\n");
    for(n=0;n<i;n++) {
        result1_d[n].print();
    }
    
    dftinv_func<double>(result1_d,result2_d,i,pi_d);
    printf("INV TRANSFORM\n");
    for(n=0;n<i;n++) {
        result2_d[n].print();
    }
    
    free(data_d);
    free(result1_d);
    free(result2_d);
  }


}

