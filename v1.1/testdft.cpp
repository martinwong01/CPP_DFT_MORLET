#include "allocate.h"
#include "dft.h"

//   DFT usage:
//   dft_func(complex<Type> *data,complex<Type> *out,int N,int Product,int sign)
//   dftinv_func(complex<Type> *data,complex<Type> *out,int N)
//
//   Product always set to 1
//   sign is 1 for forward dft, -1 for inverse (dftinv_func simply call dft_func with sign -1)



//   here 6 tests are provided. Set the if statements to run them 




int main(int argc, char *argv[]) {
    int i,j,k,n;
    FILE *fp;
   
    complex<double> *result1_d,*result2_d;
    double a_d,b_d;
    complex<double> *data_d;

    complex<float> *result1_f,*result2_f;
    float a_f,b_f;
    complex<float> *data_f;



  // test 1, double computations
  if(1 == 1) {
    data_d = allocate1D<complex<double>>(24);
    result1_d = allocate1D<complex<double>>(24);
    result2_d = allocate1D<complex<double>>(24);
    data_d[0] = complex<double>(0.,1.);
    data_d[1] = complex<double>(2.,-1.);
    data_d[2] = complex<double>(1.,3.);
    data_d[3] = complex<double>(-2.,-2.);
    data_d[4] = complex<double>(-1.,2.5);
    data_d[5] = complex<double>(-1.5,0.8);
    data_d[6] = complex<double>(20.,10.);
    data_d[7] = complex<double>(-7.,1.4);
    data_d[8] = complex<double>(-0.02,3.);
    data_d[9] = complex<double>(-4.3,-2.5);
    data_d[10] = complex<double>(102.,0.7);
    data_d[11] = complex<double>(-1.5555,0.35);
    data_d[12] = complex<double>(0.33,1.33);
    data_d[13] = complex<double>(2.66,-1.66);
    data_d[14] = complex<double>(-1.27,5.09);
    data_d[15] = complex<double>(-2.44,-2.88);
    data_d[16] = complex<double>(1.369,2.28);
    data_d[17] = complex<double>(1.5757,0.437);
    data_d[18] = complex<double>(15.3,-4.333);
    data_d[19] = complex<double>(-5.,1.1111);
    data_d[20] = complex<double>(20.2,-4.5);
    data_d[21] = complex<double>(-34.3,-15.5);
    data_d[22] = complex<double>(21.777,0.39);
    data_d[23] = complex<double>(201.55,1000.35);

    printf("test1\n");
    dft_func<double>(data_d,result1_d,24,1,1);

    printf("TRANSFORM\n");  
    for(n=0;n<24;n++) {
        cout << result1_d[n] << endl;
    }

    dftinv_func<double>(result1_d,result2_d,24);    // no need to compute roots again

    printf("INV TRANSFORM\n");
    for(n=0;n<24;n++) {
        cout << result2_d[n] << endl;
    }

    free(data_d);
    free(result1_d);
    free(result2_d);
  } 



  // test 2, float computations
  if(1 == 0) {
    data_f = allocate1D<complex<float>>(24);
    result1_f = allocate1D<complex<float>>(24);
    result2_f = allocate1D<complex<float>>(24);
    data_f[0] = complex<float>(0.,1.);
    data_f[1] = complex<float>(2.,-1.);
    data_f[2] = complex<float>(1.,3.);
    data_f[3] = complex<float>(-2.,-2.);
    data_f[4] = complex<float>(-1.,2.5);
    data_f[5] = complex<float>(-1.5,0.8);
    data_f[6] = complex<float>(20.,10.);
    data_f[7] = complex<float>(-7.,1.4);
    data_f[8] = complex<float>(-0.02,3.);
    data_f[9] = complex<float>(-4.3,-2.5);
    data_f[10] = complex<float>(102.,0.7);
    data_f[11] = complex<float>(-1.5555,0.35);
    data_f[12] = complex<float>(0.33,1.33);
    data_f[13] = complex<float>(2.66,-1.66);
    data_f[14] = complex<float>(-1.27,5.09);
    data_f[15] = complex<float>(-2.44,-2.88);
    data_f[16] = complex<float>(1.369,2.28);
    data_f[17] = complex<float>(1.5757,0.437);
    data_f[18] = complex<float>(15.3,-4.333);
    data_f[19] = complex<float>(-5.,1.1111);
    data_f[20] = complex<float>(20.2,-4.5);
    data_f[21] = complex<float>(-34.3,-15.5);
    data_f[22] = complex<float>(21.777,0.39);
    data_f[23] = complex<float>(201.55,1000.35);

    printf("test2\n");
    dft_func<float>(data_f,result1_f,24,1,1);
    
    printf("TRANSFORM\n");  
    for(n=0;n<24;n++) {
        cout << result1_f[n] << endl;
    }
    dftinv_func<float>(result1_f,result2_f,24);
    printf("INV TRANSFORM\n");
    for(n=0;n<24;n++) {
        cout << result2_f[n] << endl;
    }
    free(data_f);
    free(result1_f);
    free(result2_f);
  } 



  // test 3, large prime. read double data, use double to compute
  if(1 == 0) {
    i = 10007;
    data_d = allocate1D<complex<double>>(i); 
    result1_d = allocate1D<complex<double>>(i); 
    result2_d = allocate1D<complex<double>>(i); 
    fp = fopen("data10007","r");
    for(int y=0;y<i;y++) {
        n = fscanf(fp,"%lf %lf",&a_d,&b_d);
	data_d[y] = complex<double>(a_d,b_d);
    }
    fclose(fp);

    printf("test3\n");
    dft_func<double>(data_d,result1_d,i,1,1);
   
    printf("TRANSFORM\n");
    for(n=0;n<i;n++) {
        cout << result1_d[n] << endl;
    }
    
    dftinv_func<double>(result1_d,result2_d,i);
    printf("INV TRANSFORM\n");
    for(n=0;n<i;n++) {
        cout << result2_d[n] << endl;
    }
    free(data_d);
    free(result1_d);
    free(result2_d);
    
  }

  // test 4, large prime. read double data, use float to compute
  if(1 == 0) {
    i = 10007;
    data_f = allocate1D<complex<float>>(i); 
    result1_f = allocate1D<complex<float>>(i); 
    result2_f = allocate1D<complex<float>>(i); 
    fp = fopen("data10007","r");
    for(int y=0;y<i;y++) {
        n = fscanf(fp,"%lf %lf",&a_d,&b_d);
	data_f[y] = complex<float>((float)a_d,(float)b_d);
    }
    fclose(fp);

    printf("test4\n");
    dft_func<float>(data_f,result1_f,i,1,1);
   
    printf("TRANSFORM\n");
    for(n=0;n<i;n++) {
        cout << result1_f[n] << endl;
    }
    
    dftinv_func<float>(result1_f,result2_f,i);
    printf("INV TRANSFORM\n");
    for(n=0;n<i;n++) {
        cout << result2_f[n] << endl;
    }
    free(data_f);
    free(result1_f);
    free(result2_f);
    
  }


  // test 5, composite number, read double data, double compute  
  if(1 == 0) {
    i = 16071;
    data_d = allocate1D<complex<double>>(i); 
    result1_d = allocate1D<complex<double>>(i); 
    result2_d = allocate1D<complex<double>>(i); 
    fp = fopen("data16071","r");
    for(int y=0;y<i;y++) {
        n = fscanf(fp,"%lf %lf",&a_d,&b_d);
	data_d[y] = complex<double>(a_d,b_d);
    }
    fclose(fp);

    printf("test5\n");
    dft_func<double>(data_d,result1_d,i,1,1);
    
    printf("TRANSFORM\n");
    for(n=0;n<i;n++) {
        cout << result1_d[n] << endl;
    }
    
    dftinv_func<double>(result1_d,result2_d,i);
    printf("INV TRANSFORM\n");
    for(n=0;n<i;n++) {
        cout << result2_d[n] << endl;
    }
    
    free(data_d);
    free(result1_d);
    free(result2_d);
  }



  // test 6, composite number, read double data, float compute  
  if(1 == 0) {
    i = 16071;
    data_f = allocate1D<complex<float>>(i); 
    result1_f = allocate1D<complex<float>>(i); 
    result2_f = allocate1D<complex<float>>(i); 
    fp = fopen("data16071","r");
    for(int y=0;y<i;y++) {
        n = fscanf(fp,"%lf %lf",&a_d,&b_d);
	data_f[y] = complex<float>((float)a_d,(float)b_d);
    }
    fclose(fp);

    printf("test6\n");
    dft_func<float>(data_f,result1_f,i,1,1);
    
    printf("TRANSFORM\n");
    for(n=0;n<i;n++) {
        cout << result1_f[n] << endl;
    }
    
    dftinv_func<float>(result1_f,result2_f,i);
    printf("INV TRANSFORM\n");
    for(n=0;n<i;n++) {
        cout << result2_f[n] << endl;
    }
    
    free(data_f);
    free(result1_f);
    free(result2_f);
  }


}

