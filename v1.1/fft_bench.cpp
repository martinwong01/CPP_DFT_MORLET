#include "allocate.h"
#include "dft.h"



int main(int argc, char *argv[]) {
    int i,j,k,n;
    FILE *fp;
   
    complex<double> *result1_d,*result2_d;
    double a_d,b_d;
    complex<double> *data_d;

    complex<double> data_d_stack[131072];
    complex<double> result1_d_stack[131072]; 
    complex<double> result2_d_stack[131072]; 

    complex<float> *result1_f,*result2_f;
    float a_f,b_f;
    complex<float> *data_f;



    if(1 == 1) {
        i = 131072;
        data_d = allocate1D<complex<double>>(i); 
        result1_d = allocate1D<complex<double>>(i); 
        result2_d = allocate1D<complex<double>>(i); 
        fp = fopen("data131072","r");
        for(int y=0;y<i;y++) {
            n = fscanf(fp,"%lf %lf",&a_d,&b_d);
	    data_d[y].setrealimag(a_d,b_d);
        }
        fclose(fp);


        system("date");
        #pragma omp parallel for default(shared)
        for(int y=0;y<5000;y++)
            fft_func<double>(data_d,result1_d,i,1,1);     //   FFT
        system("date");
        result1_d[8191].print();

        free(data_d);
        free(result1_d);
        free(result2_d);
  }


}

