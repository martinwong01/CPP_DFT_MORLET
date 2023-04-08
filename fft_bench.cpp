#include <iostream>
#include <cmath>
#include <cstring>
#include "allocate.h"
#include "complex.h"
#include "dft.h"
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


  // test 3, large prime. read double data, use double to compute
  if(1 == 1) {
    i = 131072;
    data_d = new complex<double>[i]; 
    result1_d = new complex<double>[i]; 
    result2_d = new complex<double>[i]; 
    fp = fopen("data131072","r");
    for(int y=0;y<i;y++) {
        n = fscanf(fp,"%lf %lf",&a_d,&b_d);
	data_d[y].setrealimga(a_d,b_d);
    }
    fclose(fp);

    system("date");
    for(int y=0;y<1000;y++)
        dft_func<double>(data_d,result1_d,i,1,pi_d,1);
    system("date");


    result1_d[8191].print();


/*   
    printf("TRANSFORM\n");
    for(n=0;n<i;n++) {
        result1_d[n].print();
    }
*/

/*    
    dftinv_func<double>(result1_d,result2_d,i,pi_d);
    printf("INV TRANSFORM\n");
    for(n=0;n<i;n++) {
        result2_d[n].print();
    }
*/
    delete(data_d);
    delete(result1_d);
    delete(result2_d);
    
  }

}

