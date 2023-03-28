#include <iostream>
#include <fstream>
#include <cmath>
#include "complex.h"
#include "dft.h"
#include "morlet.h"
using namespace std;


int main(int argc, char *argv[]) {
    int Nx;                           
    int x,i;                               
    double dx;
    int Sx;                                //   number of scales
    int s;                                 
    int k;                                 //   index for frequency
    double param;                          //   parameter for Morlet
    double pi;
    double *data;                          //   Dim: [Nx]
    double *data2;                         //   complex part. not used 
    complex<double> **xtransform;          //   Dim: [Sx][Nx]     
    ifstream file;

    if(1 == 1) {
        pi = atan(1.)*4.;
        param = 6.;
        Nx = 10007;
        Sx = 29;
        dx = 1.;

        data = new double[Nx];
        data2 = new double[Nx];
        xtransform = new complex<double>*[Sx];
        for(s=0;s<Sx;s++) xtransform[s] = new complex<double>[Nx];
    
        // read data
        file.open("data10007");
        for(x=0;x<Nx;x++) file >> data[x] >> data2[x];
        file.close();
   
        morlet<double>(data,xtransform,Nx,Sx,param,dx,pi);

        printf("test1 data\n");
	
	for(s=0;s<Sx;s++)
	for(x=7000;x<=7000;x++) {
            printf("S X = %d %d\n",s,x);
            xtransform[s][x].print();
	}
        delete(data);
        delete(data2);
        for(s=0;s<Sx;s++) delete(xtransform[s]);
        delete(xtransform);        
    }


    if(1 == 1) {
        pi = atan(1.)*4.;
        param = 5.495528950892663;
        Nx = 144;
        Sx = 33;
        dx = 1.;

        data = new double[Nx];
        xtransform = new complex<double>*[Sx];
        for(s=0;s<Sx;s++) xtransform[s] = new complex<double>[Nx];
   
        // read data
        file.open("data144");
        for(x=0;x<Nx;x++) file >> data[x];
        file.close();

        morlet<double>(data,xtransform,Nx,Sx,param,dx,pi);

        printf("test2 data\n");
	for(int s=0;s<Sx;s++)
	for(x=0;x<=0;x++) {
	    printf("%d %d\n",s,x);	
            xtransform[s][x].print();
	}

        delete(data);
        for(s=0;s<Sx;s++) delete(xtransform[s]);
        delete(xtransform);        
    }
}
