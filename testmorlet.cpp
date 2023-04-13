#include <iostream>
#include <fstream>
#include <cmath>
#include "allocate.h"
#include "complex.h"
#include "dft.h"
#include "morlet.h"
#define align 32
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

        data = allocate1D<double>(Nx,align);
        data2 = allocate1D<double>(Nx,align);
        xtransform = allocate2D<complex<double>>(Sx,Nx,align);
    
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
        free(data);
        free(data2);
        free(xtransform[0]);        
    }


    if(1 == 1) {
        pi = atan(1.)*4.;
        param = 5.495528950892663;
        Nx = 144;
        Sx = 33;
        dx = 1.;

        data = allocate1D<double>(Nx,align);
        xtransform = allocate2D<complex<double>>(Sx,Nx,align);
   
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

        free(data);
        free(xtransform[0]);        
    }
}
