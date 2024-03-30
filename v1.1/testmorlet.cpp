#include <fstream>
#include "allocate.h"
#include "morlet.h"

// Morlet Usage:
// morlet(Type *data,complex<Type> **transform,int N,int S,Type param,Type dx)
//
// data:  real data array 


int main(int argc, char *argv[]) {
    int Nx;                           
    int x,i;                               
    double dx;
    int Sx;                                //   number of scales
    int s;                                 
    int k;                                 //   index for frequency
    double param;                          //   parameter for Morlet
    double *data;                          //   Dim: [Nx]
    double *data2;                         //   complex part. not used 
    complex<double> **xtransform;          //   Dim: [Sx][Nx]     
    ifstream file;

    if(1 == 1) {
        param = 6.;
        Nx = 10007;
        Sx = 29;
        dx = 1.;

        data = allocate1D<double>(Nx);
        data2 = allocate1D<double>(Nx);
        xtransform = allocate2D<complex<double>>(Sx,Nx);
    
        // read data
        file.open("data10007");
        for(x=0;x<Nx;x++) file >> data[x] >> data2[x];
        file.close();
   
        morlet<double>(data,xtransform,Nx,Sx,param,dx);

        printf("test1 data\n");
	
	for(s=0;s<Sx;s++)
	for(x=7000;x<=7000;x++) {
            printf("S X = %d %d\n",s,x);
            cout << xtransform[s][x] << endl;
	}
        free(data);
        free(data2);
        free(xtransform[0]);        
    }


    if(1 == 1) {
        param = 5.495528950892663;
        Nx = 144;
        Sx = 33;
        dx = 1.;

        data = allocate1D<double>(Nx);
        xtransform = allocate2D<complex<double>>(Sx,Nx);
   
        // read data
        file.open("data144");
        for(x=0;x<Nx;x++) file >> data[x];
        file.close();

        morlet<double>(data,xtransform,Nx,Sx,param,dx);

        printf("test2 data\n");
	for(int s=0;s<Sx;s++)
	for(x=0;x<=0;x++) {
	    printf("%d %d\n",s,x);	
            cout << xtransform[s][x] << endl;
	}

        free(data);
        free(xtransform[0]);        
    }
}
