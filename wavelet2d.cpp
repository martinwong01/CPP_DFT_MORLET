#include <iostream>
#include <cmath>
#include "allocate.h"          
#include "complex.h"
#include "dft.h"
#include "morlet.h"
#include "omp.h"
using namespace std;


// usage:   main.exe inputfile outputfile Ny Nx Sy Sx dy dx
// input file dimension  [Ny][Nx]
// output file dimension [Sx][Nx][Sy][Ny][2]

int main(int argc, char *argv[]) {
    int Nx,Ny;                           
    int x,y,i,j,k;                         
    double dx,dy;

    int Sx,Sy;                          // number of scales (see morlet.h) 
    int s;                              // index for scale

    double param;                       // parameter for Morlet

    double pi;

    double **data;                      // data input.                                              Dim: [Ny][Nx]
    complex<double> ***xtransform;      // for each y, the transform in x direction.                Dim: [Ny][Sx][Nx] 

    complex<float> ***cwt1;             // output in float (save disk)                              Dim: [Nx][Sy][Ny] 
    complex<float> ***cwt2;
    double **yseries1;                  // dummy array.                                             Dim: [threads][Ny]
    double **yseries2;
    complex<double> ***transform_real;  // transform of Re{xtransform[:][Sx][Nx]} in y direction   Dim: [Nx][Sy][Ny]
    complex<double> ***transform_imga;  // transform of Im{xtransform[:][Sx][Nx]} in y direction   


    FILE *fp,*fp1,*fp2;
    char cmd[500];
    int recordno;

    int Nthreads,thread; 


    // set constants
    pi = atan(1.)*4.;
//  param = 6.;
    param = 16.*pi/9.-9./32./pi;

    printf("param = %lf\n",param); 

    // read parameters
    Ny = atoi(argv[3]);
    Nx = atoi(argv[4]);
    Sy = atoi(argv[5]);
    Sx = atoi(argv[6]);
    dy = atof(argv[7]);
    dx = atof(argv[8]);
    printf("%d %d %d %d %lf %lf\n",Ny,Nx,Sy,Sx,dy,dx);



    Nthreads = omp_get_max_threads();
    printf("THREADS = %d\n",Nthreads);

    xtransform = allocate3D<complex<double>>(Ny,Sx,Nx);
    yseries1 = allocate2D<double>(Nthreads,Ny);
    yseries2 = allocate2D<double>(Nthreads,Ny);
    transform_real = allocate3D<complex<double>>(Nx,Sy,Ny);
    transform_imga = allocate3D<complex<double>>(Nx,Sy,Ny);
    cwt1 = allocate3D<complex<float>>(Nx,Sy,Ny);
    cwt2 = allocate3D<complex<float>>(Nx,Sy,Ny);
    data = allocate2D<double>(Ny,Nx);



    // read file and prepare output file
    fp = fopen(argv[1],"rb");
    for(y=0;y<Ny;y++) i = fread(data[y],sizeof(double),Nx,fp);
    fclose(fp);



    sprintf(cmd,"%s1",argv[2]);
    fp1 = fopen(cmd,"rb+");
    sprintf(cmd,"%s2",argv[2]);
    fp2 = fopen(cmd,"rb+");


    #pragma omp parallel default(shared)
    {
        #pragma omp for
        for(y=0;y<Ny;y++) {                                                                    // xtransform for each y
	    morlet<double>(data[y],xtransform[y],Nx,Sx,param,dx,pi);
        }
    }


    for(s=0;s<Sx;s++) {                                                                     // each Sx 
	i = system("date");
	printf("s= %d\n",s);

        #pragma omp parallel default(shared) private(thread,y,i)
        {	 
        #pragma omp for   
	for(x=0;x<Nx;x++) {                                                                         // each Nx
            thread = omp_get_thread_num();

            for(y=0;y<Ny;y++) {
                yseries1[thread][y] = xtransform[y][s][x].realpart();
                yseries2[thread][y] = xtransform[y][s][x].imgapart();
	    }

            morlet<double>(yseries1[thread],transform_real[x],Ny,Sy,param,dy,pi);
            morlet<double>(yseries2[thread],transform_imga[x],Ny,Sy,param,dy,pi);

	    for(i=0;i<Sy;i++)
	    for(y=0;y<Ny;y++) {
                cwt1[x][i][y].setreal((float)(transform_real[x][i][y].realpart() - transform_imga[x][i][y].imgapart()));
                cwt1[x][i][y].setimga((float)(transform_real[x][i][y].imgapart() + transform_imga[x][i][y].realpart()));
            }
	    for(i=0;i<Sy;i++)
	    for(y=0;y<Ny;y++) {
                cwt2[x][i][y].setreal((float)(transform_real[x][i][y].realpart() + transform_imga[x][i][y].imgapart()));
                cwt2[x][i][y].setimga((float)(transform_imga[x][i][y].realpart() - transform_real[x][i][y].imgapart()));
            }
	}
        }
        fwrite(cwt1[0][0],Nx*Sy*Ny*sizeof(complex<float>),1,fp1);
        fwrite(cwt2[0][0],Nx*Sy*Ny*sizeof(complex<float>),1,fp2);
    }

    fclose(fp1);
    fclose(fp2);

}
