#include <iostream>
#include <cmath>
#include "allocate.h"
#include "complex.h"
#include "dft.h"
#include "morlet.h"
#include "omp.h"
using namespace std;


int main(int argc, char *argv[]) {
    int Nx;                           
    int Ny;
    int x,y,i,j,k;                         // index for position
    double dx;
    double dy;

    int Sx;                                // number of scales in x  
    int Sy;
    int s;                                 // index for scale

    int year,month;
    int startyear = 1979;
    int endyear = 2022;
    int tiles = 12;
    int startsec,sec1,sec2;
    int startindex[endyear-startyear+1][12];
    int endindex[endyear-startyear+1][12];
    int counter;

    double a;
    double b;
    double pi;
    double param;

    // number of threads
    complex<float> ***cwt1;                     // input                                      Dim: [Nx][Sy][Ny]
    complex<float> ***cwt2;
    double ***cov;
    double *****data;
    double **data_global;

    FILE *fp,*fp1,*fp2;
    char cmd[500];
    int recordno;
    int Nthreads,thread; 
    char buffer[256];
    std::string result = "";


    // set constants
    pi = atan(1.)*4.;
    param = 16.*pi/9.-9./32./pi;

    Nthreads = omp_get_max_threads();
    printf("THREADS = %d\n",Nthreads);


    Ny = atoi(argv[3]);
    Nx = atoi(argv[4]);
    Sy = atoi(argv[5]);
    Sx = atoi(argv[6]);
    dy = atof(argv[7]);
    dx = atof(argv[8]);


    cwt1 = allocate3D<complex<float>>(Nx,Sy,Ny);
    cwt2 = allocate3D<complex<float>>(Nx,Sy,Ny);
    cov = allocate3D<double>(Nx,Sy,Ny);
    data = allocate5D<double>(endyear-startyear+1,12,tiles,Sx,Sy);
    data_global = allocate2D<double>(Sx,Sy);


    fp = popen("date --date='19790101' +%s","r");
    while(fgets(buffer,sizeof buffer,fp) != NULL) result += buffer;
    startsec = stoi(result);

    for(year=startyear;year<=endyear;year++) 
    for(month=1;month<=12;month++) {
        sprintf(cmd,"date --date='%d%02d01' +%s",year,month,"%s");
	result = "";
	fp = popen(cmd,"r");
	while(fgets(buffer,sizeof buffer,fp) != NULL) result += buffer;
	fclose(fp);
	sec1 = stoi(result);
        if(month == 12) {
            sprintf(cmd,"date --date='%d%02d01' +%s",year+1,1,"%s");
            result = "";
            fp = popen(cmd,"r");
            while(fgets(buffer,sizeof buffer,fp) != NULL) result += buffer;
            fclose(fp);
            sec2 = stoi(result);
        } else {
            sprintf(cmd,"date --date='%d%02d01' +%s",year,month+1,"%s");
            result = "";
            fp = popen(cmd,"r");
            while(fgets(buffer,sizeof buffer,fp) != NULL) result += buffer;
            fclose(fp);
            sec2 = stoi(result);
        }
	startindex[year-startyear][month-1] = (sec1-startsec)/86400;
        endindex[year-startyear][month-1] = (sec2-startsec)/86400-1;
    }


    fp1 = fopen(argv[1],"rb");
    fp2 = fopen(argv[2],"rb");
    
    for(s=0;s<Sx;s++) {
    //for(s=0;s<1;s++) {
        printf("s = %d\n",s);
        recordno++;
	
        x = fread(cwt1[0][0],Nx*Sy*Ny*sizeof(complex<float>),1,fp1);                    // f cwt1

	if(strcmp(argv[1],argv[2]) == 0)
	    memcpy(cwt2[0][0],cwt1[0][0],Nx*Sy*Ny*sizeof(complex<float>));	
	else
	    x = fread(cwt2[0][0],Nx*Sy*Ny*sizeof(complex<float>),1,fp2);
        
        for(x=0;x<Nx;x++)
	for(k=0;k<Sy;k++)	
	for(y=0;y<Ny;y++)
            cov[x][k][y] = (cwt1[x][k][y]*cwt2[x][k][y].conjugate()).getreal();
	

        counter = 0;
        for(i=0;i<tiles;i++)                                                          // 0-27.5, 30-57.5, 60-87.5, .....
	for(x=i*Nx/tiles;x<(i+1)*Nx/tiles;x++)	
        for(k=0;k<Sy;k++)
        for(year=startyear;year<=endyear;year++) 
        for(month=1;month<=12;month++)
	for(y=startindex[year-startyear][month-1];y<=endindex[year-startyear][month-1];y++) {
	//	printf("year month i s x k y = %d %d %d %d %d %d %d\n",year,month,i,s,x,k,y);
            data[year-startyear][month-1][i][s][k] += cov[x][k][y];
	}
    }
    
    fclose(fp1);
    fclose(fp2);
    printf("records = %d\n",recordno);


    for(s=0;s<Sx;s++) {
        a = pow(2.,s/4.); 
        for(k=0;k<Sy;k++) {
            b = pow(2.,k/4.);
	    for(year=startyear;year<=endyear;year++)
	    for(month=1;month<=12;month++)
	    for(i=0;i<tiles;i++)
                data[year-startyear][month-1][i][s][k] /= (a*b); 
	}
    }

    a = 0.5*dx*dy/Ny/Nx*pow(0.25*log(2.)/pi/integral_covariance(param),2.);

    for(s=0;s<Sx;s++)
    for(k=0;k<Sy;k++)
    for(year=startyear;year<=endyear;year++)
    for(month=1;month<=12;month++)
    for(i=0;i<tiles;i++)
        data[year-startyear][month-1][i][s][k] *= a;
    

    
    fp = fopen("test.bin","wb"); 
    fwrite(data[0][0][0][0],(endyear-startyear+1)*12*tiles*Sx*Sy*sizeof(double),1,fp);
    fclose(fp);
}
