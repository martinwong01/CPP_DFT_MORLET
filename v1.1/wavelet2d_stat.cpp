#include "allocate.h"
#include "morlet.h"


int main(int argc, char *argv[]) {
    int Nx;                           
    int Ny;
    int x,y,i,j,k;                         // index for position
    double dx;
    double dy;

    int Sx;                                // number of scales in x  
    int Sy;
    int s;                                 // index for scale

    int counter;

    double a;
    double b;
    double pi;
    double param;

    complex<float> ***cwt1;                // wavelet transform data
    complex<float> ***cwt2;
    double **scales_sqroot_inv;

    double **data;                         // 2D data
    double areamean;
    double variance;
    double *xmean;
    double *ymean;

    double **anomaly;                      // Equation (B1), LHS
    double **anom1;                        // Equation (B1), RHS wavefunc1, summed over all J/K
    double **anom2;                        // Equation (B1), RHS wavefunc2, summed over all J/K
    double MSE;                            // mean squared error between LHS and RHS in (B1) 
   
    double LHS1,LHS2,LHS3,LHS4,LHS;        // Equation (B2), LHS terms 
    double **cov1;                         // Equation (B2), RHS wavefunc1, summed over all x/y 
    double **cov2;                         // Equation (B2), RHS wavefunc2, summed over all x/y
    double RHS;

    FILE *fp1,*fp2;


    // set constants
    pi = acos(-1.);
    param = 16.*pi/9.-9./32./pi;


    Ny = 16071;
    Nx = 144;
    Sy = 35;
    Sx = 33;
    dy = 1.;
    dx = 1.;


    cwt1 = allocate3D<complex<float>>(Nx,Sy,Ny);
    cwt2 = allocate3D<complex<float>>(Nx,Sy,Ny);
    scales_sqroot_inv = allocate2D<double>(Sx,Sy);
    data = allocate2D<double>(Ny,Nx);
    xmean = allocate1D<double>(Ny);
    ymean = allocate1D<double>(Nx);
    anomaly = allocate2D<double>(Ny,Nx);
    anom1 = allocate2D<double>(Ny,Nx);
    anom2 = allocate2D<double>(Ny,Nx);
    cov1 = allocate2D<double>(Sx,Sy);
    cov2 = allocate2D<double>(Sx,Sy);



    for(s=0;s<Sx;s++)
    for(k=0;k<Sy;k++)
        scales_sqroot_inv[s][k] = 1./sqrt(pow(2.,s/4.)*pow(2.,k/4.)); 

    fp1 = fopen("olr_sym_0","rb");
    for(y=0;y<Ny;y++) i = fread(data[y],sizeof(double),Nx,fp1);
    fclose(fp1);

    for(x=0;x<Nx;x++) {
        for(y=0;y<Ny;y++) ymean[x] += data[y][x];
        ymean[x] = ymean[x]/Ny;
    }

    for(y=0;y<Ny;y++) {
        for(x=0;x<Nx;x++) xmean[y] += data[y][x];
        xmean[y] = xmean[y]/Nx;
    }

    for(y=0;y<Ny;y++)
    for(x=0;x<Nx;x++)
        areamean += data[y][x];
    areamean = areamean/Nx/Ny;

    variance = 0.;
    for(y=0;y<Ny;y++)
    for(x=0;x<Nx;x++)
        variance += (data[y][x] - areamean)*(data[y][x] - areamean);
    variance = variance/Ny/Nx;

    printf("data2D area_mean,variance wrt area_mean = %lf,%lf\n",areamean,variance);

    


    fp1 = fopen("outputfile1","rb");
    fp2 = fopen("outputfile2","rb");
    for(s=0;s<Sx;s++) {
        //printf("s = %d\n",s);
        x = fread(cwt1[0][0],Nx*Sy*aligned_dim<complex<float>>(Ny)*sizeof(complex<float>),1,fp1); // Note aligned_dim here
        x = fread(cwt2[0][0],Nx*Sy*aligned_dim<complex<float>>(Ny)*sizeof(complex<float>),1,fp2); // Note aligned_dim here
        for(x=0;x<Nx;x++)
        for(k=0;k<Sy;k++)
        for(y=0;y<Ny;y++) {
            anom1[y][x] += cwt1[x][k][y].getreal()*scales_sqroot_inv[s][k];
            cov1[s][k] += (cwt1[x][k][y]*cwt1[x][k][y].conjugate()).getreal();
            anom2[y][x] += cwt2[x][k][y].getreal()*scales_sqroot_inv[s][k];
            cov2[s][k] += (cwt2[x][k][y]*cwt2[x][k][y].conjugate()).getreal();
        }
    }
    fclose(fp1);
    fclose(fp2);

    a = 0.5*dx*dy/Ny/Nx*pow(0.25*log(2.)/pi/integral_covariance(param),2.);
    for(s=0;s<Sx;s++)
    for(k=0;k<Sy;k++) {
        cov1[s][k] *= scales_sqroot_inv[s][k]*scales_sqroot_inv[s][k]*a;
        cov2[s][k] *= scales_sqroot_inv[s][k]*scales_sqroot_inv[s][k]*a;
    }
   
    a = sqrt(dx*dy)*pow(0.25*log(2.)/integral_reconstruction(param),2.)/pi;
    for(y=0;y<Ny;y++)
    for(x=0;x<Nx;x++) {
        anom1[y][x] *= a;
        anom2[y][x] *= a;
    }
        

   
//  Results 
    for(y=0;y<Ny;y++)
    for(x=0;x<Nx;x++)
        anomaly[y][x] = data[y][x] - xmean[y] - ymean[x] + areamean;     // Eq (B1) LHS

    a = 0.;
    for(y=0;y<Ny;y++)
    for(x=0;x<Nx;x++)
        a += anomaly[y][x];
    a = a/Ny/Nx;
    variance = 0.;
    for(y=0;y<Ny;y++)
    for(x=0;x<Nx;x++)
        variance += (anomaly[y][x]-a)*(anomaly[y][x]-a);
    variance = variance/Ny/Nx;
    printf("Mean/Variance of (B1) LHS = %lf %lf\n",a,variance);


    MSE = 0.;
    for(y=0;y<Ny;y++)
    for(x=0;x<Nx;x++)
        MSE += (anomaly[y][x] - anom1[y][x] - anom2[y][x])*(anomaly[y][x] - anom1[y][x] - anom2[y][x]);
    MSE = MSE/Ny/Nx;
    printf("MSE of reconstruction = %lf\n",MSE);



    LHS1 = 0.;
    LHS2 = 0.;
    LHS3 = 0.;
    LHS4 = 0.;

    for(x=0;x<Nx;x++)
    for(y=0;y<Ny;y++)
        LHS1 += data[y][x]*data[y][x];
    LHS1 = LHS1/Nx/Ny;

    for(y=0;y<Ny;y++)
    for(x=0;x<Nx;x++)
        LHS2 += data[y][x]*xmean[y];
    LHS2 = LHS2/Nx/Ny;

    for(y=0;y<Ny;y++)
    for(x=0;x<Nx;x++)
        LHS3 += data[y][x]*ymean[x];
    LHS3 = LHS3/Nx/Ny;

    LHS4 = areamean*areamean;

    LHS = LHS1 - LHS2 - LHS3 + LHS4;


    RHS = 0.;
    for(s=0;s<Sx;s++)
    for(k=0;k<Sy;k++)
        RHS += (cov1[s][k] + cov2[s][k]);

    printf("(B2) LHS RHS = %lf %lf\n",LHS,RHS);
}
