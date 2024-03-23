#include "allocate.h"          
#include "morlet.h"
#include "omp.h"


// usage:   main.exe inputfile outputfile Ny Nx Sy Sx dy dx
// input file dimension  [Ny][Nx]
// output file dimension [Sx][Nx][Sy][aligned_dim(Ny)][2]
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
    complex<double> ***transform_real;  // transform of Re{xtransform[:][Sx][Nx]} in y direction    Dim: [Nx][Sy][Ny]
    complex<double> ***transform_imga;  // transform of Im{xtransform[:][Sx][Nx]} in y direction   


    FILE *fp,*fp1,*fp2;
    char cmd[500];
    __uint128_t filesize;               // linux gnu C++ has difficulty handling __uint128_t. not used 

    int Nthreads,thread; 
	
#if AVX512F > 0
    alignas(ALIGN) __m512d mw01_wavelet2d_a;
    alignas(ALIGN) __m512d mw01_wavelet2d_b;
#elif AVX > 0
    alignas(ALIGN) __m256d mw01_wavelet2d_a;
    alignas(ALIGN) __m256d mw01_wavelet2d_b;
#endif



    // set constants
    pi = acos(-1.);
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
    // Note that we have to read line by line, becoz the data array has actual size [Ny][aligned_dim(Nx)]. See allocate.h  
    fp = fopen(argv[1],"rb");
    for(y=0;y<Ny;y++) i = fread(data[y],sizeof(double),Nx,fp);
    fclose(fp);



    sprintf(cmd,"%s1",argv[2]);
    fp1 = fopen(cmd,"wb");
    sprintf(cmd,"%s2",argv[2]);
    fp2 = fopen(cmd,"wb");

    #pragma omp parallel for default(shared)
    for(y=0;y<Ny;y++)                                                           // xtransform for each y
        morlet<double>(data[y],xtransform[y],Nx,Sx,param,dx);


    for(s=0;s<Sx;s++) {                                                         // each Sx
        printf("s= %d\n",s);
        system("date");

        #pragma omp parallel for default(shared) private(thread,y,i,mw01_wavelet2d_a,mw01_wavelet2d_b)
	for(x=0;x<Nx;x++) {                                                 // each Nx
            thread = omp_get_thread_num();

            for(y=0;y<Ny;y++) {
                yseries1[thread][y] = xtransform[y][s][x].getreal();
                yseries2[thread][y] = xtransform[y][s][x].getimga();
	    }

            morlet<double>(yseries1[thread],transform_real[x],Ny,Sy,param,dy);
            morlet<double>(yseries2[thread],transform_imga[x],Ny,Sy,param,dy);

#if !defined(AVX) || AVX == 0
	    for(i=0;i<Sy;i++)
	    for(y=0;y<Ny;y++) {
                cwt1[x][i][y].setreal((float)(transform_real[x][i][y].getreal() - transform_imga[x][i][y].getimga()));
                cwt1[x][i][y].setimga((float)(transform_real[x][i][y].getimga() + transform_imga[x][i][y].getreal()));
                cwt2[x][i][y].setreal((float)(transform_real[x][i][y].getreal() + transform_imga[x][i][y].getimga()));
                cwt2[x][i][y].setimga((float)(transform_imga[x][i][y].getreal() - transform_real[x][i][y].getimga()));
            }
#elif AVX512F == 0
            for(i=0;i<Sy;i++)
            for(y=0;y<Ny;y+=2) {
                mw01_wavelet2d_a = _mm256_load_pd((double *)&transform_real[x][i][y]);    
                mw01_wavelet2d_b = _mm256_permute_pd(_mm256_load_pd((double *)&transform_imga[x][i][y]),0b0101);
                _mm_store_ps((float *)&cwt1[x][i][y],_mm256_cvtpd_ps(_mm256_addsub_pd(mw01_wavelet2d_b,mw01_wavelet2d_a)));
                _mm_store_ps((float *)&cwt2[x][i][y],_mm256_cvtpd_ps(_mm256_add_pd(_mm256_xor_pd(mw01_wavelet2d_a,_mm256_setr_pd(0.0,-0.0,0.0,-0.0)),mw01_wavelet2d_b)));
                //_mm_store_ps((float *)&cwt1[x][i][y],_mm256_cvtpd_ps(_mm256_addsub_pd(_mm256_permute_pd(_mm256_load_pd((double *)&transform_imga[x][i][y]),0b0101),_mm256_load_pd((double *)&transform_real[x][i][y]))));
                //_mm_store_ps((float *)&cwt2[x][i][y],_mm256_cvtpd_ps(_mm256_add_pd(_mm256_xor_pd(_mm256_load_pd((double *)&transform_real[x][i][y]),_mm256_setr_pd(0.0,-0.0,0.0,-0.0)),_mm256_permute_pd(_mm256_load_pd((double *)&transform_imga[x][i][y]),0b0101))));
            }
#else
            for(i=0;i<Sy;i++)
            for(y=0;y<Ny;y+=4) {
                mw01_wavelet2d_a = _mm512_load_pd((double *)&transform_real[x][i][y]);
                mw01_wavelet2d_b = _mm512_permute_pd(_mm512_load_pd((double *)&transform_imga[x][i][y]),0b01010101);
                _mm256_store_ps((float *)&cwt1[x][i][y],_mm512_cvtpd_ps(_mm512_add_pd(mw01_wavelet2d_a,_mm512_xor_pd(mw01_wavelet2d_b,_mm512_setr_pd(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0)))));
                _mm256_store_ps((float *)&cwt2[x][i][y],_mm512_cvtpd_ps(_mm512_add_pd(_mm512_xor_pd(mw01_wavelet2d_a,_mm512_setr_pd(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)),mw01_wavelet2d_b)));
                //_mm256_store_ps((float *)&cwt1[x][i][y],_mm512_cvtpd_ps(_mm512_add_pd(_mm512_load_pd((double *)&transform_real[x][i][y]),_mm512_xor_pd(_mm512_permute_pd(_mm512_load_pd((double *)&transform_imga[x][i][y]),0b01010101),_mm512_setr_pd(-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0)))));
                //_mm256_store_ps((float *)&cwt2[x][i][y],_mm512_cvtpd_ps(_mm512_add_pd(_mm512_xor_pd(_mm512_load_pd((double *)&transform_real[x][i][y]),_mm512_setr_pd(0.0,-0.0,0.0,-0.0,0.0,-0.0,0.0,-0.0)),_mm512_permute_pd(_mm512_load_pd((double *)&transform_imga[x][i][y]),0b01010101))));
            }
#endif
	}
        fwrite(cwt1[0][0],Nx*Sy*aligned_dim<complex<float>>(Ny)*sizeof(complex<float>),1,fp1);
        fwrite(cwt2[0][0],Nx*Sy*aligned_dim<complex<float>>(Ny)*sizeof(complex<float>),1,fp2);
        // Note the aligned_dim(Ny). These 1 line fwrite could save time than fwriting Nx*Sy times.
        // When reading these files, have to use same structure. See waveletcov_stat.cpp as an example.
    }

    fclose(fp1);
    fclose(fp2);

}
