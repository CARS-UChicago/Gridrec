/* gcc fftw_test.c timer.c -I ~/lib/fftw-3.1.2/include/ -L ~/lib/fftw-3.1.2/lib/ -lfftw3f -lm -o a.serial */
/* gcc -DTHREADED fftw_test.c timer.c -I ~/lib/fftw-3.1.2/include/ -L ~/lib/fftw-3.1.2/lib/ -lfftw3f_threads -lfftw3f -lpthread -lm -o a.threads */
/* nvcc -DCUDA -I /usr/local/cuda/include -I ~/NVIDIA_CUDA_SDK/common/inc -L /usr/local/cuda/lib fftw_test.cu timer.cu -lcufft -o a.cuda */

// CUDA - use cuda, otherwise fftw
// POW2 - use power of two transform, otherwise use non power of 2
// ONED | TWOD | THREED - exclusive, dimension of transform
// THREADED - (!CUDA) - use threaded fftw
// HOSTMEM - (CUDA) memory copy from host included in timing
// INPLACE - transforms are done in place
// FASTPLAN - use FFTW_ESTIMATE for FFTW plan creation rather than FFTW_MEASURE

     #include <stdlib.h>
     #include <stdio.h>
     #include <math.h>

#ifdef CUDA
     #include <cufft.h>
     #include <cutil.h>
#else
     #include <fftw3.h>
     const int N_THREADS = 4;
#endif
     double dt( char );

#ifdef ONED
#ifndef POW2
#define NNPOW2 22
     int npow2[NNPOW2] = {6,9,12,15,18,24,36,80,108,210,504,1000,1960,4725,10368,27000,75600,165375,362880,1562500,3211264,6250000}; //,12250000,25401600}; 
#endif
#endif

#ifdef TWOD
#ifdef POW2
#define NPOW2 18
     int pow2[NPOW2][2] = { {4,4},{8,4},{4,8},{8,8},{16,16},{32,32},{64,64},{16,512},{128,64},{128,128},{256,128},{512,64},{64,1024},{256,256},{512,512},{1024,1024},{2048,2048},{4096,4096}}; //,{8192,8192} };  -- this fails when creating plan
//int igrid[2][3] = { {0, 1, 2}, {3, 4, 5} };
#else
#define NNPOW2 33
     int npow2[NNPOW2][2] = { {5,5},{6,6},{7,7},{9,9},{10,10},{11,11},{12,12},{13,13},{14,14},{15,15},{25,24},{48,48},{49,49},{60,60},{72,56},{75,75},{80,80},{84,84},{96,96},{100,100},{105,105},{112,112},{120,120},{144,144},{180,180},{240,240},{360,360},{1000,1000},{1050,1050},{1458,1458},{1960,1960},{2916,2916},{4116,4116} }; // ,{5832,5832},{8400,8400},{10368,10368} };
#endif
#endif

#ifdef THREED
#ifdef POW2
#define NPOW2 13
     int pow2[NPOW2][3] = { {4,4,4},{8,8,8},{4,8,16},{16,16,16},{32,32,32},{64,64,64},{256,64,32},{16,1024,64},{128,128,128},{512,128,64},{256,128,256},{256,256,256},{512,64,1024} }; // this won't fit on card ,{512,512,512} }
#else
#define NNPOW2 28
     int npow2[NNPOW2][3] = { {5,5,5},{6,6,6},{7,7,7},{9,9,9},{10,10,10},{11,11,11},{12,12,12},{13,13,13},{14,14,14},{15,15,15},{24,25,28},{48,48,48},{49,49,49},{60,60,60},{72,60,56},{75,75,75},{80,80,80},{84,84,84},{96,96,96},{100,100,100},{105,105,105},{112,112,112},{120,120,120},{144,144,144},{180,180,180},{210,210,210},{270,270,270},{324,324,324} }; // this won't fit on card ,{420,420,420} }
#endif
#endif

     int main( int argc, char** argv)
     {

#ifdef CUDA
         cufftHandle plan;
         cufftComplex *data, *datao, *devdata, *devdatao;
#else
         fftwf_complex *in, *out;
         fftwf_plan p;
#endif
         int i,k,n_el,n_test,nx;
#ifdef TWOD
         int ny;
#endif
#ifdef THREED
         int ny,nz;
#endif
         double sec, mflops;

#ifdef CUDA

       CUT_DEVICE_INIT();
#endif

#ifdef ONED
#ifdef POW2
       for (k=2; k<24; k++) 
       {
         nx = pow(2.0,(double)k);
         n_el = nx;
#else
       for (k=0; k<NNPOW2; k++)
       {
         nx = npow2[k];
         n_el = nx;
#endif
#endif

#ifdef TWOD 
#ifdef POW2
       for (k=0; k<NPOW2; k++) 
       {
         nx = pow2[k][0];
         ny = pow2[k][1];
         n_el = nx*ny; 
#else
       for (k=0; k<NNPOW2; k++)
       {
         nx = npow2[k][0];
         ny = npow2[k][1];
         n_el = nx*ny; 
#endif
#endif

#ifdef THREED 
#ifdef POW2
       for (k=0; k<NPOW2; k++) 
       {
         nx = pow2[k][0];
         ny = pow2[k][1];
         nz = pow2[k][2];
         n_el = nx*ny*nz; 
#else
       for (k=0; k<NNPOW2; k++)
       {
         nx = npow2[k][0];
         ny = npow2[k][1];
         nz = npow2[k][2];
         n_el = nx*ny*nz; 
#endif
#endif

#ifdef CUDA
         size_t arraySize = sizeof(cufftComplex) * n_el;
         cudaMallocHost((void**) &data, arraySize);
         cudaMallocHost((void**) &datao, arraySize);
#else
         in = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n_el);
         out = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * n_el);
#endif

#ifdef THREADED
         fftwf_init_threads();
         fftwf_plan_with_nthreads(N_THREADS);
#endif

// start plan creation
#ifdef CUDA

#ifdef ONED
         cufftPlan1d(&plan, n_el, CUFFT_C2C, 1);
#endif
#ifdef TWOD
         cufftPlan2d(&plan, nx, ny, CUFFT_C2C);
#endif
#ifdef THREED
         cufftPlan3d(&plan, nx, ny, nz, CUFFT_C2C);
#endif

#else

#ifdef ONED
#ifdef INPLACE
#ifdef FASTPLAN 
         p = fftwf_plan_dft_1d(n_el, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
#else
         p = fftwf_plan_dft_1d(n_el, in, in, FFTW_FORWARD, FFTW_MEASURE);
#endif
#else
#ifdef FASTPLAN 
         p = fftwf_plan_dft_1d(n_el, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
#else
         p = fftwf_plan_dft_1d(n_el, in, out, FFTW_FORWARD, FFTW_MEASURE);
#endif
#endif
#endif

#ifdef TWOD
#ifdef INPLACE
#ifdef FASTPLAN 
         p = fftwf_plan_dft_2d(nx, ny, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
#else
         p = fftwf_plan_dft_2d(nx, ny, in, in, FFTW_FORWARD, FFTW_MEASURE);
#endif
#else
#ifdef FASTPLAN 
         p = fftwf_plan_dft_2d(nx, ny, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
#else
         p = fftwf_plan_dft_2d(nx, ny, in, out, FFTW_FORWARD, FFTW_MEASURE);
#endif
#endif
#endif

#ifdef THREED
#ifdef INPLACE
#ifdef FASTPLAN 
         p = fftwf_plan_dft_3d(nx, ny, nz, in, in, FFTW_FORWARD, FFTW_ESTIMATE);
#else
         p = fftwf_plan_dft_3d(nx, ny, nz, in, in, FFTW_FORWARD, FFTW_MEASURE);
#endif
#else
#ifdef FASTPLAN 
         p = fftwf_plan_dft_3d(nx, ny, nz, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
#else
         p = fftwf_plan_dft_3d(nx, ny, nz, in, out, FFTW_FORWARD, FFTW_MEASURE);
#endif
#endif
#endif


#endif
// end plan creation

//malloc arrays on device and populate input

#ifdef CUDA
         cudaMalloc((void**)&devdata, arraySize);
         cudaMalloc((void**)&devdatao, arraySize);
#endif

         for (i=0; i<n_el; i++)
         {
#ifdef CUDA
//           data[i][0]=(float)i;
//           data[i][1]=(float)2*i;
           data[i].x=(float)i;
           data[i].y=(float)2*i;
#else
           in[i][0]=(float)i;
           in[i][1]=(float)2*i;
#endif
         }

// start timer, transfer data to card

#ifdef CUDA
#ifndef HOSTMEM
         cudaMemcpy(devdata, data, arraySize, cudaMemcpyHostToDevice);
#endif
#endif
         dt('i');

//#ifdef HOSTMEM
//         dt('i');
//
//#ifdef CUDA
//         cudaMemcpy(devdata, data, arraySize, cudaMemcpyHostToDevice);
//#endif
//
//#else
//
//#ifdef CUDA
//         cudaMemcpy(devdata, data, arraySize, cudaMemcpyHostToDevice);
//#endif
//         dt('i');
//#endif


// calculate ffts

         for (n_test=0; n_test<N_FFTS; n_test++)
         {
#ifdef CUDA
#ifdef HOSTMEM
           cudaMemcpy(devdata, data, arraySize, cudaMemcpyHostToDevice);
#endif
#ifdef INPLACE
           cufftExecC2C(plan, devdata, devdata, CUFFT_FORWARD);
#else
           cufftExecC2C(plan, devdata, devdatao, CUFFT_FORWARD);
#endif
#ifdef HOSTMEM
#ifdef INPLACE
         cudaMemcpy(datao, devdata, arraySize, cudaMemcpyDeviceToHost);
#else
         cudaMemcpy(datao, devdatao, arraySize, cudaMemcpyDeviceToHost);
#endif
#endif
#else
//           for (i=0; i<N_FFTS ; i++)
 //          {
             fftwf_execute(p); /* repeat as needed */
 //          }
#endif
         }

// copy data from card, end timer

//#ifdef HOSTMEM
//
//#ifdef CUDA
//#ifdef INPLACE
//         cudaMemcpy(datao, devdata, arraySize, cudaMemcpyDeviceToHost);
//#else
//         cudaMemcpy(datao, devdatao, arraySize, cudaMemcpyDeviceToHost);
//#endif
//#endif
//         sec=dt('e');
//#else
//         sec=dt('e');
//#ifdef CUDA
//#ifdef INPLACE
//         cudaMemcpy(datao, devdata, arraySize, cudaMemcpyDeviceToHost);
//#else
//         cudaMemcpy(datao, devdatao, arraySize, cudaMemcpyDeviceToHost);
//#endif
//#endif
//#endif

         sec=dt('e');

#ifdef CUDA
#ifndef HOSTMEM
#ifdef INPLACE
         cudaMemcpy(datao, devdata, arraySize, cudaMemcpyDeviceToHost);
#else
         cudaMemcpy(datao, devdatao, arraySize, cudaMemcpyDeviceToHost);
#endif
#endif
#endif


//HERE
// Have to change for 2d and 3d arrays

         mflops = 5. * (double) n_el * log( (double) n_el ) / log( 2.0 ) /  sec * (double) N_FFTS / 1e+6;  

#ifdef ONED
      //   printf("calculated %d ffts in %lf seconds = %lf mflops \n", N_FFTS, sec, mflops);
         printf("%ld\t%lf\t%lf\n", n_el, sec, mflops);
#endif
#ifdef TWOD
         printf("%ldx%ld\t%lf\t%lf\n", nx, ny, sec, mflops);
#endif
#ifdef THREED
         printf("%ldx%ldx%ld\t%lf\t%lf\n", nx, ny, nz, sec, mflops);
#endif


#ifdef THREADED
         fftwf_cleanup_threads();
#endif
#ifdef CUDA
         cufftDestroy(plan);
         cudaFreeHost(data);
         cudaFreeHost(datao);
         cudaFree(devdata);
         cudaFree(devdatao);
#else
         fftwf_destroy_plan(p);
         fftwf_free(in); fftwf_free(out);
#endif
       }
     }
