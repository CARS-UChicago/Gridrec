#include <stdlib.h>
#include <fftw3.h>

#ifdef LINUX
#define CFFT1D cfft1d_
#define CFFT2D cfft2d_
#endif

#ifdef WINDOWS
#define EXPORT __declspec(dllexport)a
#else
#define EXPORT
#endif

/** Defined in fft.c (FFT routines from Numerical Recipes) **/
void four1(float data[], unsigned long nn, int isign);
void fourn(float data[], unsigned long nn[], int ndim, int isign);

/* Storage for CFFT1D */
static float *wsave;

EXPORT void fft_test1n (int argc, char *argv[])
{

        float *data = (float *)argv[0];
        unsigned long nn  = *(long *)argv[1];
        int isign  = *(int *)argv[2];

        /* Call the Numerical Recipes routine */
        four1(data-1, nn, isign);
}

EXPORT void fft_test1i (int argc, char *argv[])
{

        float *data = (float *)argv[0];
        unsigned long nn  = *(long *)argv[1];
        int isign  = *(int *)argv[2];

        /* Call the Intel Math Kernel Library routine */
        if (isign == 0) {
                /* The required storage is (3*N)/2 complex elements = 3*N floats */
                free(wsave);
                wsave = malloc(3 * nn * sizeof(float));
        }
/*        CFFT1D(data, &nn, &isign, wsave); */
}

EXPORT void fft_test1f (int argc, char *argv[])
{

        float *data = (float *)argv[0];
        int n = *(int *)argv[1];
        int isign  = *(int *)argv[2];

        static int n_prev;
        static fftwf_complex *in, *out;
        static fftwf_plan forward_plan, backward_plan;

        if (n != n_prev) {
           /* Create plans */
           if (n_prev != 0) fftwf_free(in);
           in = fftwf_malloc(sizeof(fftwf_complex)*n);
           out = in;
           printf("fft_test1f: creating plans, n=%d, n_prev=%d\n", n, n_prev);
           n_prev = n;
           forward_plan = fftwf_plan_dft_1d(n, in, out, FFTW_FORWARD, FFTW_MEASURE);
           backward_plan = fftwf_plan_dft_1d(n, in, out, FFTW_BACKWARD, FFTW_MEASURE);
        }
        memcpy(in, data, n*sizeof(fftwf_complex));
        if (isign == -1) fftwf_execute(forward_plan);
        else             fftwf_execute(backward_plan);
        memcpy(data, in, n*sizeof(fftwf_complex));
}


EXPORT void fft_test2n (int argc, char *argv[])
{

        float *data = (float *)argv[0];
        unsigned long nx  = *(long *)argv[1];
        unsigned long ny  = *(long *)argv[2];
        int isign  = *(int *)argv[3];
        unsigned long nn[2];

        nn[0] = ny;
        nn[1] = nx;

        /* Call the Numerical Recipes routine */
        fourn(data-1, nn-1, 2, isign);
}

EXPORT void fft_test2i (int argc, char *argv[])
{

        float *data = (float *)argv[0];
        int nx  = *(int *)argv[1];
        int ny  = *(int *)argv[2];
        int isign  = *(int *)argv[3];

        /* Call Intel Math Kernel Library routine */
/*        CFFT2D(data, &nx, &ny, &isign); */
}

EXPORT void fft_test2f (int argc, char *argv[])
{

        float *data = (float *)argv[0];
        int nx  = *(int *)argv[1];
        int ny  = *(int *)argv[2];
        int isign  = *(int *)argv[3];

        static int nx_prev, ny_prev;
        static fftwf_complex *in, *out;
        static fftwf_plan forward_plan, backward_plan;

        if ((nx != nx_prev) || (ny != ny_prev)) {
           /* Create plans */
           if (nx_prev != 0) fftwf_free(in);
           in = fftwf_malloc(sizeof(fftwf_complex)*nx*ny);
           out = in;
           printf("fft_test2f: creating plans, nx=%d, ny=%d, nx_prev=%d, ny_prev=%d\n", 
                  nx, ny, nx_prev, ny_prev);
           nx_prev = nx;
           ny_prev = ny;
           forward_plan = fftwf_plan_dft_2d(ny, nx, in, out, FFTW_FORWARD, FFTW_MEASURE);
           backward_plan = fftwf_plan_dft_2d(ny, nx, in, out, FFTW_BACKWARD, FFTW_MEASURE);
        }
        memcpy(in, data, nx*ny*sizeof(fftwf_complex));
        if (isign == -1) fftwf_execute(forward_plan);
        else             fftwf_execute(backward_plan);
        memcpy(data, in, nx*ny*sizeof(fftwf_complex));
}

