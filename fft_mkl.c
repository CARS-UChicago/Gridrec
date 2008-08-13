/*** File fft_mkl.c 5/5/2000
     This file calls routines from the Intel Math Kernel Library.  It 
     emulates the 1-D and n-D routines from Numerical Recipes, so that 
     gridrec requires essentially no modification

     Written by Mark Rivers
 **/
#include <stdlib.h>

/* The symbols for Linux and Windows are different.  Note that in both cases
 * we are using the routines designed to be called from Fortran, i.e. with
 * COMPLEX data type.  C does not normally support COMPLEX, but Gridrec uses
 * a C structure to emulate it.
 */

#ifdef LINUX
#define CFFT1D cfft1d_
#define CFFT2D cfft2d_
#endif

void four1(float data[], unsigned long nn, int isign)
{
   static float* wsave;
   static int n_prev;
   float scale, *p;
   int n = nn;
   int i;
   int zero = 0;

   /* Call the Intel Math Kernel Library routine */
	if ((isign == 0) || (n != n_prev)) {
		/* The required storage is (3*N)/2 complex elements = 3*N floats */
		free(wsave);
		wsave = malloc(3 * n * sizeof(float));
		n_prev = n;
	    CFFT1D(data+1, &n, &zero, wsave);
	}
	/* The Numerical Recipes routines are passed a pointer to one element
	   before the start of the array - add one */
	CFFT1D(data+1, &n, &isign, wsave);
	/* Must rescale data if isign is 1, since Numerical Recipes output is scaled by N, CFFT1D is not */
	if (isign == 1) {
		scale = (float)n;
	    for (i=0, p=data+1; i<2*n; i++) *p++ *= scale;
	}
}
	
void fourn(float data[], unsigned long nn[], int ndim, int isign)
{

	int nx = nn[2];
	int ny = nn[1];
	int i;
	float scale, *p;

	/* NOTE: This function only works for ndim=2 */
	if (ndim != 2) return;
	/* Call Intel Math Kernel Library routine */
	/* The Numerical Recipes routines are passed a pointer to one element
	   before the start of the array - add one */
	CFFT2D(data+1, &nx, &ny, &isign);
	/* Must rescale data if isign is 1, since Numerical Recipes output is scaled by N, CFFT2D is not */
	if (isign == 1) {
		scale = (float)(nx*ny);
	    for (i=0, p=data+1; i<2*nx*ny; i++) *p++ *= scale;
	}
}
