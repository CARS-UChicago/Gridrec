/*** File grid.h  12:27 PM 11/7/97 **/

#define ANSI

/**** System includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stddef.h>
#include <time.h>
#include <sys/stat.h>


/**** Macros and typedefs ****/
#define max(A,B) ((A)>(B)?(A):(B))
#define min(A,B) ((A)<(B)?(A):(B))
#define free_matrix(A) (free(*(A)),free(A))
#define abs(A) ((A)>0 ?(A):-(A))
#define pi  3.14159265359
#define TOLERANCE 0.1	/* For comparing centers of two sinugrams */
#define LTBL_DEF 512	/* Default lookup table length */

typedef struct {float r;float i;} complex;

typedef struct SG_STRUCT { 		/* Sinugram parameters */

   int n_ang;	/* No. angles in sinugram */
   int n_det;	/* No. elems (detectors) per angle */
   int geom;		/* 0 if array of angles provided;
				1,2 if uniform in half,full circle */ 
   float *angles;	/* Ptr to the array of angles, if used */
   float center;	/* Rotation axis location */
			
			} sg_struct;

typedef struct PSWF_STRUCT{ /*Prolate spheroidal wave fcn (PSWF) data */

   float C;	/* Parameter for particular 0th order pswf being used*/
   int nt; /* Degree of Legendre polynomial expansion */
   float lmbda; 	/* Eigenvalue */
   float coefs[15];	/* Coeffs for Legendre polynomial expansion */
	
			} pswf_struct;

typedef struct GRID_STRUCT { /* Parameters for gridding algorithm */

   pswf_struct *pswf;	/* Ptr to data for PSWF being used  */
   float sampl;	  	/* "Oversampling" ratio */
   float MaxPixSiz; 	/* Max pixel size for reconstruction */
   float R;		/* Region of interest (ROI) relative size */
   float X0;		/* (X0,Y0)=Offset of ROI from rotation axis, */
   float Y0;		/* in units of center-to-edge distance.  */
   char fname[16];		/*  Name of filter function   */		
   float (*filter)(float);	/* Ptr to filter function.  */
   long ltbl;		/* No. elements in convolvent lookup tables. */

			} grid_struct;

/** Global variables **/
extern int verbose;



/**** Function Prototypes ****/

/** Defined in grid_io.c **/
void get_parm(int argc,char **argv, grid_struct *A);
int data_setup(int argc, char *argv[]);
float **get_sgram(int ifile, sg_struct *A);
void rel_sgram(float **S);
void put_image(int ifile,float **image, int size);
void usage (void);
float *malloc_vector_f (long n);
complex *malloc_vector_c (long n) ;
float **malloc_matrix_f (long nr, long nc);
complex **malloc_matrix_c (long nr, long nc);
void prn_gparm(grid_struct *G);
void prn_pswf(pswf_struct *P);
void prn_sgparm(sg_struct *S);

/** Defined in pswf.c **/
void get_pswf(float C, pswf_struct **P);

/** Defined in filters.c  **/
float (*get_filter(char *name))(float);

/** Defined in grid.c  **/
void recon_init(grid_struct *gparm,sg_struct *sgparm);
void do_recon(float **sgram1, float **sgram2, float ***image1,
		float ***image2, long *im_size);
void recon_free(void);

/** Defined in fft.c (FFT routines from Numerical Recipes) **/
void four1(float data[], unsigned long nn, int isign);
void fourn(float data[], unsigned long nn[], int ndim, int isign);

