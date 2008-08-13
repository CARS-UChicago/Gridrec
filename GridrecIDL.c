/* File GridrecIDL.c
   This file is a thin wrapper layer which is called from IDL
   It calls gridrec functions
   Mark Rivers
   May, 2000

   Aug. 22, 2001 MLR  Modifications to allow compiling on non-Windows systems
*/

#ifdef WINDOWS
#define EXPORT __declspec(dllexport)a
#else
#define EXPORT
#endif

#include "grid.h"

static grid_struct GP;
static sg_struct SGP;

EXPORT void recon_init_IDL (int argc, char *argv[])
{
	long *imgsiz;
	SGP.n_ang    = *(long *) argv[0];
	SGP.n_det    = *(long *) argv[1];
	SGP.geom     = *(long *) argv[2];
	SGP.angles   =  (float *)argv[3];
	SGP.center   = *(float *)argv[4];
	get_pswf(*(float *)argv[5], &GP.pswf);
	GP.sampl     = *(float *)argv[6];
	GP.R         = *(float *)argv[7];
	GP.MaxPixSiz = *(float *)argv[8];
	GP.X0        = *(float *)argv[9];
	GP.Y0        = *(float *)argv[10];
	GP.ltbl      = *(long *) argv[11];
	GP.filter    =  get_filter((char *)argv[12]);
	imgsiz       =  (long *) argv[13];
	recon_init(&GP, &SGP, imgsiz);
	return;
}

EXPORT void do_recon_IDL (int argc, char *argv[])
{

	int n_ang  = *(long *)argv[0];
	int n_det  = *(long *)argv[1];
	int imgsiz = *(long *)argv[2];
	float **G1, **G2, **S1, **S2;
	int i;

	/* IDL passes addresses of arrays (float *), while Gridrec
	   wants a pointer to a table of the starting address of each row.
	   Need to build those tables */
	G1 = (float **) malloc((size_t) (n_ang * sizeof(float *)));
	G2 = (float **) malloc((size_t) (n_ang * sizeof(float *)));
    G1[0] = (float *) argv[3];
    G2[0] = (float *) argv[4];
    for (i = 1; i < n_ang; i++) {
		G1[i] = G1[i-1] + n_det;
		G2[i] = G2[i-1] + n_det;
    }

	S1 = (float **) malloc((size_t) (imgsiz * sizeof(float *)));
	S2 = (float **) malloc((size_t) (imgsiz * sizeof(float *)));
    S1[0] = (float *) argv[5];
    S2[0] = (float *) argv[6];
    for (i = 1; i < imgsiz; i++) {
		S1[i] = S1[i-1] + imgsiz;
		S2[i] = S2[i-1] + imgsiz;
    }

	do_recon(G1, G2, &S1, &S2);
	free(G1);
	free(G2);
	free(S1);
	free(S2);
	recon_free();
}







