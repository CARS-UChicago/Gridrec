 /* 
 * Read the 2D netCDF files in and write a 3D Volume netCDF file.
 */

/* System includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>
/* Header file for netCDF */
#include "netcdf.h"

/* Function prototypes */

int
fexist (char *filename);

float **
malloc_matrix_f (long nr, long nc);

void
free_matrix_f (float **m);

float ***
malloc_tensor_f (long nr, long nc, long nd);

void
free_tensor_f (float ***t);

unsigned short ***
malloc_tensor_us (long nr, long nc, long nd);

void
free_tensor_us (unsigned short ***t);

void
writenetCDF (int y, int x, int z);

float **
readnetCDF (char *name, int *xsize, int *ysize);

void
usage (void);

typedef enum TrueFalse { False, True } Flag;
Flag verbose, usint, subvol;

#define SCALE 1.0e+6

unsigned short ***Vol3D = NULL;
float ***Data3D = NULL;

int
main (int argc, char *argv[]) {
    char CDFfile[256];
    float **Data2D = NULL;
    long files = 0;
    int xsize = 0, ysize = 0,zsize = 0;
    int right = 0, left = 0, top = 0, bottom = 0;
    int c, x, y, z;
    extern char *optarg;
    extern int optind, opterr;
    int optindr;
    int errflg = 0;
    float min, max;

    while ((c = getopt(argc, argv, "vusl:r:t:b:")) != EOF)
        switch (c) {
        case 'v':
            verbose = True;
          break;
	case 'u':
	    usint = True;
	    break;
	case 's':
	    subvol = True;
	    break;
	case 'l':
	    left = atoi(optarg);
	    break;
	case 'r':
	    right = atoi(optarg);
	    break;
	case 't':
	    top = atoi(optarg);
	    break;
	case 'b':
	    bottom = atoi(optarg);
	    break;
        case '?':
            errflg++;
        }

    if (errflg) {
        usage ();
    } else {
        printf("\t-v\tverbose = %d\n", verbose);
        printf("\n");
    }

    zsize = files = argc - optind;
    if (files == 0) {
        fprintf (stderr, "No filename was specified.\n");
        exit(1);
    } else {
        strcpy (CDFfile, argv[optind]);
    /* Read the file names from the command line; make sure
       that they all exist. */

   if (verbose) printf("Attempting to read %d files.\n", files);
    optindr = optind;
    for (; optind < argc; optind++) {
        if (!fexist(argv[optind])) {
            fprintf (stderr, "File %s does not exist. Aborting execution.\n",
                     argv[optind]);
            exit(2);
        }
    }

    /* All the files exist, now process them */
    optind = optindr;
    for (z = 0; optind < argc; optind++, z++) {
        strcpy (CDFfile, argv[optind]);

        /* Read the netCDF file */
     if (verbose) printf ("Opening slice data file %s. \n", CDFfile);

     Data2D = readnetCDF(CDFfile, &xsize, &ysize);

     if (usint){

        min = 0.0;
        max = 0.0;

        for (x = 0;  x < xsize;  x++ ) {
              for (y = 0;  y < ysize;  y++) {
		  if (Data2D[x][y] < min)
		  min = Data2D[x][y]; 
		  else if (Data2D[x][y] > max)
		  max = Data2D[x][y]; 
              }
        }
   
 	/* Shift negative values, and rescale */

        for (x = 0;  x < xsize;  x++ ) {
               for (y = 0;  y < ysize;  y++) {
/*      	Data2D[x][y]  = max * ((Data2D[x][y]-min)/(max-min));*/
		Data2D[x][y] = (Data2D[x][y] - min) * SCALE;
               }
        }

     	if (subvol){
           xsize -= (left + right);
           ysize -= (top + bottom);
        }

       /* Allocate the 3D usint volume */
       if (xsize && ysize && zsize && z == 0) {
            Vol3D = malloc_tensor_us (xsize, ysize, zsize);
            if (!Vol3D) {
                fprintf (stderr, "Not enough space to allocate 3D array.\n");
                exit(3);
            }
        }

	for (x = 0; x < xsize; x++){
           for (y = 0; y < ysize; y++){ 
                Vol3D[x][y][z] = Data2D[x + left][y + bottom];
           }
        }
       } else {

     	if (subvol){
           xsize -= (left + right);
           ysize -= (top + bottom);
     	}

       /* Allocate the 3D float volume */
       if (xsize && ysize && zsize && z == 0) {
            Data3D = malloc_tensor_f (xsize, ysize, zsize);
            if (!Data3D) {
                fprintf (stderr, "Not enough space to allocate 3D array.\n");
                exit(3);
            }
        }
	for (x = 0; x < xsize; x++){
           for (y = 0; y < ysize; y++){ 
                Data3D[x][y][z] = Data2D[x + left][y + bottom];
           }
        }
       } /*end else*/

    } /* end of for z loop*/

      free_matrix_f(Data2D);

    } /*end of if (files = 0) else */

        /* Write the slices to a netCDF file */
        if (verbose) printf ("Writing Volume to netCDF file...\n");

    if (usint){
        writenetCDF (ysize, xsize, zsize);
        free_tensor_us(Vol3D);
    } else {
        writenetCDF (ysize, xsize, zsize);
        free_tensor_f(Data3D);
    }
     
   	exit(0);
}


/*
  **************** end of main ****************
*/


float **
readnetCDF (char *name, int *xsize, int *ysize)

{
    char CDFname[256];
    int ncid, xid, yid;
    long xdim, ydim;
    int slice_dimids[2], slice_id;
    long start[2], count[2];
    int i,j,cnt;
    float **Data = NULL;

    strcpy (CDFname, name);

    /* Open a netCDF file for reading */
    ncid = ncopen (CDFname, NC_NOWRITE);
    /*Get the dimensions of the Hyperslab*/
    xid = ncdimid (ncid, "x");
    ncdiminq (ncid, xid, (char *) 0, &xdim);
    *xsize = xdim;
    yid = ncdimid (ncid, "y");
    ncdiminq (ncid, yid, (char *) 0, &ydim);
    *ysize = ydim;
    /*Get the variable id*/
    slice_dimids[0] = xid;
    slice_dimids[1] = yid;
    slice_id = ncvarid (ncid, "slice");

    /*Allocate space for 2D array*/
    if (xdim && ydim){
    Data = malloc_matrix_f (xdim,ydim);

   /* Read the data */
    start[0] = 0;
    start[1] = 0;
    count[0] = xdim;
    count[1] = ydim;
    ncvarget (ncid, slice_id, start, count, (void *)&Data[start[0]][start[1]]);

    /* Close the netCDF file */
        ncclose (ncid);
    }

    return Data;
}

void
writenetCDF (int y, int x, int z)
{
    int ncid;
    int xid, yid, zid;
    int volume_dimids[3], volume_id;
    long start[3], count[3];

   /* Create a netCDF file for writing */
    ncid = nccreate ("volume.cdf",  NC_CLOBBER);

    /* Define the dimensions */
    xid = ncdimdef (ncid, "x", x);
    yid = ncdimdef (ncid, "y", y);
    zid = ncdimdef (ncid, "z", z);
    volume_dimids[0] = xid;
    volume_dimids[1] = yid;
    volume_dimids[2] = zid;

    if (usint){

    volume_id = ncvardef (ncid, "volume", NC_SHORT, 3, volume_dimids);
            /* Leave defined mode */
    ncendef (ncid);

    /* Write the data */
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = x;
    count[1] = y;
    count[2] = z;
    ncvarput (ncid, volume_id, start, count, (void **)&Vol3D[start[0]][start[1]][start[2]]);

    } else {

    volume_id = ncvardef (ncid, "volume", NC_FLOAT, 3, volume_dimids);
            /* Leave defined mode */
    ncendef (ncid);

    /* Write the data */
    start[0] = 0;
    start[1] = 0;
    start[2] = 0;
    count[0] = x;
    count[1] = y;
    count[2] = z;
    ncvarput (ncid, volume_id, start, count, (void **)&Data3D[start[0]][start[1]][start[2]]);
    }
        
    /* Close the netCDF file */
    ncclose (ncid);
    return;
}

void
usage (void)
{
    fprintf(stderr, "usage: VOLUME [-v] files...\n");
    exit (2);
}

int
fexist (char *filename) 
{
    /* Check for existence of a file */
    struct stat stbuf;
    return (stat (filename, &stbuf) == 0);
}

/* Memory allocation routines */


float **
malloc_matrix_f (long nr, long nc)
{
    float **m = NULL;
    long i;

    /* Allocate pointers to rows */
    m = (float **) malloc((size_t) (nr * sizeof(float *)));
    if (!m) {
        fprintf (stderr, "malloc error in malloc_matrix_f for %ld row pointers.\n", nr);
        m = NULL;
        return m;
    }
    /* Allocate rows and set the pointers to them */
    m[0] = (float *) malloc((size_t) (nr * nc * sizeof(float)));
    if (!m[0]) {
        fprintf (stderr, "malloc error in malloc_matrix_f for %ld row with %ld columns.\n", nr, nc);
        m[0] = NULL;
        free (m);
        m = NULL;
        return m;
    }
    for (i = 1; i < nr; i++) m[i] = m[i-1] + nc;

    return m;
}

void
free_matrix_f (float **m) 
{
    free (m[0]);
    free (m);
    return;
}



float ***
malloc_tensor_f (long nr, long nc, long nd)
{
    float ***t = NULL;
    long i, j;

    /* Allocate pointers to rows */
    t = (float ***) malloc((size_t) (nr * sizeof(float **)));
    if (!t) {
        fprintf (stderr, "malloc error in malloc_tensor_f for %ld row pointers.\n", nr);
        t = NULL;
        return t;
    }
    /* Allocate rows and set the pointers to them */
    t[0] = (float **) malloc((size_t) (nr * nc * sizeof(float *)));
    if (!t[0]) {
        fprintf (stderr, "malloc error in malloc_tensor_f for %ld row with %ld columns.\n", nr, nc);
        t[0] = NULL;

        free (t);
        t = NULL;
        return t;
    }
    /* Allocate rows and set the pointers to them */
    t[0][0] = (float *) malloc((size_t) (nr * nc *nd * sizeof(float)));
    if (!t[0][0]) {
        fprintf (stderr, "malloc error in malloc_tensor_f for %ld row with %ld columns and %ld deep.\n", nr, nc, nd);
        t[0][0] = NULL;
        free (t[0]);
        t[0] = NULL;
        free(t);
        t = NULL;
        return t;
    }

    for (j = 1; j < nc; j++) t[0][j] = t[0][j-1] + nd;
    for (i = 1; i < nr; i++) {
        t[i] = t[i-1] + nc;
        t[i][0] = t[i-1][0] + nc*nd;
        for (j = 1; j < nc; j++) t[i][j] = t[i][j-1] + nd;
    }
    return t;
}

void
free_tensor_f (float ***t) 
{
    free (t[0][0]);
    free (t[0]);
    free (t);
    return;
}

unsigned short ***
malloc_tensor_us (long nr, long nc, long nd)
{
    unsigned short ***t = NULL;
    long i, j;

    /* Allocate pointers to rows */
    t = (unsigned short ***) malloc((size_t) (nr * sizeof(unsigned short **)));
    if (!t) {
        fprintf (stderr, "malloc error in malloc_tensor_us for %ld row pointers.\n", nr);
        t = NULL;
        return t;
    }
    /* Allocate rows and set the pointers to them */
    t[0] = (unsigned short **) malloc((size_t) (nr * nc * sizeof(unsigned short *)));
    if (!t[0]) {
        fprintf (stderr, "malloc error in malloc_tensor_us for %ld row with %ld columns.\n", nr, nc);
        t[0] = NULL;

        free (t);
        t = NULL;
        return t;
    }
    /* Allocate rows and set the pointers to them */
    t[0][0] = (unsigned short *) malloc((size_t) (nr * nc *nd * sizeof(unsigned short)));
    if (!t[0][0]) {
        fprintf (stderr, "malloc error in malloc_tensor_us for %ld row with %ld columns and %ld deep.\n", nr, nc, nd);
        t[0][0] = NULL;
        free (t[0]);
        t[0] = NULL;
        free(t);
        t = NULL;
        return t;
    }

    for (j = 1; j < nc; j++) t[0][j] = t[0][j-1] + nd;
    for (i = 1; i < nr; i++) {
        t[i] = t[i-1] + nc;
        t[i][0] = t[i-1][0] + nc*nd;
        for (j = 1; j < nc; j++) t[i][j] = t[i][j-1] + nd;
    }
    return t;
}

void
free_tensor_us (unsigned short ***t) 
{
    free (t[0][0]);
    free (t[0]);
    free (t);
    return;
}

