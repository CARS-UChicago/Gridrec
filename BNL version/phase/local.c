 /* 
 * Read unfiltered backprojected reconstructed image, apply 2D Laplacian and 
   add together with mixing parameter * original image to obtain local correction.
 */

/* System includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>

/* Header file for netCDF */
#include <netcdf.h>
#include "local.h"

typedef enum TrueFalse { False, True } Flag;
Flag local, verbose, clip, median;

int
main (int argc, char *argv[]) {
    char CDFfile[256];
    float **Data2D = NULL;
    float **Buff2D = NULL;
    float **Medi2D = NULL;
    long files = 0;
    int xsize = 0, ysize = 0, zsize = 0;
    int x = 0, y = 0, z = 0;
    extern char *optarg;
    extern int optind, opterr;
    int optindr;
    int c, border = 1, radius = 0; 
    long xcut = 0, ycut = 0;
    float add = .01;
    int errflg = 0;
    float a[9] = {0,0,0,0,0,0,0,0,0};

    while ((c = getopt(argc, argv, "vma:b:")) != EOF)
        switch (c) {
        case 'v':
            verbose = True;
          break;
        case 'm':
            median = True;
          break;
	case 'a':
	    add = atof(optarg);
	    local = True;
	    break;
	case 'b':
	    border = atoi(optarg);
	    clip = True;
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

    Buff2D = readnetCDF(CDFfile, &xsize, &ysize);

    xcut = xsize - 2*border;
    if ((xcut % 2) == 0) (xcut)--;
    ycut = ysize - 2*border;
    if ((ycut % 2) == 0) (ycut)--;

    radius = xcut/2;
    if (verbose) printf ("radius = %d \n", radius);

/*    Data2D = malloc_matrix_f (xcut, ycut);

	for (x = border; x < xcut + border; x++) {
           for (y = border; y < ycut + border; y++) { 
		if (sqrt (pow((x - (xcut + border)/2),2) + pow((y - (ycut + border)/2),2)) <= radius){
               Data2D[x-border][y-border] = ((4 * Buff2D[x][y] - (Buff2D[x][y + 1] + Buff2D[x - 1][y] 
						+ Buff2D[x + 1][y] + Buff2D[x][y - 1])) + (add * Buff2D[x][y])); 
		} else {
		Data2D[x-border][y-border] = 0.0;
		}
            }
        }*/

    Data2D = malloc_matrix_f (xsize, ysize);

	for (x = 0; x < xsize; x++) {
           for (y = 0; y < ysize; y++) { 
		if (sqrt (pow((x - (xsize)/2),2) + pow((y - (ysize)/2),2)) <= radius){
                Data2D[x][y] = ((4 * Buff2D[x][y] - (Buff2D[x][y + 1] + Buff2D[x - 1][y] 
					+ Buff2D[x + 1][y] + Buff2D[x][y - 1])) + (add * Buff2D[x][y])); 
		} else {
		Data2D[x][y] = 0.0;
		}
            }
        }


     if (median) {
     Medi2D = malloc_matrix_f (xsize, ysize);

/*	for (x = border; x < xcut + border; x++) {
          for (y = border; y < ycut + border; y++) { 
		Medi2D[x - border][y - border] = Data2D[x][y];
            }
        }*/

	for (x = 0; x < xsize; x++) {
          for (y = 0; y < ysize; y++) { 
		Medi2D[x][y] = Data2D[x][y];
            }
        }

	for (x = border; x < xcut + border; x++) {
          for (y = border; y < ycut + border; y++) { 

                a[0] = Medi2D[x-1][y-1];
                a[1] = Medi2D[x-1][y];
                a[2] = Medi2D[x-1][y+1];
                a[3] = Medi2D[x][y-1];
                a[4] = Medi2D[x][y];
                a[5] = Medi2D[x][y+1];
                a[6] = Medi2D[x+1][y-1];
                a[7] = Medi2D[x+1][y];
                a[8] = Medi2D[x+1][y+1];
                
                shellsort (9, a);
                Data2D[x - border][y - border] = a[4];
           } 
	}
     free_matrix_f(Medi2D);
     }

     printf ("Writing local image to netCDF file...\n");
     writenetCDF (CDFfile, Data2D, xsize, ysize);

  }   /* end of for z loop*/

    free_matrix_f(Data2D);
    free_matrix_f(Buff2D);

  }   /*end of else */

        exit(0);
}


/*
  **************** end of the main ****************
*/


void
writenetCDF (char *name, float **m, int x, int y)
{
    char CDFname[256];
    char stoken[] = {"."};
    char *token;

    int ncid;
    int xid, yid;
    int slice_dimids[2], slice_id;
    long start[2], count[2];

   /* Create a netCDF file for writing */
    token = strtok(name, stoken);
    strcpy (CDFname, token);
    strcat (CDFname, ".ncdf");
    ncid = nccreate (CDFname,  NC_CLOBBER);

    /* Define the dimensions */
    yid = ncdimdef (ncid, "y", y);
    xid = ncdimdef (ncid, "x", x);
    slice_dimids[0] = yid;
    slice_dimids[1] = xid;

    /* Define the variables */
    slice_id = ncvardef (ncid, "slice", NC_FLOAT, 2, slice_dimids);
            /* Leave defined mode */
    ncendef (ncid);

    /* Write the data */
    start[0] = 0;
    start[1] = 0;
    count[0] = y;
    count[1] = x;
    ncvarput (ncid, slice_id, start, count, (void **)&m[start[0]][start[1]]);

        
    /* Close the netCDF file */

    ncclose (ncid);
    return;
}


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
    yid = ncdimid (ncid, "N1");
    ncdiminq (ncid, yid, (char *) 0, &ydim);
    *ysize = ydim;
    xid = ncdimid (ncid, "N2");
    ncdiminq (ncid, xid, (char *) 0, &xdim);
    *xsize = xdim;
    /*Get the variable id*/
    slice_dimids[0] = yid;
    slice_dimids[1] = xid;
    slice_id = ncvarid (ncid, "LBLREC");

    /*Allocate space for 2D array*/
    if (xdim && ydim)

    Data = malloc_matrix_f (ydim,xdim);

   /* Read the data */
    start[0] = 0;
    start[1] = 0;
    count[0] = ydim;
    count[1] = xdim;
    ncvarget (ncid, slice_id, start, count, (void *)&Data[start[0]][start[1]]);

    /* Close the netCDF file */
        ncclose (ncid);


    return Data;
}


void
shellsort (unsigned long n, float a[])
{
    unsigned long i, j, v, inc = 1;

    do {
        inc *=3;
        inc++;
    } while (inc <= n);

    do {
        inc /= 3;
        for (i = inc; i < n; i++) {
            v = a[i];
            j = i;
            while (a[j-inc] > v) {
                a[j] = a[j-inc];
                j -= inc;
                if (j < inc) break;
            }
            a[j] = v;
        }
    } while (inc > 1);
    return;
}

void
usage (void)
{
    fprintf(stderr, "usage: PMIStoCDF [-a air_values -i -v -w<file>] files...\n");
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



