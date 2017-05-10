 /* 
 * Read innetCDF files andchange formats to flat binary or ascii.
 */

/* System includes */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/stat.h>

/* Header file for netCDF */
#include <netcdf.h>

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

float **
opennetCDF (char *name, int *xsize, int *ysize);

void
writeascii (char *name, float **m, int xsize, int ysize);

void
writebin (char *name, float **m, int xsize, int ysize);

void
usage (void);

static int verbose = 0;
static int saveascii = 0;
static int savebinary = 0;

int
main (int argc, char *argv[]) {
    char CDFfile[256];
    float **Data2D = NULL;
    long files = 0;
    int xsize = 0, ysize = 0,zsize = 0;
    int x, y, z;
    extern char *optarg;
    extern int optind, opterr;
    int optindr;
    int c;
    float max,min;
    int errflg = 0;
    int ncid, xid, yid;
    long xdim, ydim;
    int slice_dimids[2], slice_id;
    long start[2], count[2];

    while ((c = getopt(argc, argv, "abv")) != EOF)
        switch (c) {
        case 'a':
            saveascii++;
            break;
        case 'b':
            savebinary++;
            break;
        case 'v':
            verbose++;
          break;
        case '?':
            errflg++;
        }

    if (errflg) {
        usage ();
    } else {
        printf("\t-v\tverbose = %d\n", verbose);
        printf("\t-a\tsaveascii = %d\n", saveascii);
        printf("\t-b\tsavebinary = %d\n", savebinary);
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

    Data2D = opennetCDF (CDFfile, &xsize, &ysize);

    if (saveascii) writeascii (CDFfile, Data2D, xsize, ysize);

    if (savebinary) writebin (CDFfile, Data2D, xsize, ysize);

  }   /* end of for z loop*/

        free_matrix_f(Data2D);

  }   /*end of else */

        exit(0);
}


/*
  **************** This is the end of the main program ****************
*/


void
writeascii (char *name, float **m, int xsize, int ysize)

{
    char ASCname[256];
    FILE *fo;
    int i,j,cnt;

    strcpy (ASCname, name);
    strcat (ASCname, ".asc");

    if (verbose) printf ("Saving data to file %s.\n", ASCname);

    if ((fo=fopen(ASCname,"wb")) == NULL)
        printf("can't open %s for writing .\n",ASCname);

    else{

    /* Write the data */

     	   for (j=0 ; j < ysize; j++){
        	cnt=1;
		for (i=0 ; i < xsize; i++){
	     		fprintf (fo, "% .10f\t", m[j][i]);
	     		if (cnt%5 == 0 || cnt%xsize == 0)
	     		fprintf(fo, "\n");
	     		cnt++;
		}
    	    }
        } 
    fclose(fo);
    return;
}

void
writebin (char *name, float **m, int xsize, int ysize)

{
    char BINname[256];
    FILE *fo;
    int i,j;

    strcpy (BINname, name);
    strcat (BINname, ".bin");

    if (verbose) printf ("Saving data to file %s.\n", BINname);

    if ((fo=fopen(BINname,"wb")) == NULL)
        printf("can't open %s for writing .\n",BINname);
    else{

    /* Write the data */

	fwrite (m[0], sizeof(float), xsize*ysize, fo);
    }
    fclose(fo);
    return;
}

float **
opennetCDF (char *name, int *xsize, int *ysize)

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
usage (void)
{
    fprintf(stderr, "usage: VOLCDF [-v] files...\n");
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


