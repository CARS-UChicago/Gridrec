/* 
 * Read header and data from PMIS Image File and convert
 * to netCDF format.
 *
 * The image file consists of a file header immediately
 * followed by image data, which is in 16-bit format.
 * The number of image frames, columns, and rows must be 
 * read from the file header.
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

#include "PMIStoCDF.h"

static int center = 0;
static int sinogram = 0;
static int verbose = 0;
static int saveimage = 0;
static float imagescale = 10000.;
static float deltatheta = 1.0;
static int n_airvalues = 1;

int
main (int argc, char *argv[]) {
    char PMISfile[256], whitefield[MAXWHITEFIELD][256];
    struct PMISHeader Header;
    unsigned short **Data2D = NULL, **WhiteField = NULL;
    float ***Data3D = NULL;
    int files = 0;
    int xsize = 0, ysize = 0, txsize = 0, tysize = 0, nangles = 0, angle = 0;
    int x, y;
    float min, max;
    int n_whitefield = 0;
    extern char *optarg;
    extern int optind, opterr;
    int optindr;
    int c;
    int errflg = 0;

    while ((c = getopt(argc, argv, "a:d:cisS:w:v")) != EOF)
        switch (c) {
        case 'a':
            n_airvalues = atoi(optarg);
            break;
        case 'c':
            center++;
            break;
        case 'd':
            deltatheta = atof(optarg);
            break;
        case 'i':
            saveimage++;
            break;
        case 's':
            sinogram++;
            break;
        case 'S':
            imagescale = atof(optarg);
            break;
        case 'v':
            verbose++;
            break;
        case 'w':
            if (n_whitefield < MAXWHITEFIELD) {
                strcpy(whitefield[n_whitefield], optarg);
                n_whitefield++;
            } else {
                fprintf(stderr, "Too many white field images specified\n."
                        "Maximum of %d white field images.\n", MAXWHITEFIELD);
            }
            break;
        case '?':
            errflg++;
        }

    if (errflg) {
        usage ();
    } else {
        printf("PMIStoCDF parameters\n====================\n");
        printf("\t-c\tcenter = %d\n", center);
        printf("\t-d\tdeltatheta = %f\n", deltatheta);
        printf("\t-S\timagescale = %f\n", imagescale);
        printf("\t-i\tsaveimage = %d\n", saveimage);
        printf("\t-s\tsinogram = %d\n", sinogram);
        printf("\t-v\tverbose = %d\n", verbose);
        printf("\n");
    }

    nangles = files = argc - optind;
    if (files == 0) {
        fprintf (stderr, "No filename was specified.\n");
        exit(1);
    } else {
        strcpy (PMISfile, argv[optind]);
        if (fexist(PMISfile)) {
            if (verbose) printf ("Reading PMI file %s for header information.\n", PMISfile);
            /* Read the PMIS file */
            Data2D = readPMIS (PMISfile, &Header, &xsize, &ysize);
            txsize = xsize;
            tysize = ysize;
            free_matrix_us (Data2D);
        }
        if (xsize && ysize && nangles) {
            Data3D = malloc_tensor_f (ysize, nangles, xsize);
            if (!Data3D) {
                fprintf (stderr, "Not enough space to allocate 3D array.\n");
                exit(3);
            }
        }
    }
    

    if (n_whitefield) {
        int nw;
        unsigned short **Tmp = NULL;
        
        if (verbose) printf("Attempting to read %d white field images.\n", n_whitefield);
        for (nw = 0; nw < n_whitefield; nw++) {
            if (fexist(whitefield[nw])) {
                if (verbose) printf ("Reading PMI white field file %s\n", whitefield[nw]);

                /* Read the PMIS whitefield file */
                Tmp = readPMIS (whitefield[nw], &Header, &xsize, &ysize);

                if (xsize != txsize || ysize != tysize) {
                    fprintf(stderr, "Image sizes are not consistent:\n"
                            "\tprevious: xsize ysize = %d %d\n"
                            "\tcurrent : xsize ysize = %d %d\n",
                            txsize, tysize, xsize, ysize);
                    exit(2);
                }
                    
                if (!WhiteField) {
                    WhiteField = malloc_matrix_us(ysize, xsize);
                    for (y = 0; y < ysize; y++) {
                        for (x = 0; x < xsize; x++) {
                            WhiteField[y][x] = Tmp[y][x];
                        }
                    }
                } else {
                    for (y = 0; y < ysize; y++) {
                        for (x = 0; x < xsize; x++) {
                            WhiteField[y][x] += Tmp[y][x];
                        }
                    }
                }
                free_matrix_us(Tmp);
            } else {
                fprintf (stderr, "File %s does not exist. Aborting execution.\n", whitefield[nw]);
                exit(2);
            }
        }

        /* Get the average white field value for all images */
        for (y = 0; y < ysize; y++) {
            for (x = 0; x < xsize; x++) {
                WhiteField[y][x] /= n_whitefield;
            }
        }
        if (saveimage) saveaspgm ("whitefield", WhiteField, xsize, ysize);
    }

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
    for (angle = 0; optind < argc; optind++, angle++) {
        strcpy (PMISfile, argv[optind]);

        /* Read the PMIS file */
        if (verbose) printf ("Reading PMI file %s\n", PMISfile);
        Data2D = readPMIS (PMISfile, &Header, &xsize, &ysize);
        if (xsize != txsize || ysize != tysize) {
            fprintf(stderr, "Image sizes are not consistent:\n"
                    "\tprevious: xsize ysize = %d %d\n"
                    "\tcurrent : xsize ysize = %d %d\n",
                    txsize, tysize, xsize, ysize);
            exit(2);
        }
        
        /* Write the image into the 3D file */
        for (y = 0; y < ysize; y++) {
            for (x = 0; x < xsize; x++) {
                if (WhiteField) {
                    Data3D[y][angle][x] = imagescale*Data2D[y][x]/WhiteField[y][x];
                } else {
                    Data3D[y][angle][x] = imagescale*Data2D[y][x];
                }
            }
        }

        if (saveimage) {
            for (y = 0; y < ysize; y++) {
                for (x = 0; x < xsize; x++) {
                    if (WhiteField) {
                        Data2D[y][x] = imagescale*Data2D[y][x]/WhiteField[y][x];
                    } else {
                        Data2D[y][x] = imagescale*Data2D[y][x];
                    }
                }
            }
            saveaspgm (PMISfile, Data2D, xsize, ysize);
        }
        free_matrix_us(Data2D);
    }

    if (sinogram) {
        int xl, xr;
        float air  = 0;

        /* Calculate air values */
        for (y = 0; y < ysize; y++) {
            for (angle = 0; angle < nangles; angle++) {
                for (xl = 0, xr = xsize - 1; xl < n_airvalues; xl++, xr--) {
                    air += (Data3D[y][angle][xl] + Data3D[y][angle][xr]);
                }
                air = air / (2.0 * n_airvalues);
                for (x = 0; x < xsize; x++) {
                    Data3D[y][angle][x] = - imagescale * log (Data3D[y][angle][x] / air);
                }
            }
        }
    }

    if (center) {
        float *cog = NULL;
	float *shiftvals = NULL;

        /* Calculate air values */
        cog = malloc_vector_f (nangles);
	shiftvals = malloc_vector_f (xsize);
        for (y = 0; y < ysize; y++) {
 	    int shift = 0;
            for (angle = 0; angle < nangles; angle++) {
                float sum = 0, sumx = 0;
                for (x = 0; x < xsize; x++) {
                    sum += Data3D[y][angle][x];
                    sumx += (x+1)*Data3D[y][angle][x];
                }
                cog[angle] = sumx / sum;
            }
	    /*for (angle = 0; angle < nangles; angle++) {
                printf ("angle b_cog = %f %f\n", angle*deltatheta, cog[angle]);
            }*/
	    /* Make an arbitrary shift here to see what moving the
	       center of gravity does */
            for (angle = 0; angle < nangles; angle++) {
                for (x = 0; x < xsize; x++) {
	            shiftvals[x] = Data3D[y][angle][x];
                }
                for (x = 0; x < xsize; x++) {
	            int xshift = x + shift;
	            if (xshift >= 0 && xshift < xsize) {
	                Data3D[y][angle][xshift] = shiftvals[x];
		    }
                }
            }
            for (angle = 0; angle < nangles; angle++) {
                float sum = 0, sumx = 0;
                for (x = 0; x < xsize; x++) {
                    sum += Data3D[y][angle][x];
                    sumx += (x+1)*Data3D[y][angle][x];
                }
                cog[angle] = sumx / sum;
            }
	    /*for (angle = 0; angle < nangles; angle++) {
                printf ("angle a_cog = %f %f\n", angle*deltatheta, cog[angle]);
            }*/
        }
	free_vector_f (shiftvals);
        free_vector_f (cog);
    }

    /* Now write out the file in horizontal slices and see what we get */
    max = 0.;
    min = 1000000.;
    for (y = 0; y < ysize; y++) {
        for (angle = 0; angle < nangles; angle++) {
            for (x = 0; x < xsize; x++) {
                if (Data3D[y][angle][x] > max) max = Data3D[y][angle][x];
                if (Data3D[y][angle][x] < min) min = Data3D[y][angle][x];
            }
        }
    }
    printf("min max = %f %f\n", min, max);

    /* Write the slices to a netCDF file */
    printf ("Writing slices to netCDF files...\n");
    for (y = 0; y < ysize; y++) {
        char name[256];
        sprintf(name, "slice%3.3d.nc", y);
        saveasnetCDF (name, Data3D[y], xsize, nangles);
        if (saveimage) {
            sprintf(name, "slice%3.3d", y);
            saveaspgm2f (name, Data3D[y], xsize, nangles, min, max);
        }
    }

    exit(0);
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

long
byteswap (long x, int nbytes) {
    long y = 0;
    unsigned long maskl, maskr;
    int i, n, nd2;

    if (nbytes != (nbytes/2)*2) {
        fprintf (stderr, "Fatal error:"
                " function byteswap only works for an even number of bytes.\n");
        exit(2);
    }

    n = nbytes;
    nd2 = nbytes / 2;
    maskr = ~(~0 << 8);
    maskl = maskr << 8*(nbytes-1);
    
    for (i = 0; i < nd2; i++, n-=2) {
        int yr = 0, yl = 0;
        yr = (x & maskr) << 8*(n-1);
        yl = (x & maskl) >> 8*(n-1);
        y |= yr;
        y |= yl;
        maskr <<= 8;
        maskl >>= 8;
    }
    return y;
}

unsigned short **
readPMIS (char *PMISfile, struct PMISHeader *Header, int *xsize, int *ysize)
{
    FILE *fp;
    size_t size;
    struct tm caltime;
    int arraysize, datasize;
    unsigned short **buffer = NULL, **Data = NULL;
    unsigned short d;
    unsigned int x, y, a[9];
    
    fp = fopen (PMISfile, "r");
    
    size = fread (&Header->ID, sizeof Header->ID, 1, fp);
    size = fread (&Header->nHeadSize, sizeof Header->nHeadSize, 1, fp);
    Header->nHeadSize = BYTESWAP2(Header->nHeadSize);
    size = fread (&Header->nVersion, sizeof Header->nVersion, 1, fp);
    Header->nVersion = BYTESWAP2(Header->nVersion);
    size = fread (&Header->nXorg, sizeof Header->nXorg, 1, fp);
    Header->nXorg = BYTESWAP2(Header->nXorg);
    size = fread (&Header->nYorg, sizeof Header->nYorg, 1, fp);
    Header->nYorg = BYTESWAP2(Header->nYorg);
    size = fread (&Header->nXsiz, sizeof Header->nXsiz, 1, fp);
    Header->nXsiz = BYTESWAP2(Header->nXsiz);
    size = fread (&Header->nYsiz, sizeof Header->nYsiz, 1, fp);
    Header->nYsiz = BYTESWAP2(Header->nYsiz);
    size = fread (&Header->nXbin, sizeof Header->nXbin, 1, fp);
    Header->nXbin = BYTESWAP2(Header->nXbin);
    size = fread (&Header->nYbin, sizeof Header->nYbin, 1, fp);
    Header->nYbin = BYTESWAP2(Header->nYbin);
    size = fread (&Header->szName, sizeof Header->szName, 1, fp);
    size = fread (&Header->szComment, sizeof Header->szComment, 1, fp);
    size = fread (&Header->tCreated, sizeof Header->tCreated, 1, fp);
    Header->tCreated = byteswap(Header->tCreated, 4);
    caltime = *localtime (&Header->tCreated);
    size = fread (&Header->tModified, sizeof Header->tModified, 1, fp);
    Header->tModified = byteswap(Header->tModified, 4);
    caltime = *localtime (&Header->tModified);
    size = fread (&Header->nGain, sizeof Header->nGain, 1, fp);
    Header->nGain = BYTESWAP2(Header->nGain);
    size = fread (&Header->nImages, sizeof Header->nImages, 1, fp);
    Header->nImages = BYTESWAP2(Header->nImages);

    if (verbose) {
        printf ("ID = %s\n", Header->ID);
        printf ("nHeadSize = %d\n", Header->nHeadSize);
        printf ("nVersion = %d\n", Header->nVersion);
        printf ("nXorg = %d\n", Header->nXorg);
        printf ("nYorg = %d\n", Header->nYorg);
        printf ("nXsiz = %d\n", Header->nXsiz);
        printf ("nYsiz = %d\n", Header->nYsiz);
        printf ("nXbin = %d\n", Header->nXbin);
        printf ("nYbin = %d\n", Header->nYbin);
        printf ("szName = %s\n", Header->szName);
        printf ("szComment = %s\n", Header->szComment);
        printf ("tCreated = %ld\n", Header->tCreated);
        printf ("year month day = %d %d %d\n", caltime.tm_year, caltime.tm_mon, caltime.tm_mday);
        printf ("tModified = %ldd\n", Header->tModified);
        printf ("year month day = %d %d %d\n", caltime.tm_year, caltime.tm_mon, caltime.tm_mday );
        printf ("nGain = %d\n", Header->nGain);
        printf ("nImages = %d\n", Header->nImages);
    }
    
    /* The data size is the x size times the y size times 2 for the number of bytes. */
    arraysize = Header->nXsiz * Header->nYsiz;
    datasize = arraysize * sizeof (unsigned short);
    /* The data buffer must be freed outside of this function. */
    /* buffer = malloc (datasize); */
    buffer = malloc_matrix_us(Header->nYsiz, Header->nXsiz);

    size = fread (buffer[0], sizeof d, arraysize, fp);
    if (verbose) printf ("size = %d\n", size);
    /* for (i = 0; i < arraysize; i++) buffer[i] = BYTESWAP2(buffer[i]); */
    for (y = 0; y < Header->nYsiz; y++) {
        for (x = 0; x < Header->nXsiz; x++) {
            buffer[y][x] = BYTESWAP2(buffer[y][x]);
        }
    }
    
    /* Apply the median filter - use 3x3 sampling */
    *xsize = Header->nXsiz - 2*XBORDER;
    if ((*xsize % 2) == 0) (*xsize)--;
    *ysize = Header->nYsiz - 2*YBORDER;
    if ((*ysize % 2) == 0) (*ysize)--;
    Data = malloc_matrix_us(*ysize, *xsize);

    if (verbose) printf("Original: xsize ysize = %d %d, Final: xsize ysize = %d %d\n",
           Header->nXsiz, Header->nYsiz, *xsize, *ysize);
    for (y = YBORDER; y < *ysize + YBORDER; y++) {
        for (x = XBORDER; x < *xsize + XBORDER; x++) {
            if (MEDIAN) {
                a[0] = buffer[y-1][x-1];
                a[1] = buffer[y-1][x];
                a[2] = buffer[y-1][x];
                a[3] = buffer[y][x-1];
                a[4] = buffer[y][x];
                a[5] = buffer[y][x+1];
                a[6] = buffer[y+1][x-1];
                a[7] = buffer[y+1][x];
                a[8] = buffer[y+1][x+1];
                
                shellsort (9, a);
                /* qsort is really slow ! */
                /* qsort (&a[0], 9, sizeof (int), (int (*)(const void *, const void *))numcmp); */
                /* printf("x y a[0-8] = %3d %3d %5d %5d %5d %5d %5d %5d %5d %5d %5d\n",
                   x, y, a[0], a[1], a[2], a[3], a[4], a[5], a[6], a[7], a[8]); */
                Data[y-YBORDER][x-XBORDER] = a[4];
            } else {
                Data[y-YBORDER][x-XBORDER] = buffer[y][x];
            }
        }
    }
    free_matrix_us (buffer);
    fclose (fp);
    return Data;
}


/* Sorting functions */
int
numcmp (char *s1, char *s2)
{
    unsigned short us1, us2;

    us1 = atoi(s1);
    us2 = atoi(s2);
    if (us1 < us2) {
        return -1;
    } else if (us1 > us2) {
        return 1;
    } else {
        return 0;
    }
}

void
shellsort (unsigned long n, unsigned int a[])
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

/* Memory allocation routines */

unsigned short *
malloc_vector_us (long n) 
{
    unsigned short *v = NULL;
    
    v = (unsigned short *) malloc((size_t) (n * sizeof(unsigned short)));
    if (!v) {
        fprintf (stderr, "malloc error in malloc_vector_us for length %ld.\n", n);
        v = NULL;
        return v;
    }
    return v;
}

void
free_vector_us (unsigned short *v) 
{
    free (v);
    return;
}

float *
malloc_vector_f (long n) 
{
    float *v = NULL;
    
    v = (float *) malloc((size_t) (n * sizeof(float)));
    if (!v) {
        fprintf (stderr, "malloc error in malloc_vector_f for length %ld.\n", n);
        v = NULL;
        return v;
    }
    return v;
}

void
free_vector_f (float *v) 
{
    free (v);
    return;
}

unsigned short **
malloc_matrix_us (long nr, long nc)
{
    unsigned short **m = NULL;
    long i;

    /* Allocate pointers to rows */
    m = (unsigned short **) malloc((size_t) (nr * sizeof(unsigned short *)));
    if (!m) {
        fprintf (stderr, "malloc error in malloc_matrix_us for %ld row pointers.\n", nr);
        m = NULL;
        return m;
    }
    /* Allocate rows and set the pointers to them */
    m[0] = (unsigned short *) malloc((size_t) (nr * nc * sizeof(unsigned short)));
    if (!m[0]) {
        fprintf (stderr, "malloc error in malloc_matrix_us for %ld row with %ld columns.\n", nr, nc);
        m[0] = NULL;
        free (m);
        m = NULL;
        return m;
    }
    for (i = 1; i < nr; i++) m[i] = m[i-1] + nc;

    return m;
}

void
free_matrix_us (unsigned short **m) 
{
    free (m[0]);
    free (m);
    return;
}

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

void
saveaspgm (char *name, unsigned short **m, int x, int y)
{
#define MAXGRAYS 255
    char PGMname[256];
    FILE *fp;
    int i, j, max = 0, min = 1000000;

    strcpy (PGMname, name);
    strcat (PGMname, ".pgm");
    if (verbose) printf ("Saving data to image file %s.\n", PGMname);
    fp = fopen (PGMname, "w");
    fprintf (fp, "P2\n#CREATOR: %s\n", VERSION);
    fprintf (fp, "%d %d\n%d\n", x, y, MAXGRAYS);

    for (j = 0; j < y; j++) {
        for (i = 0; i < x; i++) {
            if (m[j][i] > max) max = m[j][i];
            if (m[j][i] < min) min = m[j][i];
        }
    }
    
    for (j = 0; j < y; j++) {
        for (i = 0; i < x; i++) {
            fprintf(fp, "%3d ", MAXGRAYS * (m[j][i] - min) / (max - min));
        }
    }
    
    fclose(fp);

    return;
}

void
saveaspgm2 (char *name, unsigned short **m, int x, int y, int min, int max)
{
#define MAXGRAYS 255
    char PGMname[256];
    FILE *fp;
    int i, j;

    strcpy (PGMname, name);
    strcat (PGMname, ".pgm");
    if (verbose) printf ("Saving data to image file %s.\n", PGMname);
    fp = fopen (PGMname, "w");
    fprintf (fp, "P2\n#CREATOR: %s\n", VERSION);
    fprintf (fp, "%d %d\n%d\n", x, y, MAXGRAYS);

    for (j = 0; j < y; j++) {
        for (i = 0; i < x; i++) {
            fprintf(fp, "%3d ", MAXGRAYS * (m[j][i] - min) / (max - min));
        }
    }
    
    fclose(fp);

    return;
}

void
saveaspgm2f (char *name, float **m, int x, int y, float min, float max)
{
#define MAXGRAYS 255
    char PGMname[256];
    FILE *fp;
    int i, j;

    strcpy (PGMname, name);
    strcat (PGMname, ".pgm");
    if (verbose) printf ("Saving data to image file %s.\n", PGMname);
    fp = fopen (PGMname, "w");
    fprintf (fp, "P2\n#CREATOR: %s\n", VERSION);
    fprintf (fp, "%d %d\n%d\n", x, y, MAXGRAYS);

    for (j = 0; j < y; j++) {
        for (i = 0; i < x; i++) {
            fprintf(fp, "%3d ", (int) (MAXGRAYS * (m[j][i] - min) / (max - min)));
        }
    }
    
    fclose(fp);

    return;
}

void
saveasnetCDF (char *name, float **m, int x, int y)
{
    char CDFname[256];
    int i, j, ncid;
    int xid, angleid;
    int sinogram_dimids[2], sinogram_id;
    long start[2], count[2];

    strcpy (CDFname, name);
    if (verbose) printf ("Saving slice data to file %s.\n", CDFname);
    
    /* Create a netCDF file for writing */
    ncid = nccreate (CDFname,  NC_CLOBBER);

    /* Define the dimensions */
    /*angleid = ncdimdef (ncid, "angle", NC_UNLIMITED);*/
    angleid = ncdimdef (ncid, "angle", y);
    xid = ncdimdef (ncid, "x", x);

    sinogram_dimids[0] = angleid;
    sinogram_dimids[1] = xid;

    /* Define the variables */
    sinogram_id = ncvardef (ncid, "sinogram", NC_FLOAT, 2, sinogram_dimids);
        
    /* Leave defined mode */
    ncendef (ncid);

    /* Write the data */
    start[0] = 0;
    start[1] = 0;
    count[0] = y;
    count[1] = x;
    ncvarput (ncid, sinogram_id, start, count, (void *)&m[start[0]][start[1]]);
        
    /* Close the netCDF file */
    ncclose (ncid);


    return;
}

