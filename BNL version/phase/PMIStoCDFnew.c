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
static int removerings = 0;
static int saveimage = 0;
static float imagescale = 10000.;
static float deltatheta = 1;
static int n_airvalues = 1;
static int median = 1;


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

    while ((c = getopt(argc, argv, "a:d:m:cirsS:w:v")) != EOF)
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
        case 'r':
            removerings++;
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
        case 'm':
            median = atof(optarg);
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
        printf("\t-r\tremoverings = %d\n", removerings);
        printf("\t-m\tmedian = %d\n", median);
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

        /* Calculate air values */
        for (y = 0; y < ysize; y++) {
            for (angle = 0; angle < nangles; angle++) {
                float air  = 0;
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
	float *theta = NULL;
	float avsin = 0;
        int i;
      
        /* Calculate air values */
        cog = vector(1,nangles);
	shiftvals = malloc_vector_f (xsize);
	theta = vector(1,nangles);
	/*printf ("x = %i\n", xsize);*/
	 /*printf ("y = %i\n", ysize);*/

       for (y = 0; y < ysize; y++) {
 	    float shift = 0;
 	    int xshift = 0;	    
            float halfsize = xsize/2;
	    float a[4]=
	    {0.0,150.0,15.0,170.0};
	     a[1] = halfsize;
	    
	    printf ("y = %i\n", y);
          
	      for (angle = 0;angle < nangles; angle++) {
                float sum = 0, sumx = 0;
                   for (x = 0; x < xsize; x++) {
                      sum += Data3D[y][angle][x];
                      sumx += (x+1)*Data3D[y][angle][x];
                   }
                    cog[angle] = sumx / sum;
		    theta[angle] = angle*deltatheta;
                  /*printf ("%f\t%f\n", angle*deltatheta, cog[angle]);*/
		    if(cog[angle] < 0)
		      cog[angle] *= -1;
		    
		    avsin +=cog[angle];
		    
      		   /* if(cog[angle] < 0){
	              printf("\n  negative cog values! slice y = %i was skipped\n",y);
	              goto cogerr;
	            }*/
	      }
		  (float) avsin /= nangles;
                  printf ("%f\n", avsin);
	    

	    if (cgfit(theta,cog,nangles,a)==1){
	       printf("Fit did not converge for slice y = %i",y);
	       continue;
	    }

	    else {
	    printf("%8s %8s %8s\n","a[1]","a[2]","a[3]");
	      for (i=1;i<=3;i++) printf("%9.4f",a[i]);
	    printf("\n");
	     
	      /* Shift the center of gravity */
                for (angle = 0; angle < nangles; angle++) {
                  for (x = 0; x < xsize; x++) {
	            shiftvals[x] = Data3D[y][angle][x];
                  }
                    for (x = 0; x < xsize; x++) {
	              shift = (halfsize - a[1]);
	              /*shift = (halfsize - avsin);*/
		         xshift = (x + shift + 0.5);
	                if (xshift >= 0 && xshift < xsize) {
	                  Data3D[y][angle][xshift] = shiftvals[x];
		        }
                    }
                }
	    }

	   /* Recalculate the center of gravity for debugging*/
             /* for (angle = 0; angle < nangles; angle++) {
                float sum = 0, sumx = 0;
                for (x = 0; x < xsize; x++) {
                    sum += Data3D[y][angle][x];
                    sumx += (x+1)*Data3D[y][angle][x];
                }
                cog[angle] = sumx / sum;
                printf ("%f\t%f\n", angle*deltatheta, cog[angle]);
              }
	    cgfit(theta,cog,nangles,a);
	    printf("%8s %8s %8s\n","a[1]","a[2]","a[3]");
	    for (i=1;i<=3;i++) printf("%9.4f",a[i]);
	    printf("\n");*/
	 /* cogerr: continue;*/
        }
	free_vector_f (shiftvals);
        free_vector(theta,1,nangles);
        free_vector(cog,1,nangles);
   }

    if (removerings) {
        float *avg = NULL;
	float *diff = NULL;
	float *smooth = NULL;
	/* The smoothing width must be odd */
	int halfsmoothwidth = 5;
	/*int smoothwidth = 2*halfsmoothwidth + 1;*/

	avg = malloc_vector_f (xsize);
	diff = malloc_vector_f (xsize);
	smooth = malloc_vector_f (xsize);
        for (y = 0; y < ysize; y++) {
 	   /* int shift = 0;*/
	    for (x = 0; x < xsize; x++) {
	        avg[x] = 0.;
                diff[x] = 0.;
	        smooth[x] = 0.;
            }

	    /* For each column, sum up all of the angles, and average */
            for (angle = 0; angle < nangles; angle++) {
                for (x = 0; x < xsize; x++) {
	            avg[x] += Data3D[y][angle][x];
                }
            }

            for (x = 0; x < xsize; x++) {
	        avg[x] /= nangles;
            }

            /* Smooth the average */
	    for (x = 0; x < xsize; x++) {
	       int ismooth, width = 0;
               float sum = 0;
	       for (ismooth = -halfsmoothwidth; ismooth <= halfsmoothwidth; ismooth++) {
	           int index = ismooth + x;
	           if (index >= 0 && index < xsize) {
                      width++;
                      sum += avg[index];
                   }
               }
	       smooth[x] = sum/width;
            }

            /* Make the difference array */
            for (x = 0; x < xsize; x++) {
	        diff[x] = avg[x] - smooth[x];
            }

            /* Subtract the difference from the 3D data set */
            for (angle = 0; angle < nangles; angle++) {
                for (x = 0; x < xsize; x++) {
	            Data3D[y][angle][x] -= diff[x];
                }
            }
        }
	free_vector_f (smooth);
        free_vector_f (diff);
        free_vector_f (avg);
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
    if (verbose) printf ("size = %d\n", (int)size);
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
            if (median) {
                a[0] = buffer[y-1][x-1];
                a[1] = buffer[y-1][x];
                a[2] = buffer[y-1][x+1];
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


#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

int nrfail(char error_text[])

{
        fprintf(stderr,"Bad values encountered during matrix inversion\n");
	fprintf(stderr,"%s\n",error_text);
        return 1;
}

void nrerror(char error_text[])

/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
        exit(1);
}

float *vector(long nl, long nh)

/* allocate a float vector with subscript range v[nl..nh] */

{
	float *v;
	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)

/* allocate an int vector with subscript range v[nl..nh] */

{
	int *v;
	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)

/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */

{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */

	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */

	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */

	return m;

}

void free_vector(float *v, long nl, long nh)

/* free a float vector allocated with vector() */

{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)

/* free an int vector allocated with ivector() */

{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)

/* free a float matrix allocated by matrix() */

{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}


#undef SWAP
#undef NRANSI


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
    int ncid;
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

/*NR Least Squares MRQ fitting routine*/

#define MA 3
#define PI 3.1415926536
#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}


int cgfit(float x[],float y[],int ndata,float a[])
{
	int i,*ia,iter,itst,k,mfit=MA,angle;
	float alamda,chisq,ochisq,*sig,**covar,**alpha;
	/*static float a[MA+1]=
	{0.0,150.0,15.0,170.0};*/
	
		ia=ivector(1,MA);
		/*x=vector(1,ndata);
		y=vector(1,ndata);*/
		sig=vector(1,ndata);
		covar=matrix(1,MA,1,MA);
		alpha=matrix(1,MA,1,MA);

	    for (angle = 0; angle < ndata; angle++) {
		sig[angle]=sqrt(y[angle]);
	    }

	for (i=1;i<=mfit;i++) ia[i]=1;
	for (iter=1;iter<2;iter++){
	    alamda = -1;
	   if(mrqmin(x,y,sig,ndata,a,ia,MA,covar,alpha,&chisq,fsin,&alamda)==1) return 1;
	  /* else*/
	    k=1;
	    itst=0;

	    for(;;){
	       /* printf("\n%s %2d %17s %10.4f %10s %9.2e\n","Iteration #",k,
      	       	"chi-squared:",chisq,"alamda:",alamda);
	        printf("%8s %8s %8s\n",	"a[1]","a[2]","a[3]");
	      for (i=1;i<=3;i++) printf("%9.4f",a[i]);
	      printf("\n");*/
	      k++;
	      ochisq=chisq;
	      mrqmin(x,y,sig,ndata,a,ia,MA,covar,alpha,&chisq,fsin,&alamda);
	      if (chisq > ochisq)
		itst=0;
	      else if (fabs(ochisq-chisq)< 0.1)
		itst++;
	      if (itst < 4)
		continue;
	      else if (k == 50){
		printf("\nConvergence not obtained after 50 iterations");
		return 1;
       	      }

	      alamda=0.0;
	      mrqmin(x,y,sig,ndata,a,ia,MA,covar,alpha,&chisq,fsin,&alamda);
	     /* printf("\nUncertainties:\n");
	      for(i=1;i<=3;i++) printf("%9.4f",sqrt(covar[i][i]));
	      printf("\n");*/
	      break;
	    }
	}
	free_matrix(alpha,1,MA,1,MA);
	free_matrix(covar,1,MA,1,MA);
	free_vector(sig,1,ndata);
	/*free_vector(y,1,ndata);
	free_vector(x,1,ndata);*/
	free_ivector(ia,1,MA);
	return 0;
}
			

int mrqmin(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int), float *alamda)
{

	int j,k,l,m;
	/*int i;*/
	static int mfit;
	static float ochisq,*atry,*beta,*da,**oneda;

	if (*alamda < 0.0) {
		atry=vector(1,ma);
		beta=vector(1,ma);
		da=vector(1,ma);
		for (mfit=0,j=1;j<=ma;j++)
			if (ia[j]) mfit++;
		oneda=matrix(1,mfit,1,1);
		*alamda=0.001;
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,fsin);
	     /* printf("\nmatrix alpha\n");
	      for(i=1;i<=MA;i++){
		for(j=1;j<=MA;j++)
		  printf("%12.4f",alpha[i][j]);
		printf("\n");
	      }*/
	     /* printf("\nvector beta\n");
	      for(i=1;i<=MA;i++) printf("%12.4f",beta[i]);*/
	      /*printf("\nchi-squared: %12.4f\n",chisq);*/
		ochisq=(*chisq);
		for (j=1;j<=ma;j++) atry[j]=a[j];
	}
	for (j=0,l=1;l<=ma;l++) {
		if (ia[l]) {
			for (j++,k=0,m=1;m<=ma;m++) {
				if (ia[m]) {
					k++;
					covar[j][k]=alpha[j][k];
				}
			}
			covar[j][j]=alpha[j][j]*(1.0+(*alamda));
			oneda[j][1]=beta[j];
		}
	}
	if (gaussj(covar,mfit,oneda,1)==1)return 1;
	for (j=1;j<=mfit;j++) da[j]=oneda[j][1];
	if (*alamda == 0.0) {
		covsrt(covar,ma,ia,mfit);
		free_matrix(oneda,1,mfit,1,1);
		free_vector(da,1,ma);
		free_vector(beta,1,ma);
		free_vector(atry,1,ma);
		return;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,fsin);
	if (*chisq < ochisq) {
		*alamda *= 0.1;
		ochisq=(*chisq);
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				for (j++,k=0,m=1;m<=ma;m++) {
					if (ia[m]) {
						k++;
						alpha[j][k]=covar[j][k];
					}
				}
				beta[j]=da[j];
				a[l]=atry[l];
			}
		}
	} else {
		*alamda *= 10.0;
		*chisq=ochisq;
	}
return 0;
}


void mrqcof(float x[], float y[], float sig[], int ndata, float a[], int ia[],
	int ma, float **alpha, float beta[], float *chisq,
	void (*funcs)(float, float [], float *, float [], int))

{

	int i,j,k,l,m,mfit=0;
	float ymod,wt,sig2i,dy,*dyda;


	dyda=vector(1,ma);
	for (j=1;j<=ma;j++)
		if (ia[j]) mfit++;
	for (j=1;j<=mfit;j++) {
		for (k=1;k<=j;k++) alpha[j][k]=0.0;
		beta[j]=0.0;
	}
	
	*chisq=0.0;
	
	for (i=0;i< ndata;i++) {
		(*fsin)(x[i],a,&ymod,dyda,ma);
		sig2i=1.0/(sig[i]*sig[i]);
		dy=y[i]-ymod;
		for (j=0,l=1;l<=ma;l++) {
			if (ia[l]) {
				wt=dyda[l]*sig2i;
				for (j++,k=0,m=1;m<=l;m++)
					if (ia[m]) alpha[j][++k] += wt*dyda[m];
				beta[j] += dy*wt;
			}
		}
		
		*chisq += dy*dy*sig2i;
	}

	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_vector(dyda,1,ma);
}


int gaussj(float **a, int n, float **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	float big,dum,pivinv,temp;

	indxc=ivector(1,n);
	indxr=ivector(1,n);
	ipiv=ivector(1,n);
	for (j=1;j<=n;j++) ipiv[j]=0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if (ipiv[j] != 1)
				for (k=1;k<=n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					} else if (ipiv[k] > 1)return (nrfail("gaussj: Singular Matrix-1"));
				}
		++(ipiv[icol]);
		if (irow != icol) {
			for (l=1;l<=n;l++) SWAP(a[irow][l],a[icol][l])
			for (l=1;l<=m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;
		if (a[icol][icol] == 0.0) return (nrfail("gaussj: Singular Matrix-2"));
		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=1;l<=n;l++) a[icol][l] *= pivinv;
		for (l=1;l<=m;l++) b[icol][l] *= pivinv;
		for (ll=1;ll<=n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=1;l<=n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=1;l<=m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n;l>=1;l--) {
		if (indxr[l] != indxc[l])
			for (k=1;k<=n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,1,n);
	free_ivector(indxr,1,n);
	free_ivector(indxc,1,n);
}


void covsrt(float **covar, int ma, int ia[], int mfit)
{
	int i,j,k;
	float temp;

	for (i=mfit+1;i<=ma;i++)
		for (j=1;j<=i;j++) covar[i][j]=covar[j][i]=0.0;
	k=mfit;
	for (j=ma;j>=1;j--) {
		if (ia[j]) {
			for (i=1;i<=ma;i++) SWAP(covar[i][k],covar[i][j])
			for (i=1;i<=ma;i++) SWAP(covar[k][i],covar[j][i])
			k--;
		}
	}
}


void fsin(float x, float a[], float *y, float dyda[], int na)

{

	int i;
	float arg;

	*y=0.0;

	for (i=1;i<=na-1;i+=3) { 
		arg=((x+a[i+2])*(PI/180));
		*y += a[i]+a[i+1]*(sin(arg));
		dyda[i]=1;
		dyda[i+1]=sin(arg);
		dyda[i+2]=a[i+1]*(cos(arg))*(PI/180);

	} 

}
