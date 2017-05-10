#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/stat.h>
#include "PMIStoCDF.h"

static int verbose = 0;

main(int argc, char**argv)
{
    unsigned short **buffer;
    unsigned short **newbuf;
    struct PMISHeader head;
    struct PMISHeader* Header = &head;
    long newx, newy;
    int x,y;
    int i, j;
    struct tm caltime;

    buffer = readPMIS(argv[1], Header, &x, &y);
    newx = (x%2 == 0) ? (x-1)/2:x/2;
    newy = (y%2 == 0) ? (y-1)/2:y/2;
    newbuf = malloc_matrix_us(newy, newx);
    for (j = 0; j<newy; j++){
	for(i=0; i<newx;i++) {
	    newbuf[j][i] = (buffer[j+j][i+i] + buffer[j+j+1][i+i] +
			    buffer[j+j][i+i+1] + buffer[j+j+1][i+i+1]) >> 2;
	}
    }
    Header->nXsiz = newx;
    Header->nYsiz = newy;
    Header->nXbin *= 2;
    Header->nYbin *= 2;
    strcat(Header->szComment, " Rebinned by 2x2");
    writePMIS(argv[2], Header, newbuf);
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
    size = fread (&Header->tModified, sizeof Header->tModified, 1, fp);
    Header->tModified = byteswap(Header->tModified, 4);
    size = fread (&Header->nGain, sizeof Header->nGain, 1, fp);
    Header->nGain = BYTESWAP2(Header->nGain);
    size = fread (&Header->nImages, sizeof Header->nImages, 1, fp);
    Header->nImages = BYTESWAP2(Header->nImages);

    /* The data size is the x size times the y size times 2 for the number of bytes. */
    arraysize = Header->nXsiz * Header->nYsiz;
    datasize = arraysize * sizeof (unsigned short);
    /* The data buffer must be freed outside of this function. */
    /* buffer = malloc (datasize); */
    buffer = malloc_matrix_us(Header->nYsiz, Header->nXsiz);

    size = fread (buffer[0], sizeof d, arraysize, fp);
    /* for (i = 0; i < arraysize; i++) buffer[i] = BYTESWAP2(buffer[i]); */
    for (y = 0; y < Header->nYsiz; y++) {
        for (x = 0; x < Header->nXsiz; x++) {
            buffer[y][x] = BYTESWAP2(buffer[y][x]);
        }
    }
    *xsize = Header->nXsiz;
    *ysize = Header->nYsiz;
    fclose (fp);
    return buffer;
}

void
writePMIS(char *PMISfile, struct PMISHeader *Header, unsigned short **buffer)
{
    FILE *fp;
    int arraysize, x, y;
    
     fp = fopen (PMISfile, "w");
     for (y = 0; y < Header->nYsiz; y++) {
        for (x = 0; x < Header->nXsiz; x++) {
            buffer[y][x] = BYTESWAP2(buffer[y][x]);
        }
     }
     arraysize = Header->nXsiz * Header->nYsiz;
     Header->nHeadSize = BYTESWAP2(Header->nHeadSize);
     Header->nVersion = BYTESWAP2(Header->nVersion);
     Header->nXorg = BYTESWAP2(Header->nXorg);
     Header->nYorg = BYTESWAP2(Header->nYorg);
     Header->nXsiz = BYTESWAP2(Header->nXsiz);
     Header->nYsiz = BYTESWAP2(Header->nYsiz);
     Header->nXbin = BYTESWAP2(Header->nXbin);
     Header->nYbin = BYTESWAP2(Header->nYbin);
     Header->tCreated = byteswap(Header->tCreated, 4);
     Header->tModified = byteswap(Header->tModified, 4);
     Header->nGain = BYTESWAP2(Header->nGain);
     Header->nImages = BYTESWAP2(Header->nImages);
     fwrite(Header,  sizeof(struct PMISHeader), 1, fp);
     fwrite(buffer[0],  sizeof(unsigned short), arraysize, fp);
     fclose(fp);
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
