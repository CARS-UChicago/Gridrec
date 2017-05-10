#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "PMIStoCDF.h"

unsigned short **
readPMIS (char *PMISfile, struct PMISHeader *Header, int *xsize, int *ysize)
{
    FILE *fp;
    size_t size;
    struct tm caltime;
    int arraysize;
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

    buffer = (unsigned short **)calloc2D(Header->nYsiz, Header->nXsiz, sizeof(short));
    arraysize = Header->nYsiz * Header->nXsiz * sizeof(short);
    size = fread (buffer[0], sizeof d, arraysize, fp);
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

