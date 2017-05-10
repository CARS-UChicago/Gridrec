#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "PXLtoCDF.h"

unsigned short **
readPXL (char *PXLfile, struct PXLHeader *Header, int *xsize, int *ysize)
{
    FILE *fp;
    size_t size;
    int arraysize,cnt;
    unsigned short **buffer = NULL;
    unsigned short d;
    int x, y;
    
    fp = fopen (PXLfile, "r");
    
    size = fread (&Header->hsiz, sizeof Header->hsiz, 1, fp);
    size = fread (&Header->text, sizeof Header->text, 1, fp);
    size = fread (&Header->startx, sizeof Header->startx, 1, fp);
    size = fread (&Header->starty, sizeof Header->starty, 1, fp);
    size = fread (&Header->totalx, sizeof Header->totalx, 1, fp);
    size = fread (&Header->totaly, sizeof Header->totaly, 1, fp);
    size = fread (&Header->bpp, sizeof Header->bpp, 1, fp);
    size = fread (&Header->exp_t, sizeof Header->exp_t, 1, fp);
    size = fread (&Header->exp_n, sizeof Header->exp_n, 1, fp);
    size = fread (&Header->spc, sizeof Header->spc, 1, fp);
    size = fread (&Header->units, sizeof Header->units, 1, fp);
    size = fread (&Header->date, sizeof Header->date, 1, fp);
    size = fread (&Header->drk, sizeof Header->drk, 1, fp);
    size = fread (&Header->rad, sizeof Header->rad, 1, fp);
    size = fread (&Header->geom, sizeof Header->geom, 1, fp);
    size = fread (&Header->src, sizeof Header->src, 1, fp);
    size = fread (&Header->opt, sizeof Header->opt, 1, fp);
    size = fread (&Header->pos, sizeof Header->pos, 1, fp);
    size = fread (&Header->expt, sizeof Header->expt, 1, fp);
    size = fread (&Header->label, sizeof Header->label, 1, fp);
    size = fread (&Header->pixel_type, sizeof Header->pixel_type, 1, fp);
    size = fread (&Header->is_rgb, sizeof Header->is_rgb, 1, fp);
    size = fread (&Header->section_type, sizeof Header->section_type, 1, fp);
    size = fread (&Header->mosaic_x, sizeof Header->mosaic_x, 1, fp);
    size = fread (&Header->mosaic_y, sizeof Header->mosaic_y, 1, fp);
    size = fread (&Header->nbands, sizeof Header->nbands, 1, fp);
    size = fread (&Header->nsections, sizeof Header->nsections, 1, fp);
    size = fread (&Header->version_magic, sizeof Header->version_magic, 1, fp);
    size = fread (&Header->zspc, sizeof Header->zspc, 1, fp);
    size = fread (&Header->xspc, sizeof Header->xspc, 1, fp);
    size = fread (&Header->yspc, sizeof Header->yspc, 1, fp);
    size = fread (&Header->magic, sizeof Header->magic, 1, fp);

/*      printf ("startx = %d\n", Header->starty);
        printf ("startx = %d\n", Header->startx);
        printf ("totalx = %d\n", Header->totalx);
        printf ("totaly = %d\n", Header->totaly);
        printf ("date  = %s", Header->date);
        printf ("\n");
*/

    buffer = malloc_matrix_us(Header->totaly, Header->totalx);
    arraysize = Header->totaly * Header->totalx;
    size = fread (buffer[0], sizeof d, arraysize, fp);

/*    for (x = 0; x < Header->totalx; x++) {
        cnt = 1;
	for (y = 0; y < Header->totaly; y++) {
                 printf (" %i ,", buffer[x][y]);
                 if (cnt%10 == 0 || cnt%Header->totalx == 0)
                 printf("\n");
                 cnt++;
        }
    }
*/

    *xsize = Header->totalx;
    *ysize = Header->totaly;
    fclose (fp);
    return buffer;

}

