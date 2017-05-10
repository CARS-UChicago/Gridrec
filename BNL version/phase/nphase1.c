#include <sys/types.h>
#include <sys/stat.h>
#include <stdio.h>
#include <fcntl.h>
#include <unistd.h>
#include <stdlib.h>
#include <bstring.h>
#include "PXLtoCDF.h"
#include "libhead.h"
#define THRESHOLD 100
int Threshold = THRESHOLD;
typedef enum TrueFalse { False, True } Flag;
Flag DelFlag = False;
Flag DeSpikeFlag = True;
char *inputFile = "Input";
FILE *in;
int arysize = 20;
#define MAXANGLES 2500
float anglelist[MAXANGLES];
char afname[] = "anglelist";
float **Norm;
int nwf = 0;
float ***buffer = NULL;
int WfX, WfY;
void writeslabs(float ***b, int x, int y, int z);
Flag validdata(unsigned short **b, int xmax, int ymax);
void initslicefiles(int nfiles, int nxvals);
void correctslicefiles(int nfiles, int nangles);
void despike(unsigned short **v, int nx, int ny);
/*
void bubblein(unsigned short d, unsigned short *q, unsigned short *a);
*/

void usage() {
    printf("Look at the writeup\n");
}
void setop(char c) {
    struct PXLHeader WfHd;

    switch(c) {
	case 'D':
	    DelFlag = True;
	    break;
	case 'f':
	    inputFile = optarg;
	    break;
	case 'A':
	    arysize = atoi(optarg);
	    break;
	case 'S':
	    DeSpikeFlag = False;
	    break;
	case 't':
	    Threshold = atoi(optarg);
	    break;
	case 'w':
	    {
		unsigned short **wf;
		int x, y;
		
		if (fexist(optarg)) {
		    wf = readPXL(optarg, &WfHd, &WfX, &WfY);
		    if (DeSpikeFlag) despike(wf, WfX, WfY);
		    if (nwf == 0) {
			Norm = (float**)calloc2D(WfY, WfX, sizeof(float));
		    }
		    for (y=0;y<WfY;y++) 
			for(x=0;x<WfX;x++){
			    Norm[y][x] += wf[y][x];
			}
		    nwf++;
		    free2D((void**)wf);
		} else {
		    printf("Whitefield file %s does not exist\n", optarg);
		    exit(1);
		}
		break;
	    }
	default:
	    {
		usage();
		exit(1);
	    }
    }
}

main(int argc, char **argv)
{
    char c;
    unsigned short **view;
    struct PXLHeader viewHd;
    char fileName[100];
    float fangle;
    int xmax, ymax;
    int fcount, err, nv;
    int oxmax = 0, oymax = 0;
    int x, y, angleindex, nvalidangles;
    int nviews;
    int afd;
    int inchar;
    Flag contin;
    char fnPrefix[100], fnSuffix[20];
    int fnNumb, fnNumbMax;
    float fnDelT;  
    
    while ((c = getopt(argc, argv, "DSf:A:w:t:")) != (char)EOF) 
	setop(c);
    if (nwf != 0) {
	for(x=0;x<WfX;x++)
	    for (y=0;y<WfY;y++) {
		Norm[y][x] /= nwf;
	    }
    } else {
	printf("Must have at least one whitefield for normalization\n");
	exit(1);
    }
afd = open("Norm",O_WRONLY,0777);
write(afd,&WfX,sizeof(int));
write(afd,&WfY,sizeof(int));
write(afd,Norm[0],WfX*WfY*sizeof(float));
close(afd);
    in = fopen(inputFile, "r");
    if (in == NULL) {
	perror("Cannot open Input file");
	exit(1);
    }
    fcount = 0;
    err = 0;
    contin = False;
    while (contin || (inchar = getc(in)) != EOF) {
	if (contin) {
	    if (fnNumb >= fnNumbMax) {
		contin = False;
		continue;
	    }
	    sprintf(fileName, "%s%03d", fnPrefix, ++fnNumb);
	    fangle = fangle + fnDelT;
	} else {
	    switch(inchar) {
		case '#':
		    if (fscanf(in, "%s %d %d %g %g\n", 
		    fnPrefix, &fnNumb, &fnNumbMax, 
		    &fangle, &fnDelT) != 5) {
			fprintf(stderr, "Input line # format error\n");
			err++;
			break;
		    }
		    contin = True;
		    sprintf(fileName, "%s%03d", fnPrefix, fnNumb);
		    break;
		case '-':
		    {
			char optargx[100];
			
			fscanf(in, "%c%s\n", &inchar, optargx);
			optarg = optargx;
			setop(inchar);
			break;
		    }
		default:
		    ungetc(inchar, in);
		    if (fscanf(in, "%s %g\n", fileName, &fangle) != 2){
			fprintf(stderr, "Input file format error\n");
			err++;
			break;
		    }
	    }
	}
	if (fcount >= MAXANGLES)  {
	    fprintf(stderr, "Too many input files, increase MAXANGLES and recompile.\n");
	    fprintf(stderr, " Current value is %s\n", fcount);
	    exit(10);
	}
	if (fexist(fileName)) {
	    fcount++;
	}
	else {
	    fprintf(stderr, "File %s does not exist\n", fileName);
	    err++;
	}
    }
    if (err > 0) {
	fprintf(stderr, "Aborting because of errors\n");
	exit(1);
    }
    printf("Found %d files to process\n", fcount);
    rewind(in);
    angleindex = 0;
    nvalidangles = 0;
    contin = False;
    for (nv = 0; nv<fcount; nv++) {
	if (contin) {
	    if (fnNumb >= fnNumbMax) {
		contin = False;
		continue;
	    }
	    sprintf(fileName, "%s%03d", fnPrefix, ++fnNumb);
	    fangle = fangle + fnDelT;	    
	} else {
	    switch(inchar = getc(in)) {
		case '#':
		    if (fscanf(in, "%s %d %d %g %g\n", 
		    fnPrefix, &fnNumb, &fnNumbMax, 
		    &fangle, &fnDelT) != 5) {
			fprintf(stderr, "Input line # format error\n");
			err++;
			break;
		    }
		    contin = True;
		    sprintf(fileName, "%s%03d", fnPrefix, fnNumb);
		    break;
		case '-':
		    {
			while(getc(in) != '\n');
		    }
		default:
		    ungetc(inchar, in);
		    fscanf(in, "%s %g\n", fileName, &fangle);
	    }
	}
	printf("Processing %s view from angle %g\n", fileName, fangle);
	view = readPXL(fileName, &viewHd, &xmax, &ymax);
	if (validdata(view, xmax, ymax) == False) {
	    /* Skip this view */
	    fprintf(stderr, "View %d (%s) Rejected for invalid data\n", 
		nv, fileName);
	    free2D((void **)view);
	    continue;
	}
	if (DeSpikeFlag) despike(view, xmax, ymax);
	if (oxmax == 0) {
	    oxmax = xmax;
	    oymax = ymax;
	    printf("View is %dx%d\n", xmax, ymax);
	    initslicefiles(ymax, xmax);
	}
	else {
	    if ((xmax != oxmax) | (ymax != oymax)){
		fprintf(stderr, "Image sizes are not consistent:\n"
                    "\tprevious: xsize ysize = %d %d\n"
                    "\tcurrent : xsize ysize = %d %d\n",
                    oxmax, oymax, xmax, ymax);
		exit(2);
	    }
	}
	if (buffer == NULL) {
	    float size;
	    
	    nviews = (arysize * 1.e6)/(xmax * ymax * sizeof(float));
	    size = (float)nviews * xmax * ymax * sizeof(float);
	    printf("Allocating buffer for %d views,  total size %g bytes\n", 
	       nviews, size);
	    buffer = (float ***)calloc3D(ymax, nviews, xmax, sizeof(float));
	}
	anglelist[nvalidangles] = fangle;
	if (angleindex >= nviews) {
	    writeslabs(buffer, xmax, ymax, nviews);
	    angleindex = 0;
	}
	for (y=0;y<ymax;y++) {
	    for (x=0;x<xmax;x++) {
		buffer[y][angleindex][x] = (float)view[y][x] / Norm[y][x];
	    }
	}
	angleindex++;nvalidangles++;
	free2D((void **)view);
	if (DelFlag == True) unlink(fileName);
    }
    /* Yes, the only way angleindex can be zero at this point is if there was no valid input! */
    if (angleindex != 0) writeslabs(buffer, xmax, ymax, angleindex);
    correctslicefiles(ymax, nvalidangles);
    afd = open(afname, O_WRONLY|O_CREAT, 0777);
    write(afd, &nvalidangles, sizeof(int));
    write(afd, anglelist, nvalidangles*sizeof(float));
    close(afd);
    return(0);
}

void writeslabs(float ***b, int xmax, int ymax, int zmax)
{
    int slice;
    char fname[10];
    int fd;
    
    printf("Writing slices\n");
    for(slice=0;slice<ymax;slice++){
	sprintf(fname, "slice%03d", slice);
	fd = open(fname, O_WRONLY|O_APPEND|O_CREAT, 0777);
	write(fd, b[slice][0], xmax * zmax * sizeof(float));
	close(fd);
    }
}

void initslicefiles(int nfiles, int xmax)
{
    int i, fd;
    char fname[10];
    
    for (i=0; i<nfiles; i++) {
	sprintf(fname, "slice%03d", i);
	fd = open(fname, O_WRONLY|O_CREAT|O_TRUNC, 0777);
	/* Here we write the true x dimension, but a dummy number
	 * of angles (because we don't know how many views will be
	 * discarded).  This is fixed up later.
	 */
	write (fd, &xmax, sizeof(int));
	write (fd, &xmax, sizeof(int));
	close(fd);
    }
}

void correctslicefiles(int nfiles, int nangles)
{
    int fd, i;
    char fname[10];
    
    for (i=0; i<nfiles; i++) {
	sprintf(fname, "slice%03d", i);
	fd = open(fname, O_WRONLY, 0777);
	lseek(fd, sizeof(int), SEEK_SET);
	write(fd, &nangles, sizeof(int));
	close(fd);
    }
    
}
Flag validdata(unsigned short **b, int xmax, int ymax)
{
    int bcnt;
    int i, j;
    
    bcnt = 0;
    for(j=0;j<ymax;j++)
	for(i=0;i<xmax;i++)
	    if (b[j][i] == 4095) 
		bcnt++;
    if (bcnt < Threshold)
	return(True);
    else
	return(False);
}
/* bubblein is a bubble sort that adds an entry to a sorted list
 * and sorts it into place to keep the list sorted in increasing 
 * order.
 */
#define bubblein(v) \
    *q = v;\
    p=q;\
    do { \
	if (*p < *(p - 1)){\
	    tmp = *p;\
	    *p = *(p - 1);\
	    *(p - 1) = tmp;\
	}\
	else break;\
    } while (--p > a);\
    q++;

void despike(unsigned short **v, int nx, int ny)
{
    /* This routine applies a 3x3 median filter to the array v with
     * dimensions nx by ny
     */
    int x, y;
    unsigned short a[8];
    unsigned short *p, *q, tmp;
    unsigned short **w;
    
    w = (unsigned short **)calloc2D(ny, nx, sizeof(unsigned short));
    for (y=0; y<ny; y++)
	for (x=0; x<nx; x++) {
	    if ((x==0)|(x==(nx-1))|(y==0)|(y==(ny-1))) {
		w[y][x] = v[y][x];
	    } else {
		a[0] = v[y-1][x-1];
		q = &(a[1]);
		bubblein(v[y][x-1]);
		bubblein(v[y+1][x-1]);
		bubblein(v[y-1][x]);
		bubblein(v[y][x]);
		bubblein(v[y+1][x]);
		bubblein(v[y-1][x+1]);
		bubblein(v[y][x+1]);
		tmp = v[y+1][x+1];
		/* This avoids the last sort by seeing if the last entry
		 * will fall in the middle or in the outer pieces of the
		 * sorted array
		 */
		if ( tmp <= a[3] ) tmp = a[3];
		else if (tmp >= a[4]) tmp = a[4];
		w[y][x] = tmp;
	    }
	}
    bcopy(w[0], v[0], nx * ny * sizeof(unsigned short));
    free2D((void**)w);
}
/*
void bubblein(unsigned short d, unsigned short *q, unsigned short *a)
{
    unsigned short *p = q;
    unsigned short tmp;
    
    *p = d;
    do {
	if (*p < *(p - 1)){
	    tmp = *p;
	    *p = *(p - 1);
	    *(p - 1) = tmp;
	}
	else break;
    } while (--p > a);
}
*/
