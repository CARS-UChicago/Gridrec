/* We assume that the whitefield normalizations are all completed before 
 * this program is invoked.
 * Call:
 *  phase2 -switches slicefile
 * 
 */
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <netcdf.h>
#include "PXLtoCDF.h"
#include "libhead.h"
typedef enum TrueFalse { False, True } Flag;
Flag DelFlag, BoxFlag, CenterFlag, GraphFlag, Verbose;
extern int opting, opterr, optopt;
extern char *optarg;
int nangles;
float *angles;
int ncols;
float **data;
char *GraphFileName;
int Airvalues = 0;
int Ndiscard = 0;
float SliceCenter;
void readfiles(char *name, int *nangles, int *ncols, float ***data, float **angles);
void usage();
void convert_to_sinogram(float **data, int nc, int na);
void suppress_rings(float **data, int nc, int na);
void center_sinogram(float **data, int nc, int na);
void write_data(char *name, float **data, int nc, int na);
int cgfit(float x[],float y[],int ndata,float a[]);
float ** trim_data(float **data, int ncols, int nangles, int lcol, int hcol);

void setopt(char c) {
    switch (c) {
	case 'D':
	    DelFlag = True;
	    break;
	case 'B':
	    BoxFlag = False;
	    break;
	case 'C':
	    CenterFlag = False;
	    break;
	case 'g':
	    GraphFlag = True;
	    GraphFileName = optarg;
	    break;
	case 'v':
	    Verbose = True;
	    break;
	case 'a':
	    Airvalues = atoi(optarg);
	    break;
	case 'n':
	    Ndiscard = atoi(optarg);
	    break;
	default:
	    usage();
	    exit(1);
    }
}
main(int argc, char **argv)
{
    char c;
    float **t;
    int lcol, hcol;
    
    BoxFlag = CenterFlag = True;
    while ((c = getopt(argc, argv, "DBCg:vn:a:")) != (char)EOF)
	setopt(c);
    readfiles(argv[optind], &nangles, &ncols, &data, &angles);
    /* Until get better scheme */
    lcol = Ndiscard;
    hcol = ncols - Ndiscard;
    if (((hcol - lcol) % 2) == 0) hcol--;
    t = trim_data(data, ncols, nangles, lcol, hcol);
    free2D((void**)data);
    data = t;
    ncols = hcol - lcol;
    convert_to_sinogram(data, ncols, nangles);
    if (BoxFlag) suppress_rings(data, ncols, nangles);
    if (CenterFlag) center_sinogram(data, ncols, nangles);
    write_data(argv[optind], data, ncols, nangles);
    if (DelFlag) {
	unlink(argv[optind]);
    }
    return(0);
}

void readfiles(char *name, int *nangles, int *ncols, float ***data,
     float **angles)
{
    int fd;
    int tmp[2];
    int ret;
    
    fd = open("anglelist", O_RDONLY);
    read(fd, nangles, sizeof(int));
    *angles = (float*)calloc(*nangles, sizeof(float));
    read(fd, *angles, (*nangles) * sizeof(float));
    close(fd);
    fd = open(name, O_RDONLY);
    read(fd, &tmp, 2*sizeof(int));
    if (tmp[1] != *nangles) {
	fprintf(stderr, "angle file length (%d) does not correspond with\n", *nangles);
	fprintf(stderr, "slice angles dimension (%d)\n", tmp[1]);
	exit(2);
    }
    *ncols = tmp[0];
    *data = (float**)calloc2D(tmp[1], tmp[0], sizeof(float));
    
    ret = read(fd, (*data)[0], tmp[0] * tmp[1] * sizeof(float));
    if (ret != tmp[0] * tmp[1] * sizeof(float)) {
		printf("Ret err\n");
	    }
    /*
    for(y=0;y<(*nangles);y++)
	for (x=0;x<(*ncols);x++) {
	    ret = read(fd, &((*data)[y][x]), sizeof(float));
	    if (ret != sizeof(float)) {
		printf("Ret err\n");
	    }
	}
    */
}

void usage()
{
    fprintf(stderr, "Error in program call\n");
}

float ** trim_data(float **data, int ncols, int nangles, int lcol, int hcol)
{
    float **t;
    int x, x1, ang;
    
    t = (float**)calloc2D(nangles, hcol - lcol, sizeof(float));
    for (ang = 0; ang < nangles; ang++) {
	x1 = 0;
	for (x = lcol; x < hcol; x++) {
	    t[ang][x1++] = data[ang][x];
	}

    }
}
void convert_to_sinogram(float **data, int nc, int na)
{
    int x, y, i;
    float air;
    
    for (y=0; y<na; y++) {
	if (Airvalues > 0) {
	    air = 0;
	    for(i=0; i<Airvalues; i++){
		air += data[y][i] + data[y][nc-1-i];
	    }
	    air /= 2 * Airvalues;	    
	} 
	else {
	   air = 1; 
	}
	for (x=0; x<nc; x++)
	{
	    data[y][x] = - log(data[y][x]/air);
	}
    }
}
void suppress_rings(float **data, int nc, int na)
{
    int x, angle;
    float *avg = NULL;
    float *diff = NULL;
    float *smooth = NULL;
    /* The smoothing width must be odd */
    int halfsmoothwidth = 5;
    /*int smoothwidth = 2*halfsmoothwidth + 1;*/

    avg = calloc(nc, sizeof(float));
    diff = calloc(nc, sizeof(float));
    smooth = calloc(nc, sizeof(float));

    /* For each column, sum up all of the angles, and average */
    for (angle = 0; angle < na; angle++) {
	for (x = 0; x < nc; x++) {
	    avg[x] += data[angle][x];
	}
    }

    for (x = 0; x < nc; x++) {
	avg[x] /= na;
    }

    /* Smooth the average */
    for (x = 0; x < nc; x++) {
       int ismooth, width = 0;
       float sum = 0;
       for (ismooth = -halfsmoothwidth; ismooth <= halfsmoothwidth; ismooth++) {
	   int index = ismooth + x;
	   if (index >= 0 && index < nc) {
	      width++;
	      sum += avg[index];
	   }
       }
       smooth[x] = sum/width;
    }

    /* Make the difference array */
    for (x = 0; x < nc; x++) {
	diff[x] = avg[x] - smooth[x];
    }

    /* Subtract the difference from the data set */
    for (angle = 0; angle < na; angle++) {
	for (x = 0; x < nc; x++) {
	    data[angle][x] -= diff[x];
	}
    }
    free(smooth);
    free(diff);
    free(avg);
}
void  covsrt(float **covar, int ma, int ia[], int mfit);
/*
void  ffsin(float x, float a[], float *y, float dyda[], int na);
*/
int  gaussj(float **a, int n, float **b, int m);
void  mrqcof(float x[], float y[], float sig[], int ndata, float a[],
	int ia[], int ma, float **alpha, float beta[], float
	*chisq, void (*funcs)(float,float [],float *,float [],int));
int  mrqmin(float x[], float y[], float sig[], int ndata, float a[],
	int ia[],int ma, float **covar, float **alpha, 
	float *chisq, void (*funcs)(float, float [], float *,float [],
	int),float *alamda);


float *vector(long nl, long nh);
float **matrix(long nrl, long nrh, long ncl, long nch);
int *ivector(long nl, long nh);
void free_vector(float *v, long nl, long nh);
void free_ivector(int *v, long nl, long nh);
void free_matrix(float **m,long nrl, long nrh, long ncl, long nch);
void nrerror(char error_text[]);
int nrfail(char error_text[]);
#define PI 3.1415926536


void center_sinogram(float **data, int nc, int na)
{
    float *cog = NULL;
    float *shiftvals = NULL;
    float *theta = NULL;
    float avsin = 0;
    int x, angle;
    float shift = 0;
    int ishift;
    float fl, fu;
    float xshift;	    
    float halfsize;
    float a[4]=
	{0.0,150.0,15.0,0.0};

    cog = (float*) calloc(nangles, sizeof(float));
    shiftvals = (float*) calloc(nc, sizeof(float));
    theta = (float*) calloc(nangles, sizeof(float));

 
    
    a[1] = halfsize = nc/2;
    do {       
    for (angle = 0;angle < na; angle++) {
	float sum = 0, sumx = 0;
    
	for (x = 0; x < nc; x++) {
	  sum += data[angle][x];
	  sumx += (x) * data[angle][x];
	}
	cog[angle] = sumx / sum;
	theta[angle] = angles[angle];
	/*
	 * ?????????????????????????????
	 */
	if(cog[angle] < 0)
	  cog[angle] = -cog[angle];
	avsin +=cog[angle];
    }
    avsin /= na;
    if (cgfit(theta,cog,na,a)==1){
       printf("Centering Fit did not converge\n");
       return;
    } else {
	if (GraphFlag) {
	    FILE *xvfile;
	    xvfile = fopen(GraphFileName, "w");
	    for (angle=0;angle<na; angle++)
	    fprintf(xvfile, "%g %g %g\n", theta[angle], cog[angle], 
		a[1] + a[2] * sin((double)theta[angle] * (double)(PI/180)
		    + atan((double)a[3])));
	    fclose(xvfile);
	}
	SliceCenter = a[1] + 1.0;   /* For 1 based indexing in RECLBL */
      /* Shift the center of gravity */
      /*  
       * The following centering code is not needed and should not be used if
       * we are letting RECLBL take care of it.
       */
       /*
	shift = (halfsize - a[1]);
	for (angle = 0; angle < na; angle++) {
	    for (x = 0; x < nc; x++) {
		shiftvals[x] = data[angle][x];
	    }
	    for (x = 0; x < nc; x++) {
		xshift = x - shift;
		ishift = xshift;
		if (ishift < 0 | ishift + 1 >= nc) continue;
		fu = xshift - ishift;
		fl = 1.0 - fu;
		data[angle][x] = fl * shiftvals[ishift] + fu * shiftvals[ishift+1];
	    }
	}
	*/
    }
    }while (fabsf(shift) > 0.5);
    free(shiftvals);
    free(theta);
    free(cog);
}

void write_data(char *inname, float **data, int nc, int na)
{
    char name[256];
    
    strcpy(name, inname);
    strcat(name, ".nc");
    saveasnetCDF(name, data, nc, na);
}

void
saveasnetCDF (char *name, float **m, int x, int y)
{
    char CDFname[256];
    int ncid;
    int xid, angleid, centerid;
    int sinogram_dimids[2], sinogram_id;
    int anglevarid;
    long start[2], count[2];

    strcpy (CDFname, name);
    printf ("Saving slice data to file %s.\n", CDFname);
    
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
    anglevarid = ncvardef(ncid, "angle", NC_FLOAT, 1, &angleid);
    centerid = ncvardef(ncid, "center", NC_FLOAT, 0, NULL);
    
    /* Leave defined mode */
    ncendef (ncid);

    /* Write the data */
    start[0] = 0;
    start[1] = 0;
    count[0] = y;
    count[1] = x;
    ncvarput (ncid, anglevarid, start, count, (void *)angles);
    ncvarput (ncid, sinogram_id, start, count, (void *)&m[start[0]][start[1]]);
    ncvarput1 (ncid, centerid, NULL, &SliceCenter);
    
    /* Close the netCDF file */
    ncclose (ncid);


    return;
}

/*NR Least Squares MRQ fitting routine*/

#define MA 3
#define PI 3.1415926536
#define EPS 1.0e-6
#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}

void ffsin(float x, float a[], float *y, float dyda[], int na)
{
    double arga, tmp;

    if (na != 3){
	printf("ffsin called with %d parameters,  expecting only 3!\n", 
	    na);
	exit(10);
    }
    arga = (double)x*(PI/180) + atan((double)a[3]);
    tmp = sin(arga);
    *y = a[1]+a[2]*tmp;
    dyda[1]=1;
    dyda[2]=tmp;
    dyda[3]=a[2] * cos(arga) / (1.0 + a[3]*a[3]);
}

#define MAXITER 30
int cgfit(float x[],float y[],int ndata,float a[])
{
	int i,*ia,iter,itst,k,mfit=MA,angle;
	float alamda,chisq,ochisq,*sig,**covar,**alpha, olamda;
	
	ia=ivector(1,MA);
	sig=vector(1,ndata);
	covar=matrix(1,MA,1,MA);
	alpha=matrix(1,MA,1,MA);

	for (angle = 0; angle < ndata; angle++) {
	    sig[angle]=sqrt(y[angle]);
	}

	for (i=1;i<=mfit;i++) ia[i]=1;
	alamda = -1;
	if(mrqmin(x,y,sig,ndata,a,ia,MA,covar,alpha,&chisq,ffsin,&alamda)==1) 
	    return 1;
	k=1;
	itst=0;

	do {
	  if (Verbose) {
	      printf("Iteration # %2d chi-squared: %10.4f lamda: %9.2e\n",
		k, chisq, alamda);
	      printf("Coefficients %g %g %g\n", a[1], a[2], a[3]);
	  }
	  /*
	  {
	    FILE *xvfile;int angle;
	    static int count = 0;
	    char fn[100];
	    
	    sprintf(fn, "cogtrace%d.xmgr", count++);
	    xvfile = fopen(fn, "w");
	    for (angle=0;angle<ndata; angle++)
	    fprintf(xvfile, "%g %g %g\n", x[angle], y[angle], 
		a[1] + a[2] * 
		    sin((double)(x[angle] * (double)(PI/180) + atan((double)a[3]))));
	    fclose(xvfile);
	  }
	  */
	  ochisq=chisq;
	  olamda = alamda;
	  mrqmin(x,y,sig,ndata,a,ia,MA,covar,alpha,&chisq,ffsin,&alamda);
	  if (alamda > olamda) 
	    itst = 0;
	  else
	    if ((ochisq - chisq) < 0.01)
		itst++;
	} while ((itst < 2) & (k++ < MAXITER));
	if (k < MAXITER){
	    alamda=0.0;
	    mrqmin(x,y,sig,ndata,a,ia,MA,covar,alpha,&chisq,ffsin,&alamda);
	    
	    printf("Centering Coefficients %g +- %.3f %g +- %.3f %g +- %.3f\n",
		a[1], sqrt(covar[1][1]), 
		a[2], sqrt(covar[2][2]),
		a[3], sqrt(covar[3][3]));
	} else {
	    printf("\nConvergence not obtained after MAXITER iterations\n");
	    return 1;	    
	}
	free_matrix(alpha,1,MA,1,MA);
	free_matrix(covar,1,MA,1,MA);
	free_vector(sig,1,ndata);
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
		mrqcof(x,y,sig,ndata,a,ia,ma,alpha,beta,chisq,ffsin);
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
		return 0;
	}
	for (j=0,l=1;l<=ma;l++)
		if (ia[l]) atry[l]=a[l]+da[++j];
	mrqcof(x,y,sig,ndata,atry,ia,ma,covar,da,chisq,ffsin);
	if (*chisq <= ochisq * (1.0 + EPS)) {
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
		(*ffsin)(x[i],a,&ymod,dyda,ma);
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
	*chisq /= ndata;
	for (j=2;j<=mfit;j++)
		for (k=1;k<j;k++) alpha[k][j]=alpha[j][k];
	free_vector(dyda,1,ma);
}


gaussj(float **a, int n, float **b, int m)
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
	return(0);
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
