/* File gridmain.c 2:04 PM 11/6/97 */

#include "grid.h"

/** static function prototypes **/
static int isconsis(sg_struct *A, sg_struct *B);
static int isnew(sg_struct *A);

int main(int argc,char **argv,char **env){

	int nfiles;  /* No. of sinugrams to be processed */
	int ifile;   /* No. of next file to be read in */
	long im_size; /* No rows or columns in output image */
	int first=1;  /* Flags 1st time thru for-loop */
	float **sgram1, **sgram2, **image1, **image2; 
	sg_struct  sg_parm1,sg_parm2;
	grid_struct grid_parm;
	pswf_struct *pswf;

/**** Get parameters for gridding algorithm **/

	get_parm(argc,argv,&grid_parm);

	if(verbose) prn_gparm(&grid_parm);

/**** Set up to read the input data files ***/

	nfiles=data_setup(argc,argv);

/**** 	Read the data files two at a time, process them as complex data,
	and write out the image pairs that result. If the two sinugrams 
	are incompatible, or if only one remains, double up on the one 
	file. ***/

	for(ifile=0;ifile<nfiles;ifile++){
		sgram1=get_sgram(ifile,&sg_parm1);

		if(ifile+1==nfiles) sgram2=sgram1;
		else
		{
			sgram2=get_sgram(ifile+1,&sg_parm2);

		if(!isconsis(&sg_parm1,&sg_parm2)){
				rel_sgram(sgram2);
				sgram2=sgram1;
			}		
		}

	
		if(isnew(&sg_parm1))  /* Don't reinit if same params */
		{
			if(!first)
				recon_free(); 
			else
				first=0;

			recon_init(&grid_parm,&sg_parm1);
		}


		do_recon(sgram1,sgram2,&image1,&image2,&im_size);

		put_image(ifile,image1,im_size);
		if(sgram2!=sgram1){
			put_image(ifile+1,image2,im_size);
			ifile++;
		}

		free_matrix(image1);
		free_matrix(image2);

	} /* End for-loop over input files */


/**** Release allocated memory and quit  ***/

	recon_free();
	exit(0);

}	/****End main() ***/


static int isconsis(sg_struct *A, sg_struct *B){

	if(A->n_ang!=B->n_ang) return 0;
	if(A->n_det!=B->n_det) return 0;
	if(abs(A->center-B->center) > TOLERANCE) return 0;
	if(A->geom!=B->geom) return 0;
	if(A->geom==0 && A->angles!=B->angles)return 0;
	return 1;

}	/**** End isconsis() ***/


static int isnew(sg_struct *A){

	static int n_ang=0, n_det=0, geom=-1;
	static float center=0.0;
	static float* angles=NULL;

	if(
		A->n_ang==n_ang && 
		A->n_det==n_det && 
		abs(A->center-center)< TOLERANCE &&
		A->geom==geom &&
		(geom!=0 || A->angles==angles)
	)
		return 0;
	else{
		n_ang=A->n_ang;
		n_det=A->n_det;		
		center=A->center;
		geom=A->geom;
		if(geom==0) angles=A->angles;
		else angles=NULL;
		return 1;
	}
}   /**** End isnew() ***/




