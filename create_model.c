   /* test code for coordinate transformation */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <Python/Python.h>

#include "fitsio.h"
#include "make_vis.h"
#include "constants.h"

void dump_vis(struct model_fft *vis, int npix){
  
  int i;
  for(i=0;i<npix;i++)
    fprintf(stdout,"%g \t %g \t %g \t %g\n",vis->u[i], vis->v[i], vis->real[i],vis->imag[i]);
  return;
}

 /* helper function to return the z coordinate of the z'=0 (i.e., disk)plane */

void define_coordinate_transform(double **combined_transform, double inc, double pa){
  double **inc_transform;
  double **pa_transform;  /* matrices that hold relvant coordinate transforms*/
  int i,j,k;

  inc_transform =(double **)malloc(sizeof(double*)*3);
  pa_transform =(double **)malloc(sizeof(double*)*3);

  for(i = 0;i<3;i++){
      inc_transform[i] = (double*) malloc(sizeof(double)*3); 
      pa_transform[i] = (double*) malloc(sizeof(double)*3); 
      combined_transform[i] = (double*) malloc(sizeof(double)*3); 
      for(j=0;j<3;j++){
	inc_transform[i][j]=0;
	pa_transform[i][j] =0;
	combined_transform[i][j]=0;
      }
  }

  inc_transform[0][0] = (pa_transform[2][2] = 1);
  inc_transform[1][1] = (inc_transform[2][2] = cos(inc*DTOR));
  inc_transform[1][2] = -1*(inc_transform[2][1] = -sin(inc*DTOR));
  pa_transform[0][0] = (pa_transform[1][1] = cos(pa*DTOR));
  pa_transform[0][1] = -1*(pa_transform[1][0] = -sin(pa*DTOR));
  
  for(i = 0; i<3;i++)
    for(j = 0; j<3;j++)
      for(k = 0;k<3;k++)
	combined_transform[i][j] += inc_transform[i][k]*pa_transform[k][j];

  for(i = 0; i<3;i++){
    free(inc_transform[i]);
    free(pa_transform[i]);
  }
  if(DEBUG){
    for(i = 0; i<3;i++){
      for(j = 0; j<3;j++)
	fprintf(stderr,"%lf ",combined_transform[i][j]);
      fprintf(stderr,"\n");
    } 
  }

  free(inc_transform);
  free(pa_transform);
  return;
}

void pad_image(double *image, double *padded_image,int pix, int padfactor){
  int i,j ;
  for(i = 0; i < pix*pix*padfactor*padfactor;i++){
    padded_image[i]=0;
    int x = (i % (pix * padfactor));
    int y = (i / (pix * padfactor));
    
    if((x <=  pix/2) && (y <= pix/2)){
      padded_image[i]  = image[x + y*pix];
    }
  
  
    else if((x <=  pix/2) && (y > pix*padfactor - pix/2)){
      padded_image[i]  = image[x + (y - pix*(padfactor-1))  * pix];
    }
    
    else if((x > pix*padfactor - pix/2) && (y < pix/2)){
      padded_image[i]  = image[x-pix*(padfactor-1) + y * pix];
    }

    else if((x > pix*padfactor - pix/2) && (y > pix*padfactor - pix/2)){
      padded_image[i]  = image[x-pix*(padfactor-1) + (y - pix*(padfactor-1))  * pix];
    }
    else
      padded_image[i] = 0;

  }
  return;
}

void resort_image(double *array, int npix){
  int i;
  for(i = 0; i < npix*npix/2;i++){ /*only go over half the array, dummy */
    /* convert from "natural" order to fft order */    
    double tempim;
    int inew;
    int old_xpix  =   (i % npix - npix/2);
    int old_ypix  =   i/npix - npix/2;
    int new_xpix  =  (old_xpix >= 0) ? old_xpix : old_xpix + npix; /*maybe? */
    int new_ypix  =  (old_ypix >= 0) ? old_ypix : old_ypix + npix; /*maybe? */
    inew = new_xpix  + new_ypix * npix;

    tempim            = array[inew]; 
    array[inew] = array[i];
    array[i]    = tempim;
  }

  return;
}


double return_zconst_plane(double x, double y, double zval, double inc, double PA){  
  double cPA = cos(PA*DTOR), sPA = sin(PA*DTOR), cinc = cos(inc*DTOR), sinc = sin(inc*DTOR);
  double det  = cPA*cPA*cinc + sPA*sPA*cinc;
  double yval = 1./det * (( -1*sPA * cinc * (x - sinc*sPA*zval)) + cPA * (y + sinc*cPA*zval)); 
  return sinc * yval + cinc * zval; 
}

void create_grid(double *x, double *y, double *z, int npix, double inc, double pa){
  int i;
  for(i = 0; i < npix*npix*NZ;i++){
    int zind = (NZ/2 - i /(npix*npix));
    double zmax = 100*AU;
    double zmin = 0.0001*AU;
    double pix_scale = IMSIZE/npix;
    double zval = zind > 0 ? pow(10,(log10(zmax) - log10(zmin))/(NZ/2)*zind + log10(zmin)): 
      -1*pow(10,(log10(zmax) - log10(zmin))/(NZ/2)*(abs(zind)) + log10(zmin));  
    x[i] = ((i % npix) - npix/2)*DIST*pix_scale;  
    y[i] = (((i / npix)%npix - npix/2))*DIST*pix_scale;  
    z[i] = return_zconst_plane(x[i],y[i],zval,inc,pa); /*currently only the disk plane */
  }

}

void initialize_image(double *im, int npix){
  int i;
  for(i = 0; i < npix*npix;i++)
    im[i] = 0;
  return;
}
void compute_temp_tau(double *x, double *y, double *z, double **combined_transform, 
		      double *params, double *temp, 
		      double *tau, int npix, double subpixels){
  int i,j,k,l;

  double model_x=0;
  double model_y=0; 
  double model_z=0;
  double model_r;
  double model_phi;
  double h, d_tau, alpha=0;

  double sigma_0   = params[0];
  double rc        = params[1];
  double gamma     = params[2];
  double q_disk    = params[3];

  double pix_scale = IMSIZE/npix;

  for(i = 0; i < npix*npix*NZ;i++){
    double dz = i < npix*npix*(NZ-1) ? (z[i] - z[i+ npix*npix]) : 0 ; /*stupid kludge for now */ 
    model_x = x[i]*combined_transform[0][0] +
      y[i]*combined_transform[0][1] +
      z[i]*combined_transform[0][2];
    
    model_y = x[i]*combined_transform[1][0] +
      y[i]*combined_transform[1][1] +
      z[i]*combined_transform[1][2];
    
    model_z = x[i]*combined_transform[2][0] +
      y[i]*combined_transform[2][1] +
      z[i]*combined_transform[2][2];
  
    /* convert to r,theta,z */ 
    model_r = sqrt(model_x*model_x + model_y*model_y);
    model_phi= atan2(model_y,model_x) ;
    model_phi = model_phi > 0 ? model_phi : model_phi + PI;
    
    if((model_r > 0 )){ // 2*SUBPIXREGION*AU)){
      temp[i]      = T10*pow(model_r/(10*AU),-1*q_disk);
      h         =  sqrt(KB*temp[i]*AD_GAMMA/MP) * (1/sqrt((GG * MSTAR/(pow(model_r,3)))));
      alpha     = 0.1 *sigma_0/(sqrt(2*PI)*h)*exp(-model_z*model_z/(2*h*h))*pow(model_r / rc, -1*gamma)*exp(-1.0*pow(model_r/rc,2-gamma));
    }

    else {
      double avg_temp = 0;
      double avg_h    = 0;
      double avg_alpha= 0;
      double local_temp;
      double local_h;
      double local_alpha;
      double subpix = subpixels*subpixels;
      
      for(k = -(subpixels-1)/2; k <= (subpixels -1)/2; k++){
	for(l = -(subpixels-1)/2; l <= (subpixels -1)/2; l++){
	  model_x = (x[i] + k*pix_scale*DIST/subpixels)*combined_transform[0][0] +
	    (y[i] + l*pix_scale*DIST/subpixels)*combined_transform[0][1] +
	    z[i]*combined_transform[0][2];
	  
	  model_y = (x[i] + k*pix_scale*DIST/subpixels)*combined_transform[1][0] +
	    (y[i] + l*pix_scale*DIST/subpixels)*combined_transform[1][1] +
	    z[i]*combined_transform[1][2];
	  
	  model_z = (x[i] + k*pix_scale*DIST/subpixels)*combined_transform[2][0] +
	    (y[i] + l*pix_scale*DIST/subpixels)*combined_transform[2][1] +
	    z[i]*combined_transform[2][2];
	  
	  model_r = sqrt(model_x*model_x + model_y*model_y);	  
	  
	  if(model_r != 0) {
	    local_temp = T10*pow(model_r/(10*AU),-1*q_disk);	  
	    avg_temp += local_temp;
	  }
	  
	  if(model_r > SUBPIXREGION*AU) 
	    continue;
	  if(model_r == 0){
	    subpix -=1;
	    continue;
	  }
	  local_h    = sqrt(KB*local_temp*AD_GAMMA/MP) * (1/sqrt((GG * MSTAR/(pow(model_r,3)))));
	  local_alpha = 0.1 *sigma_0/(sqrt(2*PI)*local_h)*exp(-model_z*model_z/(2*local_h*local_h))*
	    pow(model_r / rc, -1*gamma)*exp(-1.0*pow(model_r/rc,2-gamma)); 
	  
	  
	  if((model_r != 0) &&  (local_temp < DUST_SUB_T))
	    avg_alpha +=  local_alpha;
	  
	}
      }
      
      temp[i] = (subpix != 0) ? avg_temp/subpix : 10;
      alpha = (subpix != 0) ? avg_alpha/subpix : 0;
    }
    
    d_tau  = alpha*dz;
    if((i % npix == (int)(npix/2))  && ((i /(int)(npix) % npix) == (int)(npix/2)))
      {} //      fprintf(stderr,"dtau = %g\n",d_tau); 
    if(i / (npix*npix) == 0)
      tau[i] = d_tau;
    else
      tau[i] = tau[i - npix*npix] + d_tau;    
  }
  return;
}


double total(double *arr, int npix){
  int i= 0;
  double total = 0;
  for(i=0;i<npix*npix;i++)
    total = total+arr[i];
  return total;
}


void integrate_radiative_transfer(double *im, double *proj_tau, double *temp, double *tau, int npix){
  int i;
  double pix_scale = IMSIZE/npix;
  double flux;
  for(i = npix*npix*NZ-1; i >= npix*npix;i--){ 
    double source_fn = temp[i] > 0 ? 2*PLANCK*NU*NU*NU/(CC*CC) * 
      1.0/(exp(PLANCK*NU/(KB*temp[i]))-1) : 0 ;
    
      im[i % (npix*npix)] =im[i % (npix*npix)] 
	+ source_fn*exp(-1*tau[i])*(tau[i] - tau[i-npix*npix]);
      
      if(i > npix*npix*NZ - npix*npix){
	proj_tau[npix*npix*NZ-i] = tau[i];
      }
  }
  

  for(i = 0; i < npix*npix;i++)
    im[i] = im[i] * pow(pix_scale*DTOR/3600.0,2)*TOJY;

  return;
}




int main(int argc, char *argv[]){
  
  if(argc != 5){
    fprintf(stderr, "Usage %s: [npix] [imfilename] [taufilename] [subpix]\n\nThis is a test code to produce 2D images of disks and envelopes in YSO systems.");
    exit(-1);
  }

  /* image / vis size declarations */

  int npix    = atoi(argv[1]);
  int size[2] = {(int)npix,(int)npix};
  int dim     = 2;
  int cube[3] = {(int) NZ,(int)npix,(int)npix};
  int cubed   = 3;
  double pix_scale = IMSIZE/npix;


  struct model_fft   *modelvis;
  struct image       *modelim;

  /* subgrid for high gradient regions */

  double subpixels = atof(argv[4]); 

  /* model parameterization */
  double disk_param[6] = {2E1,15*AU,1,0.5,0.0,0.0}; 
  double env_param[3]  = {1E-18,0.5,1.0};
  double inc = disk_param[4];
  double pa  = disk_param[5];

  /* arrays to hold grids */

  double *xi,*yi,*zi;
  double *temp;       /* temperature */
  double *tau;        /* optical depth CUBE -- optical depth at given slice of cube. */
  double *tauplane;   /* the first plane ot the optical depth cube. tau through source */


  /* counters */
  int i;
  
  /* coordinate transformation */

  double **combined_transform;

  /* allocate the structs that hold the image and the FFT of the image */

  modelvis = (struct model_fft*) malloc(sizeof(struct model_fft));
  modelvis->real = (double *)malloc(sizeof(double)*npix*npix);
  modelvis->imag = (double *)malloc(sizeof(double)*npix*npix);
  modelvis->u = (double *)malloc(sizeof(double)*npix*npix);
  modelvis->v = (double *)malloc(sizeof(double)*npix*npix);
  
  modelim  = (struct image*) malloc(sizeof(struct image));
  modelim->im = (double *)malloc(sizeof(double)*npix*npix);
  modelim->npix = (int)(npix);
  modelim->pix_scale = pix_scale;

  /*setup image grid*/

  xi = (double *)malloc(sizeof(double)*npix*npix*NZ);
  yi = (double *)malloc(sizeof(double)*npix*npix*NZ);
  zi = (double *)malloc(sizeof(double)*npix*npix*NZ);
  
  temp = (double*) malloc(sizeof(double)*npix*npix*NZ);
  tau = (double*) malloc(sizeof(double)*npix*npix*NZ);
  tauplane = (double*)malloc(sizeof(double)*npix*npix);

  combined_transform =(double **)malloc(sizeof(double*)*3);

  /* check to make sure that we have enough memory */

  if(modelvis == NULL || modelim == NULL || modelvis->real==NULL || modelvis->imag==NULL ||
     modelvis->u == NULL || modelvis->v == NULL || modelim->im == NULL || combined_transform == NULL 
      || xi == NULL || yi == NULL || zi==NULL || 
     temp==NULL || tau==NULL ||  tauplane==NULL){
    fprintf(stderr, "Memory allocation failure. Run on a better machine, dummy.\n");
    exit(-1);
  }

  /* initialize and define geometry */
  
  if(DEBUG)
    fprintf(stderr,"Defining coordinate system...");
  define_coordinate_transform(combined_transform,inc,pa);
  create_grid(xi,yi,zi,npix,inc,pa);
  if(DEBUG)
    fprintf(stderr,"done.\n");
  initialize_image(modelim->im, npix);
  if(DEBUG)
    fprintf(stderr,"Computing temperature and optical depth from model...");
  compute_temp_tau(xi, yi, zi, combined_transform, disk_param, temp, tau, npix,subpixels);  
  if(DEBUG)
    fprintf(stderr,"done.\n");
  
  /* free memory we don't need anymore */
  free(xi); free(yi); free(zi);
  
  if(DEBUG)
      fprintf(stderr,"Integrating radiative transfer eqn...");
  integrate_radiative_transfer(modelim->im,tauplane,temp,tau,npix);

  /* free memory we don't need anymore */
  free(temp);free(tau);//;free(tauplane);
  if(DEBUG)
    fprintf(stderr,"done.\n");

  /* deposit images and shiz into fits files */
  initialize_python();
  if(DEBUG)
    fprintf(stderr,"Exporting images to fits...."); 
  
  make_fits_model_im("before.fits",modelim->im,size,dim,pix_scale);   
  make_fits_model_im("optical_depth.fits",tauplane,size,dim,pix_scale);   
  resort_image(modelim->im, npix);
  make_fits_model_im("after.fits",modelim->im,size,dim,pix_scale);   
  if(DEBUG)
    fprintf(stderr,"done.\n");  
  /* now that we've got the FFT ready, let's do it */
  if(DEBUG)
    fprintf(stderr,"Computing FT...");
  
  compute_fft(modelim,size,dim,modelvis);
  make_fits_model_im("padded_real.fits",modelvis->real,size,dim,pix_scale);   
  make_fits_model_im("padded_imag.fits",modelvis->imag,size,dim,pix_scale);   
  

  if(DEBUG)
    fprintf(stderr,"done.\n");
  if(DEBUG)
    fprintf(stderr,"Testing interpolating FT...");

  free(modelvis->real);free(modelvis->imag);free(modelvis->u);free(modelvis->v);
  free(modelvis);

  {
    int pad = 4;
    struct image *padded_image;
    int sizepad[2]={(int)(npix*pad),(int)(npix*pad)};    
    struct model_fft *modelvis;

    padded_image = (struct image*)malloc(sizeof(struct image));
    padded_image->im=  (double *)malloc(sizeof(double)*npix*npix*pad*pad);
    padded_image->pix_scale = pix_scale;
    modelvis = (struct model_fft*) malloc(sizeof(struct model_fft));
    modelvis->real = (double *)malloc(sizeof(double)*npix*npix*pad*pad);
    modelvis->imag = (double *)malloc(sizeof(double)*npix*npix*pad*pad);
    modelvis->u = (double *)malloc(sizeof(double)*npix*npix*pad*pad);
    modelvis->v = (double *)malloc(sizeof(double)*npix*npix*pad*pad);

    pad_image(modelim->im,padded_image->im,npix,pad);
    free(modelim->im);    free(modelim);    
    make_fits_model_im("after.fits",modelim->im,size,dim,pix_scale);   
    compute_fft(padded_image,sizepad,dim,modelvis);

    dump_vis(modelvis,npix*npix*pad*pad);

    make_fits_model_im("padded_real.fits",modelvis->real,sizepad,dim,pix_scale);   
    make_fits_model_im("padded_imag.fits",modelvis->imag,sizepad,dim,pix_scale);   
    free(modelvis->real);free(modelvis->imag);free(modelvis->u);free(modelvis->v);
    free(modelvis);
    free(padded_image->im);free(padded_image);    
  }
  if(DEBUG)
    fprintf(stderr,"done!\n");
  //quit_python(); /* actually causes a seg fault on my machine if this is uncommented...*/

   return 0;
}




