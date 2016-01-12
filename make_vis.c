#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "fftw3.h"
#include "make_vis.h"
#include "constants.h"
/* essentially just a function to compute an fft */

void compute_fft(struct image *im, int *dims, int dim, struct model_fft *fft){

  fftw_plan plan;
  fftw_complex *in_im, *out_vis;
  
  int i;
  int npixels = 1;
  double du = 1.0/(im->pix_scale * dims[0] * DTOR * 1.0/3600.);  

  for(i = 0; i < dim; i++)
    npixels = npixels * dims[i];
  
  //Allocate arrays.
  in_im = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npixels);
  out_vis = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * npixels);
  plan = fftw_plan_dft_2d(dims[0], dims[1], in_im, out_vis, FFTW_FORWARD, FFTW_MEASURE);
  //Fill in arrays with the pixelcolors.
  fprintf(stderr,"Number of pixels = %d\n",npixels);
  for(i=0;i<npixels;i++){
    in_im[i][0] = im->im[i];
    in_im[i][1] = 0;
  }

  //Forward plans.
  //  fprintf(stderr,"%d\t%d\n",dims[0],dims[1]);

  //Forward FFT.
  fftw_execute(plan);
  for(i=0;i<npixels;i++){
    //    if(out_vis[i][0] != 0)
    //fprintf(stderr,"I'm not 0!");
    fft->real[i]=out_vis[i][0]; 
    fft->imag[i]=out_vis[i][1];
    /* maybe try to compute u and v */
    fft->u[i]  = (i % dims[0] <= dims[0]/2) ? (i % dims[0])  * du : (i % dims[0] - dims[0])*du;
    fft->v[i]  = (i / dims[0] <= dims[1]/2) ? (i / dims[0])  * du : (i / dims[0] - dims[1])*du;  
  }  
  
  fftw_free(in_im);
  fftw_free(out_vis);
   
  return;

}

