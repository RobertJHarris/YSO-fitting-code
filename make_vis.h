#ifndef _MAKE_VIS_H
#define _MAKE_VIS_H 

struct image{
  double *im;
  double npix;
  double pix_scale;
};


struct model_fft {
  double *real;
  double *imag;
  double *u;
  double *v;
};


void compute_fft(struct image *im, int *dims, int dim, struct model_fft *fft);
#endif
