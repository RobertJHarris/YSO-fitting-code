#ifndef _FITSIO_H
#define _FITSIO_H 

void initialize_python();
void quit_python();
void init_numpy();
void  *set_path();
void make_fits_model_im(char *filename,double *im, int *size, int dim, double pix_scale);
#endif
