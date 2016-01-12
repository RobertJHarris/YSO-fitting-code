#include <stdlib.h>
#include <math.h>
#include "fftw3.h"

#define N 16
int main (void)
{
  fftw_complex in[N], out[N], in2[N2], out2[N2];        /* double [2] */
  fftw_plan p, q;
  int i;
  int half;
 
  half=(N/2+1);
  /* prepare a cosine wave */
  for (i = 0; i < N; i++)
    {
      //in[i][0] = cos ( 3 * 2 * M_PI * i / N);
      in[i][0] = (i > 3 && i< 12 )?1:0;
      in[i][1] = (i > 3 && i< 12 )?1:0;
      //in[i][1] = 0;
    }

  /* forward Fourier transform, save the result in 'out' */
  p = fftw_plan_dft_1d (N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
  fftw_execute (p);
  for (i = 0; i < N; i++)
    printf ("input: %3d %+9.5f %+9.5f I   %+9.5f %+9.5f I\n", i, in[i][0], in[i][1],out[i][0],out[i][1]);
  fftw_destroy_plan (p);

  for (i = 0; i<N; i++) {out2[i][0]=0.;out2[i][1]=0.;}

  for (i = 0; i<half; i++) {
    out2[i][0]=2.*out[i][0];
    out2[i][1]=2.*out[i][1];
  }
  for (i = half;i<N;i++) {
    out2[N+i][0]=2.*out[i][0];
    out2[N+i][1]=2.*out[i][1];
  }



  /* backward Fourier transform, save the result in 'in2' */
  printf ("\nInverse transform:\n");
  q = fftw_plan_dft_1d (N2, out2, in2, FFTW_BACKWARD, FFTW_ESTIMATE);
  fftw_execute (q);
  /* normalize */
  for (i = 0; i < N2; i++)
    {
      in2[i][0] /= N2;
      in2[i][1] /= N2;
    }
  for (i = 0; i < N2; i++)
    printf ("recover: %3d %+9.1f %+9.1f I\n",
            i, in2[i][0], in2[i][1]);
  fftw_destroy_plan (q);

  fftw_cleanup ();
  return 0;
}
