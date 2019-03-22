/* A program for 1d interpolation
 * Compile with:
 * gcc -o exercise1a exercise1a.c -g -Wall -lm
 */


#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979
#endif

/* function to interpolate */
double u(double x)
{
  return sin(2.*M_PI*x); /* a smooth function */
#if 0 /* a possible way to remove a few lines of code with comments included */
  return sqrt(fabs(x));  /* a not quite so smooth function */
#endif
#if 0
  if (fabs(x)>1e-10)
    return x*sin(1/x);   /* a realy unpleasant function */
  else
    return 0.;
#endif
}

/* the program always starts with this function */
int main(int argc, char ** argv) 
{
  /* domain definition */
  double a = -0.5;
  double b =  1.5;
  /* discretization */
  const int N = 150;
  /* array for coordinates of grid points */
  double x[N];
  /* array for discrete values */
  double uh[N];
  /* compute grid width */
  double h = (double)(b-a) / (double)(N-1);
  /* for output */
  double error = 0.;
  FILE *out;

  /* loop variable */
  int i;

  /* store value at grid points */
  for (i=0;i<N;++i)
  { 
   //Hello world 
    x[i] = a+h*i;
    uh[i] = u(x[i]);
  }

  /* output function and interpolation at cell centers into file 
     also compute interpolation error */
  out = fopen("interpol1a.gnu", "w");
  if (out == NULL) 
  {
    printf("Can't open input file interpol.gnu!\n");
    exit(1);
  }

  for (i=0;i<N-1;++i)
  {
    double xplushalf = x[i]+0.5*h;
    double ux = u(xplushalf);
    double uhx = 0.5*(uh[i]+uh[i+1]);
    fprintf(out,"%f %f %f\n",xplushalf,ux,uhx);
    error += h*fabs(ux-uhx);
  }
  printf("error: %f\n",error);

  fclose(out);
  return 0;
}

