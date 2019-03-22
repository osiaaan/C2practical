/* A program for 1d interpolation
 * Compile with:
 * g++ -o exercise1a exercise1a.cc -g -Wall
 */


#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
using namespace std;

/* function to interpolate */
double u(double x)
{
#if 1
  return sin(2.*M_PI*x); // a smooth function 
#endif
#if 0 // a possible way to remove a few lines of code with comments included 
  return sqrt(fabs(x));  // a not quite so smooth function
#endif
#if 0
  if (fabs(x)>1e-10)
    return x*sin(1/x);   // a realy unpleasant function 
  else
    return 0.;
#endif
}

/* the program always starts with this function */
int main(int argc, char ** argv) 
{
  // domain definition 
  const double a = -0.5;
  const double b =  1.5;
  // discretization 
  const int N = 150;
  // array for coordinates of grid points
  double x[N];
  // array for discrete values
  double uh[N];
  // compute grid width 
  double h = (double)(b-a) / (double)(N-1);

  // store value at grid points
  for (int i=0;i<N;++i)
  {
    x[i] = a+h*i;
    uh[i] = u(x[i]);
  }

  /* output function and interpolation at cell centers into file 
     also compute interpolation error */
  double error = 0.;
  ofstream out("interpol1a.gnu");
  if (!out) 
  {
    cout << "Can't open input file interpol.gnu!\n";
    return 1;
  }

  for (int i=0;i<N-1;++i)
  {
    double xplushalf = x[i]+0.5*h;
    double ux = u(xplushalf);
    double uhx = 0.5*(uh[i]+uh[i+1]);
    out << xplushalf << " " << ux << " " << uhx << endl;
    error += h*fabs(ux-uhx);
  }
  cout << "error: " << error << endl;

  return 0;
}

