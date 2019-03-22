/* A program for 1d interpolation
 * Compile with:
 * g++ -o exercise1c exercise1c.cc -g -Wall
 */

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <vector>
using namespace std;

/* function to interpolate */
struct Problem
{
  virtual double u(double x) const = 0;
};
struct Smooth : public Problem
{
  double u(double x) const 
  {
    return sin(2.*M_PI*x); // a smooth function 
  }
};
struct C0 : public Problem
{
  double u(double x) const 
  {
    return sqrt(fabs(x));  // a not quite so smooth function
  }
};
struct ReallyBad : public Problem
{
  double u(double x) const 
  {
    if (fabs(x)>1e-10)
      return x*sin(1/x);   // a realy unpleasant function 
    else
      return 0.;
  }
};

class Grid
{
  double a_;
  double b_;
  double N_;
  double h_;
  vector<double> x_;
  public:
  Grid(double a, double b, double N)
  : a_(a), b_(b), N_(N),
    x_(N)
  {
    // compute grid width 
    h_ = (double)(b-a) / (double)(N-1);
    // store coordinates of grid points
    for (int i=0;i<N;++i)
    {
      x_[i] = a_+h_*i;
    }
  }
  double h() const
  {
    return h_;
  }
  double size(int c) const
  {
    assert( c==0 || c==1 );
    return (c==0)?N_-1:N_;
  }
  double x(int i) const
  {
    assert( 0<= i && i<N_ );
    return x_[i];
  }
  double baryCenter(int i) const
  {
    assert( 0<= i && i<N_-1 );
    return x(i)+0.5*h();
  }
};

/* the program always starts with this function */
int main(int argc, char ** argv) 
{
  // domain definition 
  double a = -0.5;
  double b =  1.5;
  // read in command line parameters
  if (argc<3) 
  {
    cout << "Usage: "<< argv[0] << " <problem> <N>" << endl;
    return 1;
  }
  int prob = atoi( argv[1] );
  // discretization (does not need to be const)
  int N = atoi( argv[2] );

  Problem *problem=0;
  switch (prob)
  {
    case 0: problem = new Smooth();
            break;
    case 1: problem = new C0();
            break;
    case 2: problem = new ReallyBad();
            break;
    default: std::cerr << "Wrong problem number " << prob
                       << ". Should be less than 3." << std::endl;
             return 1;
  }
  // construct grid
  Grid grid(a,b,N);
  // array for discrete values
  vector<double> uh(N);

  // store value at grid points
  for (int i=0;i<grid.size(1);++i)
  {
    uh[i] = (*problem).u(grid.x(i));
  }

  /* output function and interpolation at cell centers into file 
     also compute interpolation error */
  double error = 0.;
  ofstream out("interpol1c.gnu");
  if (!out) 
  {
    cout << "Can't open input file interpol.gnu!\n";
    return 1;
  }

  // iterate over x - access to uh is as before for the moment
  for (int i=0;i<grid.size(0);++i)
  {
    double ux = problem->u(grid.baryCenter(i));
    double uhx = 0.5*(uh[i]+uh[i+1]);
    out << grid.baryCenter(i) << " " << ux << " " << uhx << endl;
    error += grid.h()*abs(ux-uhx);
  }
  cout << "error: " << error << endl;

  delete problem;

  return 0;
}

