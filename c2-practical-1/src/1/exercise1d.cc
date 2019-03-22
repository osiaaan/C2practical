/* A program for 1d interpolation
 * Compile with:
 * make exercise1d
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

// forward declaration:
// Element needs to make Grid friend and Grid stores Element - 
// both need to know each other...
class Grid;

struct Element
{
  Element() : h_(-1) {}
  double area() const
  {
    return h_;
  }
  double xM() const
  {
    return xM_;
  }
  double vertex(int i) const
  {
    assert( i==0 || i==1 );
    if (i==0)
      return xM_-0.5*h_;
    else
      return xM_+0.5*h_;
  }
  int index() const
  {
    return index_;
  }
  int index(int i) const
  {
    assert( i==0 || i==1 );
    return (i==0)? index_:index_+1; 
  }
  private:
  friend class Grid; // Grid needs to access the setup method
  void setup(int index,double h,double xM)
  {
    index_ = index; 
    h_ = h;
    xM_ = xM;
  }
  int index_;
  double h_;
  double xM_;
};

class Grid
{
  double a_;
  double b_;
  double N_;
  double h_;
  typedef vector<Element> ElementStorage;
  ElementStorage element_;
  public:
  typedef ElementStorage::const_iterator ElementIterator;
  Grid(double a, double b, double N)
  : a_(a), b_(b), N_(N),
    element_(N-1)
  {
    // compute grid width 
    h_ = (double)(b-a) / (double)(N-1);
    // store coordinates of grid points
    for (unsigned int i=0;i<element_.size();++i)
    {
      element_[i].setup(i,h_,a_+h_*(i+0.5));
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
  const ElementIterator begin() const
  {
    return element_.begin();
  }
  const ElementIterator end() const
  {
    return element_.end();
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
  vector<double> uh(grid.size(1));
  // store value at grid points
  Grid::ElementIterator endIt = grid.end();
  for ( Grid::ElementIterator it = grid.begin(); it != endIt; ++it )
  {
    const Element &element = *it;
    uh[element.index(0)] = (*problem).u(element.vertex(0));
    uh[element.index(1)] = (*problem).u(element.vertex(1));
  }

  /* output function and interpolation at cell centers into file 
     also compute interpolation error */
  double error = 0.;
  ofstream out("interpol1d.gnu");
  if (!out) 
  {
    cout << "Can't open input file interpol.gnu!\n";
    return 1;
  }

  // iterate over x - access to uh is as before for the moment
  for ( Grid::ElementIterator it = grid.begin(); it != endIt; ++it )
  {
    const Element &element = *it;
    double ux = problem->u(element.xM());
    double uhx = 0.5*(uh[element.index(0)]+uh[element.index(1)]);
    out << element.xM() << " " << ux << " " << uhx << endl;
    error += element.area()*abs(ux-uhx);
  }
  cout << "error: " << error << endl;

  delete problem;

  return 0;
}

