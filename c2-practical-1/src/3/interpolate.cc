#include <config.h>
#include <cassert>

#include "../common/grid.hh"
#include "../common/vtkout.hh"
#include "../common/vector.hh"

#include "../problem/problem3.hh"

using namespace std;

//////////////////////////////////////////////////////////////
//
// calculate the linear interpolation of a given function
//
//////////////////////////////////////////////////////////////
void interpolate(const GridType& grid,
                 const Problem& problem,
                 Vector<double> &uh)
{
  //////////////////////////////////////////////////////////////////////////
  //
  //  Assignment 1:
  //
  //  Compute the linear lagrange interpolation of problem.u(x)
  //  1) iterate over the elements K of the grid
  //  2)   on one element K iterate over the verticies p
  //  3)     store the value u(p) in the correct position of the vector uh
  //         using the correct local to global mapping
  //
  //////////////////////////////////////////////////////////////////////////

  const ElementIteratorType end = grid.end();
  for( ElementIteratorType it = grid.begin(); it != end; ++it )
  {
    const ElementType& element = *it;
    for(int j = 0; j < element.nCorners ; ++j)
    {
      //inedx(i,j) = global number of object j, which has co-dimesion i
      uh[element.index(2,j)] = problem.u(element[j]);
    }
  }

  //EXERCISE("Compute linear lagrange interpolation");
}
//////////////////////////////////////////////////////////////
//
// calculate the error between the exact solution and a
// linear discrete function
//
//////////////////////////////////////////////////////////////
void error(const GridType& grid,
           const Problem& problem,
           const Vector<double> &uh,
           double &l2error,
           double &h1error)

{ // compute H1 error
  l2error=0.;
  h1error=0.;

  // loop over all elements
  const ElementIteratorType endit = grid.end();
  for(ElementIteratorType it = grid.begin();
      it != endit; ++it)
  {
    // get element
    const ElementType& element = *it;

    double locl2error=0;
    double loch1error=0.;

    // coordinate of bary center
    LocalCoordType bary(1./3.);
    GlobalCoordType x = element.global( bary );

    //////////////////////////////////////////////////////////////////////////
    //
    //  Assignment 1:
    //
    //  Store L2 error in locl2error:
    //  1) evaluate u at the bary center x
    //  2) evaluate uh at the bary center (using coordinate in the reference
    //     element is easiest)
    //  3) compute the approximation
    //
    //////////////////////////////////////////////////////////////////////////


    //EXERCISE("Store int_K(uh-u)^2(x) in locl2error");

    double ux = problem.u(x);
    double uhx = problem.u(element[0]);
    uhx += problem.u(element[1]);
    uhx += problem.u(element[2]);
    uhx /= 3.0;

    locl2error = element.area()*(ux - uhx)*(ux - uhx);

    //////////////////////////////////////////////////////////////////////////
    //
    //  Assignment 1:
    //
    //  Store |duh(x)-du(x)|^2 in loch1error:
    //  1) evaluate du at the bary center x
    //  2) compute the gradient of uh (is constant in this case)
    //     i)  compute the gradient in local coordinates
    //     ii) compute the gradient in global coordinates by using the chain rule
    //  3) compute the approximation
    //
    /////////////////////////////////////////////////////////////////////////a
  //  EXERCISE("Store int_K|duh(x)-du(x)|^2 in loch1error");


  GlobalCoordType du(0);
  problem.du(x,du);

  //vector of coefficients
  vector<double> alpha(3);
  for(int i = 0 ; i < alpha.size() ; ++i)
  {
    alpha[i] = uh[ element.index(2,i) ];
  }

  //gradient of the reference triangle basis
  LocalCoordType x_(0);
  LocalCoordType y_(0);
  LocalCoordType z_(0);

  x_[0] = -1;
  x_[1] = -1;
  y_[0] = 1;
  y_[1] = 0;
  z_[0] = 0;
  z_[1] = 1;

  //gradient of the local triangle basis functions
  GlobalCoordType gradx_(0);
  GlobalCoordType grady_(0);
  GlobalCoordType gradz_(0);

  element.gradientGlobal(x_,gradx_);
  element.gradientGlobal(y_,grady_);
  element.gradientGlobal(z_,gradz_);

  //computing grad u_h = \sum_{i=1}^3 \alpha_i grad \phi_K^i
  GlobalCoordType duh(0);
  for(int i = 0; i < 2 ; ++ i)
  {
    gradx_[i] *= alpha[0];
  }
  for(int i = 0; i < 2 ; ++ i)
  {
    grady_[i] *= alpha[1];
  }
  for(int i = 0; i < 2 ; ++ i)
  {
    gradz_[i] *= alpha[2];
  }

  duh += gradx_ + grady_ + gradz_;

    double square = 0.0;
    for(int i = 0; i < 2; ++i)
      {
        square += (duh[i] - du[i])*(duh[i] - du[i]);
      }
    //H1 norm error!
    loch1error = element.area()*square;

    //add local error to global error
    l2error += locl2error;
    h1error += loch1error;
  }
  l2error = sqrt(l2error);
  h1error = sqrt(h1error);
}

void compute(GridType& grid,
             const Problem& problem,
             const int refines)
{
  const int dimension = GridType::dimension;
  double l2error[refines];
  double h1error[refines];
  double eocl2[refines];
  double eoch1[refines];

  double hold=1./sqrt(grid.size(dimension));


  cout << "uh" << "\t  "
       << "hnew" << " "
       << " \t\t "
       << "l2error[n]" << " \t " << "eocl2[n]"
       << " \t\t "
       << "h1error[n]" << " \t " << "eoch1[n]"
       << " \t\t "
       << "inequality const"
       << endl;
  // start computation
  for (int n=0; n < refines; ++n)
  {
    double hnew=1./sqrt(grid.size(dimension));

    // approximate solution
    Vector<double> uh( grid.size(dimension) );

    // compute approximate solution
    interpolate(grid,problem,uh);

    // compute the approximation error
    error( grid, problem, uh, l2error[n], h1error[n] );

    if (n>0)
    {
      const double hfraction=hold/hnew;
      eocl2[n]=(log(l2error[n-1])-log(l2error[n]))/
                log(hfraction);

      eoch1[n]=(log(h1error[n-1])-log(h1error[n]))/
                log(hfraction);
    }
    else
    {
      eocl2[0]=-1;
      eoch1[0]=-1;
    }

    cout << uh.size() << "\t  "
         << hnew
         << " \t\t "
         << l2error[n] << " \t " << eocl2[n]
         << " \t\t "
         << h1error[n] << " \t " << eoch1[n]
         << " \t\t "
         << (l2error[n] + hnew*h1error[n])/(hnew*hnew)
         << endl;

    // new h is old h now
    hold = hnew;

    // optput data
    Output out( grid, n );
    out.addVertexData( uh, "I_h(u)" );
    out.write();

    // refine all elements of grid twice to half grid width
    grid.globalRefine( GridType :: dimension );

  } // end of loop for (n=..)
}

/**************************************************************/
int main(int argc, char ** argv)
{
  // read in command line parameters
  string gridname ( "../problem/cube.dgf" );
  int loops = 2;
  int prob = 0;

  if (argc<4)
  {
    cout << "Usage: "<< argv[0] << " <grid file> <loops> <problem> " << endl;
    cout << "                         loops > 0" << endl;
    cout << "                         problem = 0 (Sinus)" << endl;
  }
  else
  {
    gridname = argv[1];
    loops  = atoi( argv[2] );
    prob     = atoi( argv[3] );
  }

  Problem *problem=0;
  switch (prob)
  {
    case 0: problem = new ProblemSin();
            break;
    case 1: problem = new Problem1();
            break;
    case 2: problem = new Problem2();
            break;
    case 3: problem = new Problem3();
            break;
    case 4: problem = new Problem4();
            break;
    case 5: problem = new Problem5();
            break;
    default: std::cerr << "Wrong problem number " << prob
                       << ". Should be between less than 1." << std::endl;
             return 1;
  }

  // initialize grid
  GridType grid ( gridname.c_str() ); // read the macro grid file
  // make initial grid
  grid.globalRefine( 4 );

  compute( grid, *problem, loops );

  delete problem;

  return 0;
}
