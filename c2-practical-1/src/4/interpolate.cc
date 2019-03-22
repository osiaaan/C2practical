#include <config.h>
#include <cassert>

#include "../common/grid.hh"
#include "../common/vtkout.hh"
#include "../common/vector.hh"

#include "../problem/problem4.hh"

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
  //
  //////////////////////////////////////////////////////////////////////////
  // set all dofs to zero 
  uh.clear();

  // start grid traversal 
  const ElementIteratorType endit = grid.end();
  for(ElementIteratorType it = grid.begin();
      it != endit; ++it) 
  {
    // get reference to element 
    const ElementType& element = *it; 

    // for all vertices do interpolation 
    for(int i=0; i < element.nCorners; ++i) 
    {
      uh[ element.index(2,i) ] = problem.u( element[i] );
    }
  }
}
//////////////////////////////////////////////////////////////    
//
// calculate the error of between the exact solution and a 
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
    //  Store L2 error in locl2error
    //
    //////////////////////////////////////////////////////////////////////////
    double uValue = problem.u( x );

    double uhValue = 0;
    for (int i=0;i<element.nCorners;++i)
    {
      uhValue += uh[ element.index(2,i) ];
    }
    uhValue /= double(element.nCorners);

    // add to local L2 error 
    locl2error = element.area() * (uValue-uhValue) * (uValue-uhValue);

    //////////////////////////////////////////////////////////////////////////
    //
    //  Assignment 1: 
    //
    //  Store |duh(x)-du(x)|^2 in loch1error
    //
    //////////////////////////////////////////////////////////////////////////
    GlobalCoordType duValue;
    problem.du( x, duValue );

    LocalCoordType duhValueLocal;
    duhValueLocal[0] = uh[ element.index(2,1) ] - uh[ element.index(2,0) ];
    duhValueLocal[1] = uh[ element.index(2,2) ] - uh[ element.index(2,0) ];
    GlobalCoordType duhValue;
    element.gradientGlobal( duhValueLocal, duhValue );

    // H1-error
    loch1error = element.area() * (duValue - duhValue).two_norm2();

    // add local error to global error 
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
  double l2error[refines];
  double h1error[refines];
  double eocl2[refines];
  double eoch1[refines];

  double hold=1./sqrt(grid.size(2));

  // start computation
  for (int n=0; n < refines; ++n) 
  {
    double hnew=1./sqrt(grid.size(2));

    // approximate solution
    Vector<double> uh( grid.size(2) );       
  
    // compute approximate solution
    interpolate(grid,problem,uh);
 
    // compute the approximation error
    error( grid, problem, uh, l2error[n], h1error[n] );
    
    //////////////////////////////////////////////////////////////////////////
    //
    //  Assignment 1: 
    //
    //  Compute numerical convergence rates for l2error and h1error
    //
    //////////////////////////////////////////////////////////////////////////
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

    // output eoc table (suitable for using in tex table)
    if (n==0)
      cout << "dofs \t&\t h \t&\t $L^2$ \t& eoc \t&\t $H^1$ \t& eoc \\\\ \n";
    cout << uh.size() << "\t& "
         << hnew
         << " \t&\t " 
         << l2error[n] << " \t& " << eocl2[n]
         << " \t&\t "
         << h1error[n] << " \t& " << eoch1[n]
         << "\\\\"
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
    cout << "                         problem = 0..2" << endl;
  }
  else 
  {
    gridname = argv[1]; 
    loops    = atoi( argv[2] );
    prob     = atoi( argv[3] );
  }

  Problem *problem=0;
  switch (prob)
  {
    case 0: problem = new ProblemSin();
            break;
    case 1: problem = new Problem2019_Ass1_1();
            break;
    case 2: problem = new Problem2019_Ass1_2_3(1./M_PI);
            break;
    case 3: problem = new Problem2019_Ass1_2_3(0.5);
            break;
    case 4: problem = new Problem2019_Ass1_4();
            break;
    case 5: problem = new Problem2019_Ass1_5(sqrt(2.)/2.);
            break;
    case 6: problem = new Problem2019_Ass1_5(0.5);
            break;
    default: std::cerr << "Wrong problem number " << prob
                       << ". Should be less than 7." << std::endl;
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
