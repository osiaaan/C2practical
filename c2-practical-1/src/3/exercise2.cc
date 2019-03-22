// needed by Dune 
#include <config.h>
#include <cassert>

// include grid 
#include "../common/grid.hh"
#include "../common/vtkout.hh"
#include "../common/vector.hh"

double area( const ElementType& element )
{
  //////////////////////////////////////////////////////////////////////////
  //
  //  Exercise 2: 
  //
  //  Return the area of the element given as parameter to this function.
  //  Use the determinant of the reference mapping to compute the area:
  //    2|T| = |det DF_T|
  //
  //  Use the following method on the Element class:
  //    operator[](int): element[i] returns the coordinate vector
  //                     (element[i][0], element[i][1]) 
  //                     of the ith vertex.
  //
  //////////////////////////////////////////////////////////////////////////
  double area = 0;

  double a = element[1][0] - element[0][0];
  double b = element[2][0] - element[0][0];
  double c = element[1][1] - element[0][1];
  double d = element[2][1] - element[0][1];
  area = 0.5*abs( a*d - b*c );

  // alternative 1: element.area(); 
  // alternative 2: 0.5*element.integrationElement();

  assert( area > 0.0 );
  return area;
}
GlobalCoordType baryCenter( const ElementType& element )
{
  //////////////////////////////////////////////////////////////////////////
  //
  //  Exercise 2: 
  //
  //  Return the baryCenter of the element given as parameter to this function.
  //  Use the average of the coordinates of the three vertices:
  //    omega_T = sum_{i=1}^3 v^T_i
  //
  //  Use the following method on the Element class:
  //    operator[](int): element[i] returns the coordinate vector
  //                     (element[i][0], element[i][1]) 
  //                     of the ith vertex.
  //
  //////////////////////////////////////////////////////////////////////////
  
  // setup return value and initialize with zero
  // GlobalCoordType is a vector so that is has an operator[] to access the
  // entries in the vector; but is also has operator+=, operator/= and so on and so
  // forth...
  GlobalCoordType bary(0);

  const int numberOfVertices = element.nCorners;
  for ( int i=0; i<numberOfVertices; ++i )
    bary += element[i];
  bary /= double(numberOfVertices);

  // alternative: LocalCoordType localBary(1./3.);
  //              return element.global( localBary );
  
  return bary;
}

// In this program we want to compute the area of a domain Omega_h defined by a grid
// and store a function using a piecewise constant representation
int main(int argc, char ** argv) 
{
  // read in command line parameters (giving defaults)
  string gridname ( "../problem/cube.dgf" );
  int startRefines = 8;

  if (argc<3) 
  {
    cout << "Usage: "<< argv[0] << " <grid file> <start refines> " << endl;
    cout << "                         start refines > 0" << endl;
    cout << "Continuing using default values: " 
         << argv[0] << " " << gridname << " " << startRefines << std::endl;
  }
  else 
  {
    gridname = argv[1]; 
    startRefines  = atoi( argv[2] );
  }
  GridType grid( gridname );
  grid.globalRefine( startRefines );

  // vector for storing piecewise constant data
  // grid.size(0) returns number of element (codim zero entites) of the grid
  // ... <double> is a template argument which defines the type of the data to store 
  // ... this is an advanced C++ feature and can be ignored here
  Vector<double> uh( grid.size(0) );

  // variable for storing area of domain
  double domainArea = 0;

  // iterate over all elements of the grid 
  const ElementIteratorType end = grid.end();
  for(ElementIteratorType it = grid.begin(); it != end; ++it) 
  {
    // show element number and area of element 
    const ElementType& element = * it;
    // add area of element to total area 
    domainArea += area( element );
    // store value of function in barycenter of element
    uh[ element.index(0,0) ] = sin( 2.*M_PI*baryCenter( element ).two_norm() );
  }

  // print the area to screen
  std::cout << "Total area: " << domainArea << std::endl;

  // output grid and data to file 
  Output output(grid,1);
  output.addCellData( uh, "function" );
  output.write();

  return 0;
}

