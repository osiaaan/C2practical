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
  //  Use only the following method on the Element class:
  //    operator[](int): element[i] returns the coordinate vector
  //                     (element[i][0], element[i][1])
  //                     of the ith vertex.
  //
  //////////////////////////////////////////////////////////////////////////
  double area = 0;
  double determinant = (element[1][0] - element[0][0])*(element[2][1] - element[1][1]) -
                       (element[2][0] - element[1][0])*(element[1][1] - element[0][1]);
//EXERCISE("Area computation");

  area = 0.5 * determinant;
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
  //  Use only the following method on the Element class:
  //    operator[](int): element[i] returns the coordinate vector
  //                     (element[i][0], element[i][1])
  //                     of the ith vertex.
  //
  //////////////////////////////////////////////////////////////////////////

  // setup return value and initialize with zero
  // GlobalCoordType is a vector so that is has an operator[] to access the
  // entries in the vector; but is also has operator+=, operator/= and so on and so
  // forth...



  double x = (1.0/3.0)*(element[0][0]+element[1][0]+element[2][0]);
  double y = (1.0/3.0)*(element[0][1]+element[1][1]+element[2][1]);

  GlobalCoordType bary(0);

  bary[0] = x;
  bary[1] = y;

  //EXERCISE("Compute barycenter");

  return bary;
}


GlobalCoordType midpoint( const ElementType& element, int i )
{
  GlobalCoordType mid(0);

  double x = (1.0/2.0)*(elemnt[i][0] + element[(i+1)%3][0]);
  double y = (1.0/2.0)*(elemnt[i][1] + element[(i+1)%3][1]);

  mid[0] = x;
  mid[1] = y;

  return mid;
}


// In this program we want to compute the area of a domain Omega_h defined by a grid
// and store a function using a piecewise constant representation
int main(int argc, char ** argv)
{
  // create grid; the parameter is the name of the grid file
  // you might want to make the grid file a command line parameter (this is an instance
  // of the Grid class)
  GridType grid( "../problem/cube.dgf" );
  // the method globalRefine on the grid refines the grid globally a given number of times
  // you might want to make the number of refinement steps a command line parameter
  grid.globalRefine(12);

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
  // for (const auto &element : grid)
  {
    // the iterator ''points'' to the grid element, to get the element we need to write...
    const ElementType& element = * it;
    // the element is an instance of the Element class which gives access to the
    // coordinates of the element corners and the local to global numbering for all
    // subentities, e.g., element.index() gives the global number of this element,
    // element.vertexIndex(i) the global number of the ith vertex....

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

  //////////////////////////////////////////////////////////////////////////
  //
  //  Exercise 2*:
  //
  //  Compute the error of the approximation and compute the EOC.
  //  You need a higher order quadrature rule to compute the error, e.g.,
  //  e_h = sum_T |T|/3 sum_{i=1}^3 |u(m^T_i)-u_T|
  //  is an approximation of ||u-u_h||_1, where
  //  u_T is the value of the piecewise constant approximation u_h on the
  //  element T and m^T_i is the midpoint of the ith edge of T
  //
  //  To solve this problem you still we only need the
  //    operator[](int)
  //  method on the Element class and the globalRefine method on the Grid class.
  //
  //////////////////////////////////////////////////////////////////////////

  double eh = 0.0;


  vector<double> erro_vec;

  for(ElementIteratorType it = grid.begin() ; it != grid.end(); ++it )
  {
    element.index();
    for(int i = 0 ; i < 3; ++i)
    {
      error_vec.push_back(sin( 2.*M_PI*))
    }
  }
  return 0;
}
