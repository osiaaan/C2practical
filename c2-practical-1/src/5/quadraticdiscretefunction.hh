#ifndef QUADRATICDISCRETEFUNCTION_HH_INCLUDED
#define QUADRATICDISCRETEFUNCTION_HH_INCLUDED

//- include vector implementation 
#include "../common/vtkout.hh"
#include "../common/vector.hh"
#include "../common/discretefunction.hh"

//////////////////////////////////////////////////////////////////////////
//
//  Exercise 5: 
//
//  Implement quadratic shape function following the interface defined by
//  the implementation of the linear shape functions.
//
//////////////////////////////////////////////////////////////////////////


/** \brief QuadraticBaseFunction implements the
    \ref BaseFunctionSet for a piecewise quadratic Lagrange space.
*/
class QuadraticBaseFunction {
public:
  enum {locNrDof = 6};

  /** \brief constructor creating a QuadraticBaseFunction set
      \param[in] grid  grid that the base function set is
       based on
  */
  QuadraticBaseFunction(const GridType& grid) :
    _edgeoffset(grid.size(2)),
    _size(grid.size(2)+grid.size(1))
    // _edgeoffset(grid.nVertices()),
    // _size(grid.nVertices()+grid.nEdges())
   {
     _point[0][0]=0.;
     _point[0][1]=0.;
     _point[1][0]=1.;
     _point[1][1]=0.;
     _point[2][0]=0.;
     _point[2][1]=1.;

     _point[3][0]=0.5;
     _point[3][1]=0.;
     _point[4][0]=0.;
     _point[4][1]=0.5;
     _point[5][0]=0.5;
     _point[5][1]=0.5;
  }
  
  /** \brief global number of degrees of freedom (DoF)
      \return number of DoFs
  */
  int size() const
  {
    return _size;
  }

  /** \brief local to global DoF number mapping
      \param[in]  element Element the contains local DoFs
      \param[in]  vx local DoF number
 
      \return global DoF number
  */
  int map(const ElementType& element,
          const int vx) const 
  {
    EXERCISE("Implement the correct local to global mapping here");
    return -1;
  }

  /** \brief returns true if local dof is located at
       boundary
      \param edge edge that is checked
      \param vx local dof number 

      \return \b true if dof is on edge
       boundary, \b false otherwise
 */
  bool onBnd(const int edge, const int vx) const 
  {
    if ( edge == 0 )
      return (vx < 2 || vx-3==0);
    else
      return ( (vx < 3 && (vx == 2 || vx-1 == edge)) || vx-3 == edge);
  }

  /** \brief evaluate base function
      \param[in]  i   number of base function to
      evaluate
      \param[in]  lambda  local coordinate
       where base function is evaluated

      \return evaluation of the i-th base
       function
  */
  double evaluate(const int i,
                  const LocalCoordType& lambda) const
  {
    EXERCISE("Return the shape function value here");
    return 0.;
  }

  /** \brief evaluate gradient of base function on the reference element 
      \param[in]  i   number of base function to
       evaluate
      \param[in]  lambda  local coordinate
       where the gradient of the
       base function is evaluated

     \return evaluation of the i-th base function gradient in local  coordinates 
  */
  LocalCoordType gradientLocal(const int i,
                               const LocalCoordType& lambda) const
  {
    EXERCISE("Return the shape function gradient here");
    return LocalCoordType(0);
  }
  
  /** \brief return lagrange point 
      \param[in] i number of lagrange point 

      \return const reference to lagrange point 
  */
  const LocalCoordType& point( const int i) const 
  {
    return _point[i];
  }

protected:
  LocalCoordType _point[locNrDof];
  LocalCoordType _gradphi[locNrDof];
  int _edgeoffset;
  int _size;
};

#endif
