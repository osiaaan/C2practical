#ifndef CRdiscrete_HH_INCLUDED
#define CRdiscrete_HH_INCLUDED

//- include vector implementation
#include "../common/vtkout.hh"
#include "../common/vector.hh"
#include "../common/discretefunction.hh"


/** \brief CRbasefunction implements the
    \ref BaseFunctionSet for a piecewise quadratic Lagrange space.
*/
class CRBaseFunction {
public:
  enum {locNrDof = 3};

  /** \brief constructor creating a CRbasefunction set
      \param[in] grid  grid that the base function set is
       based on
  */


  CRBaseFunction(const GridType& grid) :
    _size(grid.size(grid.dimension -1))
  {
     _point[0][0]=0.5;
     _point[0][1]=0.;
     _point[1][0]=0.;
     _point[1][1]=0.5;
     _point[2][0]=0.5;
     _point[2][1]=0.5;

     _gradphi[0][0]=0.;
     _gradphi[0][1]=-2.;
     _gradphi[1][0]=-2.;
     _gradphi[1][1]=0.;
     _gradphi[2][0]=2.;
     _gradphi[2][1]=2.;
  }

  /** \brief number of degrees of freedom (DoF)
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
    return element.index( element.dimension - 1, vx );
  }

  /** \brief returns true if local dof is located at
       boundary
      \param edge edge that is checked
      \param vx local vertex number of edge

      \return \b true if dof is on
       boundary, \b false otherwise
 */
 bool onBnd(const int edge, const int vx) const
 {
   return (edge == vx);
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
    double bary[3]={1.-lambda[0]-lambda[1],lambda[0],lambda[1]};
    return 1. - 2. * bary[ 2-i ];
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
    return _gradphi[i];
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
  int _size;
};

#endif