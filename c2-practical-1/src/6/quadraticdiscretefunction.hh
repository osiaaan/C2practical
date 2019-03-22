#ifndef QUADRATICDISCRETEFUNCTION_HH_INCLUDED
#define QUADRATICDISCRETEFUNCTION_HH_INCLUDED

//- include vector implementation 
#include "../common/vtkout.hh"
#include "../common/vector.hh"
#include "../common/discretefunction.hh"


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
    _edgeoffset(grid.size(grid.dimension)),
    _size(_edgeoffset+grid.size(grid.dimension-1))
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
    if (vx<3)
      return element.index(element.dimension,vx);
    else
      return _edgeoffset+element.index(element.dimension-1,vx-3);
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
    double bary[3]={1.-lambda[0]-lambda[1],lambda[0],lambda[1]};
    if (i<3)
    {
      return bary[i]*(2.*bary[i]-1.);
    }
    else 
    {
      int j=i-3;
      return 4.*bary[(3-j)%3]*bary[(4-j)%3];
    }
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
    LocalCoordType ret;
    const double bary[3]={1.-lambda[0]-lambda[1],lambda[0],lambda[1]};
    const double dbary[3][2]={{-1.,-1},{1.,0.},{0.,1.}};
    if (i<3) {
      ret[0] =
        dbary[i][0]*(2.*bary[i]-1.)+bary[i]*2.*dbary[i][0];
      ret[1] =
        dbary[i][1]*(2.*bary[i]-1.)+bary[i]*2.*dbary[i][1];
    } else {
      int j=i-3;
      ret[0] =
        4.*dbary[(3-j)%3][0]*bary[(4-j)%3]+
        4.*bary[(3-j)%3]*dbary[(4-j)%3][0];
      ret[1] =
        4.*dbary[(3-j)%3][1]*bary[(4-j)%3]+
        4.*bary[(3-j)%3]*dbary[(4-j)%3][1];
    }
    return ret;
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
