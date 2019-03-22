#ifndef QUADRATURE5_HH_INCLUDED
#define QUADRATURE5_HH_INCLUDED

//- system includes 
#include <vector> 
#include "quadrature.hh"

/** @addtogroup Quadratures 
    @{
*/

/** \brief Quadrature implementation for order = 5
 */
class QuinticQuadrature : public Quadrature 
{
  void addQuadraturePoint(const double px, const double py,
                          const double weight)
  {
    LocalCoordType point;
    point[0] = px;
    point[1] = py;
    Quadrature::addQuadraturePoint( point, weight );
  }
public:
  //! constructor creating a quadrature 
  QuinticQuadrature() 
  {
    const double a1 = (6.0+sqrt(15.0))/21.0;
    const double a2 = (6.0-sqrt(15.0))/21.0;
    const double b1 = (9.0-2.0*sqrt(15.0))/21.0;
    const double b2 = (9.0+2.0*sqrt(15.0))/21.0;
    const double w1 = (155.0+sqrt(15.0))/2400.0;
    const double w2 = (155.0-sqrt(15.0))/2400.0;
    addQuadraturePoint( 1./3.,1./3., 9./80.);
    addQuadraturePoint( a1,a1, w1 );
    addQuadraturePoint( b1,a1, w1 );
    addQuadraturePoint( a1,b1, w1 );
    addQuadraturePoint( a2,a2, w2 );
    addQuadraturePoint( b2,a2, w2 );
    addQuadraturePoint( a2,b2, w2 );
  }
};

/** 
  @}
*/

#endif
