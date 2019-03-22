#ifndef QUADRATURE_HH_INCLUDED
#define QUADRATURE_HH_INCLUDED

//- system includes 
#include <vector> 

/** @addtogroup Quadratures 
    @{
*/

/** \brief Quadrature interface class. To implement a new quadrature 
    derive from this class and use the method \b addQuadraturePoint to add
    new quadrature points. 
*/
class Quadrature 
{
  // no copying 
  Quadrature(const Quadrature& ); 

protected:  
  /** \brief creating empty quadrature */
  Quadrature() 
   : points_(), 
     weights_() 
  {
  }

  /** \brief add (point,weight) pair to quadrature 
      \param[in] point local coordinates of quadrature point 
      \param[in] weight weight for quadrature point 
  */
  void addQuadraturePoint(const LocalCoordType& point,
                          const double weight)
  {
    points_.push_back( point );
    weights_.push_back( weight );
  }

public:
  /** \brief return number of quadrature points
      \return number of quadrature points */
  int nop() const 
  {
    assert( points_.size() == weights_.size() );
    return points_.size();
  }

  /** \brief return weight for quadrature point 
      \param[in] qp number of quadrature point that weight is returned for 

      \return weight for quadrature point \b qp 
  */
  double weight(const int qp) {
    return weights_[ qp ];
  }

  /** \brief return local coordinates for quadrature point 
      \param[in] qp number of quadrature point that local coordinates are 
             returned for 

      \return local coordinates for quadrature point \b qp 
  */
  const LocalCoordType& point(const int qp) const 
  {
    assert( qp >= 0 && qp < (int)points_.size() );
    return points_[ qp ];
  }

private:
  std::vector<LocalCoordType> points_;
  std::vector<double> weights_;
};


/** \brief Quadrature implementation for order = 1 (center of gravity)
 */
class PointQuadrature : public Quadrature 
{
public:
  //! constructor creating a quadrature 
  PointQuadrature() 
  {
    LocalCoordType point (1./3.);
    double weight = 0.5; 
    addQuadraturePoint( point, weight );
  }
};

/** \brief Quadrature implementation for order = 2 (center of edges)
 */
class LineQuadrature : public Quadrature 
{
public:
  //! constructor creating a quadrature 
  LineQuadrature() 
  {
    for(int i=0; i<3; ++i) 
    {
      LocalCoordType point (0.5);
      if( i == 0 ) point[1] = 0;
      if( i == 1 ) point[0] = 0;
      double weight = 1.0/3.0; 
      addQuadraturePoint( point, weight );
    }
  }
};
/** \brief Quadrature implementation for order = 3 
 */
class ThirdQuadrature : public Quadrature 
{
public:
  //! constructor creating a quadrature 
  ThirdQuadrature() 
  {
    for(int i=0; i<3; ++i) 
    {
      LocalCoordType point (0.);
      if (i == 1 ) point[0] = 1.;
      if (i == 2 ) point[1] = 1.;
      double weight = 3.0/120.; 
      addQuadraturePoint( point, weight );
    }
    for(int i=0; i<3; ++i) 
    {
      LocalCoordType point (0.5);
      if( i == 0 ) point[1] = 0;
      if( i == 1 ) point[0] = 0;
      double weight = 8.0/120.0; 
      addQuadraturePoint( point, weight );
    }
    LocalCoordType point (1./3.);
    double weight = 27./120.; 
    addQuadraturePoint( point, weight );
  }
};

/** 
  @}
*/

#endif
