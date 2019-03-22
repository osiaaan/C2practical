#ifndef DISCRETEFUNCTION_HH_INCLUDED
#define DISCRETEFUNCTION_HH_INCLUDED

//- include vector implementation
#include "../common/vtkout.hh"
#include "../common/vector.hh"

/**
    @addtogroup DiscreteFunction
    @{
*/

/** \brief LinearBaseFunction implements the
    base function set for a piecewise linear Lagrange space.
    This class also represents the interface for all base functions
    sets as required as template arguments for the DiscreteFunction class.
*/
class LinearBaseFunction {
public:
  /** \brief Number of base functions
   */
  enum {locNrDof = 3};

  /** \brief constructor creating a LinearBaseFunction set
      \param[in] grid  grid that the base function set is based on
  */
  LinearBaseFunction(const GridType& grid) :
    _size(grid.size(2))
  {
    _point[0][0]=0.;
    _point[0][1]=0.;
    _point[1][0]=1.;
    _point[1][1]=0.;
    _point[2][0]=0.;
    _point[2][1]=1.;
    _gradphi[0][0]=-1.;
    _gradphi[0][1]=-1.;
    _gradphi[1][0]=1.;
    _gradphi[1][1]=0.;
    _gradphi[2][0]=0.;
    _gradphi[2][1]=1.;
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
    return element.index( 2, vx );
  }

  /** \brief returns true if local dof is located at boundary
      \param edge edge that is checked
      \param vx local vertex number of edge

      \return \b true if dof is on boundary, \b false otherwise
  */
  bool onBnd(const int edge, const int vx) const
  {
    if ( edge == 0 )
      return (vx < 2);
    else
      return ( vx == 2 || vx-1 == edge );
  }

  /** \brief evaluate base function
      \param[in]  i   number of base function to evaluate
      \param[in]  lambda  local coordinate where base function is evaluated

      \return evaluation of the i-th base function
  */
  double evaluate(const int i,
                  const LocalCoordType& lambda) const
  {
    double bary[3]={1.-lambda[0]-lambda[1],
              lambda[0],lambda[1]};
    return bary[i];
  }

  /** \brief evaluate gradient of base function
      \param[in]  i   number of base function to evaluate
      \param[in]  lambda  local coordinate where the gradient of the
                          base function is evaluated

      \return evaluation of the i-th base function gradient
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
  //! array with lagrange points in referece element
  LocalCoordType _point[locNrDof];
  //! array with gradients of basis functions (not really required)
  LocalCoordType _gradphi[locNrDof];
  //! dimension of discrete function space (set in constructor depending on
  //! grid
  int _size;
};

/** \brief Discrete Function implementation. A Discrete Function
    combines the DoF storage with a base function set
    \param BaseSet implementation class of the base function set
*/
template <class BaseSet>
class DiscreteFunction : public Vector<double>
{
  typedef Vector<double> BaseType;
  DiscreteFunction(const DiscreteFunction& );
public:
  //! type of base function set
  typedef BaseSet BaseFunctionSetType;

  /** \brief constructor creating discrete function
      \param[in]  grid  Grid for this discrete function
  */
  explicit DiscreteFunction(const GridType& grid)
   : BaseType(),
     grid_(grid),
     base_(grid)
  {
    // resize vector with number of DoFs
    this->resize( base_.size() );
  }

  /** \brief evaluate discrete function
      \param[in]  element element the discrete function is evaluated on
      \param[in]  lambda local coordinate in element for evaluation

      \return f(x)
   */
  double evaluate(const ElementType& element,
                  const LocalCoordType& lambda) const
  {
    double ret = 0;
    for (int i=0; i<BaseFunctionSetType::locNrDof; ++i)
    {
      ret += (*this)[base_.map(element,i)] * evaluate(element, i,lambda);
    }
    return ret;
  }

  /** \brief evaluate the gradient of this discrete function
      \param[in] element Element on which contains the coordinate to
             evaluate
      \param[in]  lambda local coodrinate of the evaluation point

      \return gradient of discrete function at given point
  */
  GlobalCoordType gradient(const ElementType& element,
                           const LocalCoordType& lambda) const
  {
    GlobalCoordType ret(0);
    GlobalCoordType grad(0);
    for (int i=0;i<BaseFunctionSetType::locNrDof;++i)
    {
      grad  = gradient(element,i,lambda);
      grad *= (*this)[base_.map(element,i)];
      ret  += grad;
    }
    return ret;
  }

  /** \brief evaluate one base function
      \param[in] element the element on which to evaluate the base function
      \param[in] i       number of base function to evaluate
      \param[in] lambda  local coordinate where base function is evaluated

      \return evaluation of the i-th base function
  */
  double evaluate(const ElementType& element,
                  const int i,
                  const LocalCoordType& lambda) const
  {
    return base_.evaluate(i, lambda);
  }

  /** \brief evaluate gradient of base function on element
      \param[in]  element  element to evaluate gradient on
      \param[in]  i   number of base function to evaluate
      \param[in]  lambda  local coordinate  where the gradient of the
       base function is evaluated

     \return evaluation of the i-th base function gradient
  */
  GlobalCoordType gradient(const ElementType& element,
                           const int i,
                           const LocalCoordType& lambda) const
  {
    GlobalCoordType ret;
    element.gradientGlobal( base_.gradientLocal(i,lambda),
                            ret);
    return ret;
  }

  /** \brief returns reference to base function set */
  const BaseFunctionSetType& baseFunctionSet () const {
    return base_;
  }

  /** \brief DoF mapping from element and local DoF number to global DoF
      \param[in]  element element the DoF is located on
      \param[in]  localDof local DoF number

      \return reference to DoF
  */
  double& dof(const ElementType& element, const int localDof)
  {
    return this->operator []( base_.map(element , localDof ) );
  }

  /** \brief return reference to corresponding grid
  */
  const GridType& grid () const { return grid_; }

protected:

  //! reference to grid
  const GridType& grid_;
  //! base function set
  BaseFunctionSetType base_;
};

/**
    @}
*/
#endif
