#ifndef DISCRETEVECTORFIELD2D_HH_INCLUDED
#define DISCRETEVECTORFIELD2D_HH_INCLUDED

//- include vector implementation
#include "vtkout.hh"
#include "../common/vector.hh"

/**
    @addtogroup DiscreteVectorfield
    @{
*/

/** \brief Discrete Vectorfield implementation. A Discrete Vectorfield
    combines the DoF storage for each dimension with a base function set
    \param BaseSet implementation class of the base function set
*/
template <class BaseSet>
class DiscreteVectorfield : public Vector<double>
{
  typedef Vector<double> BaseType;
  DiscreteVectorfield(const DiscreteVectorfield& );
public:
  //! type of base function set
  typedef BaseSet BaseFunctionSetType;

  /** \brief constructor creating discrete vectorfield
      \param[in]  grid  Grid for this discrete vectorfield
  */
  explicit DiscreteVectorfield(const GridType& grid)
   : BaseType(),
     grid_(grid),
     base_(grid)
  {
    // resize vector with dow times number of DoFs
    this->resize( grid_.dimension * base_.size() );
  }

// BEGIN EDIT

  /** \brief evaluate discrete-vector function
      \param[in]  element element the discrete function is evaluated on
      \param[in]  lambda local coordinate in element for evaluation

      \return f(x)
   */
  LocalCoordType evaluate(const ElementType& element,
                  const LocalCoordType& lambda) const
  {
    LocalCoordType vec(0);
    for (int i=0; i<BaseFunctionSetType::locNrDof; ++i)
    {
      /*
      Note to self:
      Our vector field is of the following form
      A =  |* *|
           |* *|
           |* *|
      However it will be represented as a vector
      | * * * * * * |
      Where evey even entry represents a value from the first column of A
      and every odd entry represents a value from the second column of A
      */
      vec[0] += (*this)[base_.map(element,i)] * evaluate(element, i,lambda);
      vec[1] += (*this)[base_.map(element,i) + base_.size()] * evaluate(element, i,lambda);
    }
    return vec;
  }

  /** \brief evaluate the gradient of this discrete function
      \param[in] element Element on which contains the coordinate to
             evaluate
      \param[in]  lambda local coodrinate of the evaluation point

      \return gradient of discrete function at given point
  */
  JacobianType gradient(const ElementType& element,
                           const LocalCoordType& lambda) const
  {
    GlobalCoordType ret1(0);
    GlobalCoordType ret2(0);
    GlobalCoordType grad1(0);
    GlobalCoordType grad2(0);

    JacobianType M(0);

    for (int i=0;i<BaseFunctionSetType::locNrDof;++i)
    {
      //Note to self:
      //Going through the gradients of the basis function
      //and multiplying them by the evaluation of the
      //the u_h at that point

      //base_.map(element,i) gives the global coordinate of vetex i in
      //'element'
      grad1  = gradient(element,i,lambda);
      grad1 *= (*this)[base_.map(element,i)];
      ret1  += grad1;

      grad2  = gradient(element,i,lambda);
      //base_ rep. the nodal basis on the reference triangle
      //We evaluate them at the point and then map to triangle
      grad2 *= (*this)[base_.map(element,i) + base_.size()];
      ret2  += grad2;
    }

    //Adding what is the Jacobian Matrix!
    for(int j = 0; j < 2 ; ++j)
    {
      M[0][j] = ret1[j];
    }
    for(int j = 0; j < 2 ; ++j)
    {
      M[1][j] =ret2[j];
    }

    return M;
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
  double& dof(const ElementType& element, const int localDof, int factor)
  {
    // Evaluates u_h[@ global numbering of locDof]

    return this -> operator[] (base_.map(element, localDof) + base_.size()*factor );
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
