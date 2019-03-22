#ifndef GRID_HH_INCLUDED
#define GRID_HH_INCLUDED

//- system includes
#include <string>

//- local includes -- types from dune that are used here
#include "../common/types.hh"


// name space for grid stuff
namespace GRID {

  /**
      @addtogroup Grid
      @{
  */

//-----------------------------------------------------------------------
//  Forward declarations
//-----------------------------------------------------------------------
class Grid;

template <class GridImp>
class ElementIterator;

template <class GridImp>
class Element;

template <class GridImp, class IteratorImp, class ItemImp>
class Iterator;

template <class GridImp>
class IntersectionIterator;

template <class GridImp>
class Intersection;

//------------------------------------------------------------------------
//  Implementations
//------------------------------------------------------------------------
/** \brief The Element class provides all information needed on single
    elements.

    \param GridImp implementation of grid, here the Grid class.
*/
template <class GridImp>
class Element
{
public:
  //! type of this class
  typedef Element < GridImp > ThisType ;
  //! entity type that is used to build this element
  typedef typename GridImp :: EntityType EntityType;

  //! dimension of the element
  enum {
    dimension = EntityType :: dimension,          //!< dimension of the grid
    dimensionworld = EntityType :: dimensionworld //!< dimension of the grid
  };

  //! number of vertices that are used to build this element
  enum {
    nCorners = dimension + 1 //!< number of corners of the element
  };

private:
  // p0,p1,...,p(dim+1): corner coordinates of element
  GlobalCoordType p_[ nCorners ];

  // these are the rows in the reference map F_T
  JacobianType A_;

  // inverse matrix of A (transposed)
  JacobianType invA_T_;

  //! store det(DF_T) of the entity
  double det_;

  //! store 1/det(DF_T) of the entity
  double detinv_;

  //! store area of the entity
  double area_;

  //! store element index
  int elidx_;
  //! store edge indicies
  int edgeidx_[ nCorners ];
  //! store vertex indicies
  int vtxidx_[ nCorners ];
  //! store if edge on boundary
  mutable bool bnd_[ nCorners ];

public:
  /** \brief type of intersection iterator */
  typedef IntersectionIterator<GridImp> IntersectionIteratorType;

private:
  //! reference to grid
  const GridImp& grid_;

  //! internal intersection iterator implementation
  typedef typename GridImp :: IntersectionIteratorImp IntersectionIteratorImp;
  IntersectionIteratorImp iIt_;
  IntersectionIteratorImp iEnd_;

  typedef typename GridImp :: EntityPointerImp  EntityPointerImp;
  EntityPointerImp entityPointer_;

  //! false if bnd is not called
  mutable bool initBnd_;

  // friendship to ElementIterator
  friend class Iterator< GridImp, typename GridImp :: IteratorImp , ThisType > ;
  // friendship to Intersection
  friend class Intersection< GridImp >;

  // friendship to grid
  friend class Grid ;

  /** \brief return reference to stored entity pointer (needed for marking elements)
      \return reference to entity pointer
  */
  const EntityType& entity() const
  {
    return *entityPointer_;
  }

protected:
  /** \brief update to next element
      \param[in] enPtr dune entity that is wrapped
  */
  void update(const EntityType& entity )
  {
    // for safety reasons
    entityPointer_ = entity.template subEntity<0>(0) ;

    // copy coordinates
    typedef typename EntityType :: Geometry Geometry;
    const Geometry& geo = entity.geometry();
    for( int i=0; i<nCorners; ++i )
    {
      p_[i] = geo.corner(i);
    }

    // get intersection information
    iIt_  = grid_.gridView_.ibegin( entity );
    iEnd_ = grid_.gridView_.iend  ( entity );

    // get element index
    elidx_ = grid_.elementIndex( entity );

    // get edge and vertex index
    for(int i=0; i<nCorners; ++i)
    {
      edgeidx_[i] = grid_.edgeIndex( entity , i );
      vtxidx_ [i] = grid_.vertexIndex( entity , i );
    }

    // setup reference mapping matrix
    for(int i=0; i<dimensionworld; ++i)
      for(int j=0; j<dimension; ++j)
        A_[ i ][ j ] = p_[ j+1 ][ i ] - p_[ 0 ][ i ];

    // setup det
    det_ = geo.integrationElement(LocalCoordType(0));
    assert(det_ > 0);

    // calculate area of element
    area_ = 0.5 * det_;

    // calculate inverse det
    detinv_ = 1.0 / det_;

    // calculate inverse matrix using Cramers rule
    invA_T_ = geo.jacobianInverseTransposed(LocalCoordType(0));

    initBnd_ = false;
  }
  void initBnd() const
  {
    if (!initBnd_)
    {
      const IntersectionIteratorType enditi= iend();
      for(IntersectionIteratorType iti = ibegin(); iti != enditi; ++iti)
      {
        const Intersection<GridImp>& intersection = *iti;
        const int nummer = intersection.numberInSelf();
        bnd_[nummer] = intersection.boundary();
      }
      initBnd_ = true;
    }
  }

private:
  //! no copying allowed
  Element(const Element& );
public:
  /** \brief constructor ceating empty element */
  Element(const GridImp& grid)
    : grid_( grid ),
      iIt_( grid.initIntersection() ),
      iEnd_( grid.initIntersection() ),
      entityPointer_( grid.entityPointer() ),
      initBnd_(false)
  {
  }
  Element(const GridImp& grid,
          const EntityPointerImp& enPtr)
    : grid_( grid ),
      iIt_( grid.initIntersection() ),
      iEnd_( grid.initIntersection() ),
      entityPointer_( enPtr ),
      initBnd_(false)
  {
    update(*enPtr);
  }

  /** \brief return intersection iterator pointing to the first
      intersection */
  IntersectionIteratorType ibegin () const
  {
    return IntersectionIteratorType( grid_, iIt_ , iEnd_ );
  }

  /** \brief return intersection iterator pointing behind the last intersection */
  const IntersectionIteratorType iend () const
  {
    return IntersectionIteratorType( grid_, iEnd_ , iEnd_ );
  }

  /** \brief return coordinate of i-th corner
      \param i number of corner of element

      \return coordinate vector of i-th corner
  */
  const GlobalCoordType& operator [] (const int i) const
  {
    assert( i >= 0 && i < nCorners );
    return p_[i];
  }


  /** \brief return global index of subentites given local index
      \param codim codimension of subentity
      \param local local number of subentity
      \return global index
  */
  int index(int codim, int local) const
  {
    switch (codim)
    {
    case 0 : assert(local==0); return elidx_;
    case 1 : assert(0<=local && local < nCorners); return edgeidx_[ local ];
    case 2 : assert(0<=local && local < nCorners); return vtxidx_[ local ];
    }
    std::cout << "Wrong codimension " << codim << "used - aborting" << std::endl;
    abort();
  }

  /** \brief mapping from local to global coordinate
      \param local local coordinate \f$ \lambda \f$ (in reference element)

      \return global coordinate \f$ x \f$ (in element), i.e. \f$ x = F(\lambda) \f$
  */
  GlobalCoordType global(const LocalCoordType& local) const
  {
    // map local to global coordinates
    GlobalCoordType X ( p_ [ 0 ] );

    // multiply local with matrix A_
    for(int i=0; i<dimensionworld; ++i)
      X[i] += A_[i] * local;

    return X;
  }

  /** \brief mapping of gradient from reference element to element coordinates
      \param dphi \f$ \nabla \hat{\varphi}\f$ gradient with respect to the
             reference element
      \param gradient return value for the gradient with respect to the element:
              \f$\nabla \varphi = J_{E}^{-T} \nabla \hat{\varphi}\f$,
              where \f$J_E = DF_E\f$.
  */
  void gradientGlobal(const LocalCoordType& dphi,
                      GlobalCoordType& gradient) const
  {
    // map a vector given on reference element to world coordinates
    // apply dphi to inverse, transposed of A
    for(int i=0; i<dimensionworld; ++i)
    {
      gradient[ i ] = invA_T_[ i ] * dphi;
    }
  }

  /** \brief returns area of element
      \return area of element
  */
  double area() const
  {
    return area_;
  }

  /** \brief returns determinant of reference mapping
      \return \f$|det(J_{E})|\f$ the determinant of the reference mapping
  */
  double integrationElement() const
  {
    // return the integration factor
    return 2.0 * area_;
  }

  /** \brief returns if an edge is is part of the boundary
      \param edge edge number in element
      \return true if on boundary
  */
  bool onBoundary( int edge ) const
  {
    initBnd();
    return bnd_[edge];
  }

};

/** \brief The Intersection class provides all information for an
    intersection.

    \param GridImp implementation of grid, here the Grid class.
*/
template <class GridImp>
class Intersection
{
  typedef Intersection < GridImp > ThisType ;
  typedef typename GridImp :: EntityType EntityType;
  typedef typename GridImp :: IntersectionImp IntersectionImp;
  typedef typename GridImp :: IntersectionIteratorImp IntersectionIteratorImp;

  typedef Element< GridImp > ElementType;

  enum { dimension = EntityType :: dimension  };

private:
  // current element
  ElementType element_;

  // true if intersection is with boundary
  bool boundary_;
  // true if intersection is with other element
  bool neighbor_;

  int numberInSelf_;
  int numberInNeighbor_;

  LocalCoordType xInside_[2];
  LocalCoordType xOutside_[2];

  GlobalCoordType normal_;

  // friendship to ElementIterator
  friend class Iterator< GridImp, typename GridImp :: IntersectionIteratorImp , ThisType > ;

  /** \brief update to next element
      \param iIt dune intersection iterator that is wrapped
  */
  void update(const IntersectionImp& intersection)
  {
    boundary_ = intersection.boundary();
    neighbor_ = intersection.neighbor();
    numberInSelf_ = intersection.indexInInside();
    numberInNeighbor_ = intersection.indexInOutside();
    FaceCoordType tmp(0.5);
    normal_ = intersection.integrationOuterNormal(tmp);

    xInside_[0] = intersection.geometryInInside().corner(0);
    xInside_[1] = intersection.geometryInInside().corner(1);
    if (neighbor_)
    {
      xOutside_[0] = intersection.geometryInOutside().corner(0);
      xOutside_[1] = intersection.geometryInOutside().corner(1);
    }

    if( neighbor_ )
    {
      element_.update( *(intersection.outside()) );
    }
  }

private:
  // no copying allowed
  Intersection(const Intersection& );
public:
  /** \brief constructor creating intersection
      \param grid reference to Grid
  */
  Intersection(const GridImp& grid)
    : element_( grid )
  {
  }

  /** \brief returns true if intersection is with boundary
      \return \b true if intersection is with boundary, \b false
      otherwise
   */
  bool boundary () const { return boundary_; }

  /** \brief returns true if intersection is with an element
      \return \b true if intersection is with an element, \b false
      otherwise
   */
  bool neighbor () const { return neighbor_; }

  /** \brief returns local edge/face number in element the intersection
      iterator was started

      \return local edge/face number in self element
  */
  int numberInSelf() const { return numberInSelf_; }

  /** \brief returns local edge/face number in neighbor element, the
      return is only valid if neighbor() returns \b true

      \return local edge/face number in neighbor element
  */
  int numberInNeighbor() const
  {
    assert( neighbor() );
    return numberInNeighbor_;
  }

  /** \brief return coordinate in reference element of inside entity for some local
   * point on intersection

     \return point in local coordinates in reference element
  */
  LocalCoordType geometryInSelf( const FaceCoordType &x) const
  {
    LocalCoordType ret(xInside_[0]);
    ret *= (1.-x[0]);
    ret.axpy(x[0],xInside_[1]);
    return ret;
  }

  /** \brief return coordinate in reference element of outside entity for some local
   * point on intersection

     \return point in local coordinates in reference element
  */
  LocalCoordType geometryInNeighbor( const FaceCoordType &x) const
  {
    assert( neighbor() );
    LocalCoordType ret(xOutside_[0]);
    ret *= (1.-x[0]);
    ret.axpy(x[0],xOutside_[1]);
    return ret;
  }

  /** \brief return outer normal to intersection

     \return outer normal
  */
  const GlobalCoordType& outerNormal() const
  {
    return normal_;
  }

  /** \brief return neighbor element, the return is only valid if
      neighbor() returns \b true

      \return current neighbor element
  */
  const ElementType& outside () const
  {
    assert( neighbor() );
    return element_;
  }
};

//--------------------------------------------------------------------
//  Iterator Interface class
//--------------------------------------------------------------------
/** \brief Iterator interface for this code. This interface is used by
    the ElementIterator as well as for the IntersectionIterator.

    \param GridImp type of grid implementation, here Grid.
    \param IteratorImp type of the Dune iterator implementation, i.e.
           LeafIterator or LeafIntersectionIterator.
    \param ItemImp type of iterated item, i.e. Element or Intersection.

*/
template <class GridImp, class IteratorImp, class ItemImp>
class Iterator
{
  //! type of this class
  typedef Iterator<GridImp, IteratorImp, ItemImp> ThisType;

public:
  //! type of iterated entity
  typedef ItemImp ItemType;

  /** \brief constructor creating iterator */
  Iterator(const GridImp& grid,
               const IteratorImp& it,
               const IteratorImp& end)
    : grid_(grid),
      it_(it),
      end_(end),
      item_(grid_)
  {
    update();
  }

  /** \brief constructor copying iterator */
  Iterator(const ThisType& other)
    : grid_(other.grid_),
      it_(other.it_),
      end_(other.end_),
      item_(grid_)
  {
    update();
  }

  /** \brief assigment operator
    \param other another instance of an iterator
    \return reference to this object
  */
  ThisType& operator = (const ThisType& other)
  {
    it_ = other.it_;
    update();
    return *this;
  }

  /** \brief increment iterator
   *  \return reference to this
   */
  ThisType& operator ++ ()
  {
    assert( it_ != end_ );
    ++it_;
    update();
    return *this;
  }

  /** \brief equality operator
      \param other another instance of this operator
      \return true if this == other
  */
  bool operator == (const ThisType& other) const
  {
    return (it_ == other.it_);
  }

  /** \brief inequality operator
      \param other another instance of this operator
      \return true if this != other
  */
  bool operator != (const ThisType& other) const
  {
    return (it_ != other.it_);
  }

  /** \brief return reference to current entity
      \return reference to current entity
   */
  const ItemType& operator * () const
  {
    assert( it_ != end_ );
    return item_;
  }

  /** \brief return reference to current entity
      \return reference to current entity
   */
  const ItemType* operator -> () const
  {
    assert( it_ != end_ );
    return & item_;
  }

protected:
  /** \brief update current entity */
  void update()
  {
    if( it_ != end_ )
    {
      item_.update( *it_ );
    }
  }

  //! reference to grid
  const GridImp& grid_;

  //! my iterator from dune
  IteratorImp it_;
  //! end iterator
  const IteratorImp end_;

  //! item implementation
  ItemType item_;
};

//-------------------------------------------------------------
// Element Iterator
//-------------------------------------------------------------
/** \brief Element Iterator class.
    \param GridImp implementation of grid, here Grid.
*/
template <class GridImp>
class ElementIterator
: public Iterator< GridImp,
                       typename GridImp :: IteratorImp,
                       Element< GridImp >
                     >
{
  //! type of this class
  typedef ElementIterator<GridImp> ThisType;
  typedef Iterator< GridImp,
                        typename GridImp :: IteratorImp,
                        Element< GridImp >
                      > BaseType;

  //! internal iterator type
  typedef typename GridImp :: IteratorImp IteratorImp;

public:
  /** \brief constructor creating iterator */
  ElementIterator(const GridImp& grid,
                  const IteratorImp& it,
                  const IteratorImp& end)
    : BaseType( grid, it, end )
  {
  }

  /** \brief constructor copying iterator */
  ElementIterator(const ThisType& other)
    : BaseType( other )
  {
  }

};

//-------------------------------------------------------------
// Intersection Iterator
//-------------------------------------------------------------
/** \brief Intersection Iterator class.
    \param GridImp implementation of grid, here Grid.
*/
template <class GridImp>
class IntersectionIterator
: public Iterator< GridImp,
                       typename GridImp :: IntersectionIteratorImp,
                       Intersection< GridImp >
                     >
{
  //! type of this class
  typedef IntersectionIterator<GridImp> ThisType;
  typedef Iterator< GridImp,
                        typename GridImp :: IntersectionIteratorImp,
                        Intersection< GridImp >
                      > BaseType;

  //! internal iterator type
  typedef typename GridImp :: IntersectionIteratorImp  IntersectionIteratorImp;

public:
  /** \brief constructor creating iterator */
  IntersectionIterator(const GridImp& grid,
                       const IntersectionIteratorImp& it,
                       const IntersectionIteratorImp& end)
    : BaseType( grid, it, end )
  {
  }

  /** \brief constructor copying iterator */
  IntersectionIterator(const ThisType& other)
    : BaseType( other )
  {
  }

};


//--------------------------------------------------------------------------
/** \brief The class Grid is the interface for a grid in this practical
    code. It provides methods to get the current size of the grid and
    also method to refine the grid.
*/
//--------------------------------------------------------------------------
//- --Grid
class Grid
{
  //! type of this class
  typedef Grid ThisType;

  //! type of grid part implementation
  typedef Dune :: GridPartView GridViewImp;
  //! type of dune grid implementation
  typedef GridViewImp :: Grid GridImp;
  //! type of index set
  typedef GridViewImp :: IndexSet IndexSetImp;

public:
  // internal iterator type
  typedef GridViewImp :: Codim <0> :: Iterator IteratorImp;
  // intenal entity type
  typedef GridImp :: Codim<0> :: Entity EntityType;
  // internal entity pointer type
  typedef GridImp :: Codim<0> :: EntityPointer EntityPointerImp;
  // internal iterator type
  typedef GridViewImp :: IntersectionIterator IntersectionIteratorImp;

  // internal iterator type
  typedef IntersectionIteratorImp :: Intersection IntersectionImp;

  friend class GRID::Element<ThisType> ;
  friend class GRID::ElementIterator<ThisType> ;
  friend class GRID::IntersectionIterator<ThisType> ;

public:
  /** \brief dimension of the grid */
  enum {
    dimension = GridImp :: dimension //!< dimension of the grid
  };

  /** \brief codimension --> element type */
  enum { elements = 0, //!< element identifier
         faces = 1,    //!< face identifier
         edges = dimension - 1, //!< edge identifer
         vertices = dimension   //!< vertex identifier
  };

  /** \brief type of the exported element iterator */
  typedef ElementIterator<  ThisType > ElementIteratorType;

  /** \brief cosntructor taking name of macro grid file
      \param name filename of macro grid file in DGF format
  */
  Grid(const std::string name)
    : gridPointer_(name)
    , gridView_( * gridPointer_ )
    , indexSet_( gridView_.indexSet() )
  {
  }

protected:
  /** \brief type of element */
  typedef ElementIteratorType :: ItemType ElementType;

  /** \brief return global vertex number of local vertex on entity
      \param entity Entity the information is needed for
      \param localVertex local vertex number

      \return global vertex number
  */
  int vertexIndex(const EntityType& entity, const int localVertex) const
  {
    return indexSet_.subIndex(entity, localVertex, vertices);
  }

  /** \brief return global edge number of local edge on entity
      \param entity Entity the information is needed for
      \param localEdge local edge number

      \return global edge number
  */
  int edgeIndex(const EntityType& entity, const int localEdge) const
  {
    return indexSet_.subIndex (entity, localEdge, edges);
  }

  /** \brief return global face number of local face on entity
      \param entity Entity the information is needed for
      \param localFace local face number

      \return global face number
  */
  int faceIndex(const EntityType& entity, const int localFace) const
  {
    return indexSet_.subIndex (entity, localFace, faces );
  }

  /** \brief return index of entity
      \param entity Entity the index information is needed for

      \return index of entity
  */
  int elementIndex(const EntityType& entity) const
  {
    return indexSet_.index( entity );
  }

public:
#if 0
  /** \brief return number of vertices
      \return number of vertices of the grid
  */
  int nVertices() const {
    return this->size( vertices );
  }

  /** \brief returns number of faces of the grid
      \return number of faces of the grid (same as nEdges in 2d)
   */
  int nFaces() const
  {
    return this->size( faces );
  }

  /** \brief returns number of edges of the grid
      \return number of edges of the grid (same as nFaces in 2d)
   */
  int nEdges() const
  {
    return this->size( edges );
  }

  /** \brief returns the number of elements
      \return number of elements of the grid
  */
  int nElements() const {
    return this->size(elements);
  }
#endif

  /** \brief returns number of elements for given codimension
      \param codim codimenion size information is required for
  */
  int size ( const int codim ) const
  {
    return indexSet_.size( codim );
  }

  /** \brief mark element for refinement or coarsening
      \param marker 1 for refinement and -1 for coarsening
      \param element Element that should be marked
  */
  void mark(const int marker, const ElementType& element)
  {
    gridPointer_->mark( marker , element.entity() );
  }

  /** \brief do adaptation of a previously marked grid */
  void adapt()
  {
    gridPointer_->preAdapt();
    gridPointer_->adapt();
    gridPointer_->postAdapt();
  }

  /** \brief globalRefine refines all elements of the grid */
  void globalRefine()
  {
    globalRefine( 1 );
  }

  /** \brief do global refinement several times
      \param refines number of global refinements to be done
  */
  void globalRefine(const int refines)
  {
    gridPointer_->globalRefine(refines);
  }

  /** \brief return iterator pointing to the first element of the grid
      \return iterator pointing to the first element of the grid
  */
  ElementIteratorType begin() const
  {
    return ElementIteratorType ( *this,
                                 gridView_.begin<0> (),
                                 gridView_.end<0> ()
                               );
  }

  /** \brief return iterator pointing behind the last element of the grid
      \return iterator pointing behind the last element of the grid
  */
  const ElementIteratorType end() const
  {
    return ElementIteratorType ( *this,
                                 gridView_.end<0> (),
                                 gridView_.end<0> ()
                               );
  }

  /** \brief conversion operator to turn into a Dune grid (needed for
      visualization )
  */
  const GridImp& operator() ()
  {
    return * gridPointer_;
  }

  //! \brief return dune grid view for visualization output
  const GridViewImp& gridView() const
  {
    return gridView_;
  }

  //! \brief return index set for visualization output
  const IndexSetImp& indexSet() const
  {
    return indexSet_;
  }

protected:
  /** \brief initializer method for intersection iterators since they
      might not have a default constructor */
  IntersectionIteratorImp initIntersection() const
  {
    return gridView_.iend( * (gridView_.begin<0> () ) );
  }


  /** \brief initializer method for entity pointers since they
      might not have a default constructor */
  EntityPointerImp entityPointer() const
  {
    return gridView_.begin<0> ()->subEntity<0>(0);
  }

  //! grid pointer holding instance of dune grid
  Dune::GridPtr<GridImp> gridPointer_;
  //! grid part
  GridViewImp gridView_;
  //! reference to index set
  const IndexSetImp& indexSet_;

}; // end class Grid

  /**
     @}
  */

} // end namespace GRID

/**
    @addtogroup Grid
    @{
*/

//---------------------------------------------------------------------
//  global typedefs for further use in this code
//---------------------------------------------------------------------
//! type of grid, here GRID::Grid
typedef GRID :: Grid GridType;
//! type of the grids element iterator, here GRID::ElementIterator
typedef GridType :: ElementIteratorType ElementIteratorType;
//! type of element, here GRID::Element
typedef ElementIteratorType :: ItemType ElementType;
//! type of intersection iterator, here GRID::IntersectionIterator
typedef ElementType :: IntersectionIteratorType IntersectionIteratorType;
//! type of intersection, here GRID::Intersection
typedef IntersectionIteratorType :: ItemType IntersectionType;

//---------------------------------------------------------------------
//  Output for visualization with Paraview
//---------------------------------------------------------------------
/** \brief output method for writing \b Paraview vtu file containing grid
    \param[in] grid  reference to grid which is written
    \param[in] name optional parameter for naming output file
*/
inline static void output(GridType& grid, const std::string name = "grid" )
{
  Dune :: OutputType output(grid.gridView());
  output.write( name.c_str() );
}

/**
    @}
*/
#endif // GRID_HH_INCLUDED
