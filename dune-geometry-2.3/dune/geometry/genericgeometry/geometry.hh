// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GENERICGEOMETRY_GEOMETRY_HH
#define DUNE_GENERICGEOMETRY_GEOMETRY_HH

#warning This header and the code it contains is deprecated.  If you need functionality \
         similar to BasicGeometry, please use the MultiLinearGeometry class.

#include <dune/common/typetraits.hh>
#include <dune/common/nullptr.hh>

#include <dune/geometry/genericgeometry/mappingprovider.hh>
#include <dune/geometry/genericgeometry/geometrytraits.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    /** \addtogroup GenericGeometry
     *
     *  \section General
     *
     *  Based on a recursive definition of the reference elements, a generic
     *  implementation of Dune::Geometry is provided. The class used for the
     *  implementation of the Dune::Geometry engine is
     *  GenericGeometry::BasicGeometry.
     *
     *  The BasicGeometry class takes a template argument Traits specifying
     *  details of the reference mapping implementation and some performance
     *  settings. A default implementation for this class is
     *  GenericGeometry::DefaultGeometryTraits. The traits class must contain
     *  the same types as this default implementation.
     *
     *  To conform with the Dune::Geometry engine, two further classes are
     *  provided: GenericGeometry::Geometry and GenericGeometry::LocalGeometry.
     *  To use these classes instead of GenericGeometry::BasicGeometry, the
     *  traits classes
     *  \code
     *  template< class Grid> GenericGeometry::GlobalGeometryTraits<Grid>
     *  template< class Grid> GenericGeometry::LocalGeometryTraits<Grid>
     *  \endcode
     *  have to be specialized. These classes are simply passed as Traits
     *  argument to GenericGeometry::BasicGeometry.
     *
     *  The reference mapping for a given topology type is given by
     *  Mapping<Topology>::type in the traits class. Here, Topology is one of
     *  the generic topology classes GenericGeometry::Point,
     *  GenericGeometry::Prism, GenericGeometry::Pyramid.
     *  An interface for the mapping is provided by GenericGeometry::Mapping.
     *  The implementation of this interface must have constructors taking a
     *  single argument. The constructor of GenericGeometry::BasicGeometry
     *  looks as follows:
     *  \code
     *  template< class CoordVector >
     *  BasicGeometry ( const GeometryType &type, const CoordVector &coords );
     *  \endcode
     *  Its first argument, <em>type</em>, specifies the type of the reference
     *  element (as a Dune::GeometryType). The second argument, <em>coords</em>
     *  is passed directly to the constructor of the mapping implementation.
     *  The most prominent implementation of GenericGeometry::Mapping is
     *  GenericGeometry::CornerMapping. It provides a polynomial interpolation
     *  of the entity's corners with minimal degree. In this case,
     *  <em>coords</em> represents the entity's corners.
     *
     *  \section Simple Usage
     *  To add first order Lagrange type geometries to a grid implementation
     *  the following steps suffice:
     *  - Overload the traits classes
     *    \code
            template<>
            struct GenericGeometry::GlobalGeometryTraits< MyGrid >
            : public GenericGeometry::DefaultGeometryTraits
                     <MyGrid::ctype,MyGrid::dimension,MyGrid::dimworld>
            {};
            template<>
            struct GenericGeometry::LocalGeometryTraits< MyGrid >
            : public GenericGeometry::DefaultGeometryTraits
                     <MyGrid::ctype,MyGrid::dimension,MyGrid::dimworld>
            {};
     *    \endcode
     *    Note that these classes are default implementations which should cover
     *    all cases but are in a specific situation far from the optimal choice.
     *    For example, an increase
     *    in efficiency can be achieved for grids with a fixed element type
     *    (set hybrid to false and set the topologyId variable)
     *    or for grids with only affine transformations - which in the case of
     *    the local geometries is often true - the last template
     *    argument (default false) can be used to switch to mappings which are
     *    assumed to always be affine (no checking done).
     *  - Add to the GridFamily::Traits::Codim<codim> structure:
     *    \code
            typedef Dune :: Geometry
              < dimension-codim, dimensionworld, const MyGrid,
                Dune :: GenericGeometry :: Geometry > Geometry;
            typedef Dune :: Geometry
              < dimension-codim, dimension, const MyGrid,
                Dune :: GenericGeometry :: LocalGeometry > LocalGeometry;
     *    \endcode
     *  - Both geometries can be build by calling the constructor taking
     *    a DUNE grid type and an instance of an arbitrary class with a method
            \code
              const FieldVector<  ctype, dimensionworld >& operator[](unsigned int i);
            \endcode
     *    The references returned must remain valid during the whole life span
     *    of the geometry.
     *  - In MyGrid::Entity<0> the following methods can then be easily implemented:
     *    - geometry(): this requires the knowledge of the dune geometry type of the entity
     *      and the coordinates of the corner points.
     *    - geometryInFather(): The corner points for each child in the
     *      reference element of the father can be used to construct the local geometry -
     *      note that this geometry is mostly affine and these geometries can be
     *      precomputed and stored.
     *    .
     *  - For the Dune::Intersection class the geometries the following implementations for the geometries can be used:
     *    - intersectionGlobal(): can be implemented in the same way as the geometry of the
     *      entity using the coordinates of the corners of the intersection.
     *      Alternatively, in the case of a conform intersection,
     *      the class GenericGeometry::Geometry provides a possibility
     *      to construct traces of a given geometry, e.g., a reference mapping
     *      restricted to a codimension one subentities of the reference
     *      element. This is achieved by calling the constructor on the
     *      GenericGeometry::Geometry class (with the codim template equal to
     *      one) passing a codimension zero geometry implementation and the number of the
     *      codimension one subentity.
     *      \code
            GenericGeometry::Geometry<myGridDim-1,myWorldDim,MyGrid>
                             (inside->geometry(),numberInSelf());
     *      \endcode
     *    - intersectionInSelf()/intersectionInNeighbor():
     *      A similar strategy as described above for the intersectionGlobal
     *      can also be used for the geometry mapping to the codimension zero
     *      reference element. Either the corners of the intersection in
     *      local coordinates can be used in the construction of the local
     *      geometries, or (for conform intersections) the traces can be used,
     *      passing an identity mapping as codimension zero geometry.
     *      The GenericGeometry::GenericReferenceElement provides these
     *      mappings directly via the template method
     *      GenericGeometry::GenericReferenceElement::mapping.
     *      The return value of this method can be directly used to construct
     *      a GenericGeometry::Geometry instance:
     *      \code
            typedef GenericReferenceElementContainer<ctype,myGridDim> RefElementContType;
            RefElementContType refElemCont;
            const RefElementContType::value_type& refElem=refElemCont(insideGeometryType);
            GenericGeometry::Geometry<myGridDim-1,myGridDim,MyGrid>(refElem.mapping(numberInSelf()));
     *      \endcode
     *    - integrationOuterNormal(): the generic geometry implementation provides a method
     *      to compute the integration outer normals, so that the following code
     *      fragment can be used:
            \code
               typedef typename Grid :: template Codim< 0 > :: Geometry Geometry;
               const Geometry &geo = inside()->geometry();
               FieldVector< ctype, dimension > x( intersectionSelfLocal().global( local ) );
               return Grid :: getRealImplementation( geo ).normal( numberInSelf(), x );
            \endcode
     *    .
     *  - To add geometries for subentitiies of codim>0
     *    given a entity en of codimension zero and the subentity number subNr :
     *    - geometry: the geometry can be constructed by the following line of code
            \code
            GenericGeometry::Geometry<myGridDim-codim,myWorldDim,MyGrid>
                             (en.geometry(),subNr);
            \endcode
     *    .
     *  .
     *
     */

    // BasicGeometry
    // -------------

    /** \ingroup GenericGeometry
     *  \brief   generic implementation of DUNE geometries
     *
     *  This class is provides a generic implementation of a DUNE geometry.
     *
     *  Parameters shared by all codimensions are summarized in one class
     *  parameter called Traits. As a default traits class, the class
     *  DefaultGeometryTraits can be used.  Alternatively, the user can
     *  provide hand-written traits classes (which may, if that helps,
     *  derive from DefaultGeometryTraits).  Such classes have to provide
     *  the following fields:
     *  \code
     *  template< My_Template_Parameters >
     *  struct MyGeometryTraits
     *  {
     *    // ctype is the type used for coordinate coefficients
     *    typedef DuneCoordTraits< ctype > CoordTraits;
     *
     *    // Dimension of the space the geometry maps into
     *    static const int dimWorld = ...;
     *
     *    //   hybrid   [ true if reference element type is a run-time parameter ]
     *    static const bool hybrid = ...;
     *
     *    //   topologyId [ reference element type, only needed if it is not a run-time parameter ]
     *    // In this example: a dim-dimensional simplex
     *    // static const unsigned int topologyId = SimplexTopology< dim >::type::id;
     *
     *    // explained below
     *    template< class Topology >
     *    struct Mapping
     *    {
     *      typedef CornerMapping< CoordTraits, Topology, dimWorld > type;
     *    };
     *
     *    // Caching behavior
     *    struct Caching
     *    {
     *      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
     *    };
     *  };
     *  \endcode
     *
     *  The structure specifying the reference mapping is
     *  Traits::Mapping::type. An example implementation
     *  is the GenericGeometry::CornerMapping which defines
     *  the simple mapping taking corners of the reference
     *  elements to corner of the entity in space.
     *
     *  The central reference mapping specified by Traits::Mapping::type
     *  requires a constructor taking a single argument.
     *  The GenericGeometry::BasicGeometry has a constructor with one template
     *  argument which is passed on to the constructor of the reference mapping.
     *  The interface for the this class is GenericGeometry::Mapping.
     *
     *  To increase the efficiency of the geometry
     *  implementation, different strategies for
     *  the caching of parts of the geometry data
     *  is provided. The specifics are given
     *  by the structure Traits::Caching. Possible
     *  values are:
     *  - ComputeOnDemand:    use caching if method called using barycenter
     *  - PreCompute:         use caching in constructor using barycenter
     *  .
     *
     *  \note This class cannot be used directly as an implementation of
     *        Dune::Geometry. Its template parameter list differs from what
     *        is expected there from the engine.
     *        One of the following derived classes
     *        can be used instead:
     *        - Dune::GenericGeometry::Geometry
     *        - Dune::GenericGeometry::LocalGeometry
     *        .
     */
    template< int mydim, class Traits >
    class BasicGeometry
    {
      typedef typename Traits :: CoordTraits CoordTraits;

      /** \brief Be friend with other instantiations of the same class */
      template< int, class > friend class BasicGeometry;

    public:

      /** \brief The dimension of the parameter space of this geometry */
      static const int mydimension = mydim;

      /** \brief The dimension of the world space of this geometry */
      static const int coorddimension = Traits :: dimWorld;

      /** \brief Type used for coordinate components */
      typedef typename CoordTraits :: ctype ctype;

      /** \brief Type used for parameter coordinates */
      typedef FieldVector< ctype, mydimension > LocalCoordinate;

      /** \brief Type used for world coordinates */
      typedef FieldVector< ctype, coorddimension > GlobalCoordinate;

    private:
      dune_static_assert( (0 <= mydimension), "Geometry dimension must be nonnegative." );

      template< bool >
      struct Hybrid
      {
        typedef VirtualMappingFactory< mydimension, Traits > MappingFactory;
      };

      template< bool >
      struct NonHybrid
      {
        static const int topologyId = Traits::template hasSingleGeometryType< mydimension >::topologyId;
        typedef typename GenericGeometry::Topology< topologyId, mydimension >::type Topology;
        typedef GenericGeometry::NonHybridMappingFactory< Topology, Traits > MappingFactory;
      };

      static const bool hybrid = !Traits::template hasSingleGeometryType< mydimension >::v;

    protected:
      typedef typename conditional< hybrid, Hybrid< true >, NonHybrid< false > >::type::MappingFactory MappingFactory;
      typedef typename MappingFactory::Mapping Mapping;

    public:
      /** \brief Type used for Jacobian matrices
       *
       * \note This is not a FieldMatrix but a proxy type that can be assigned
       *       to a FieldMatrix.
       */
      typedef typename Mapping::JacobianTransposed JacobianTransposed;
      /** \brief Type used for Jacobian matrices
       *
       * \note This is not a FieldMatrix but a proxy type that can be assigned
       *       to a FieldMatrix.
       */
      typedef typename Mapping::JacobianInverseTransposed Jacobian;
      // for cenvencience, Jacobian is the name of the type in the geometry interface
      typedef Jacobian JacobianInverseTransposed;

    public:
      /** \brief Default constructor
       */
      BasicGeometry ()
        : mapping_( nullptr )
      {}

      /** \brief Constructor using a GeometryType and a list of corner coordinates */
      template< class CoordVector >
      BasicGeometry ( const GeometryType &type, const CoordVector &coords )
      {
        assert(type.dim() == mydim);
        mapping_ = MappingFactory::construct( type.id(), coords, mappingStorage_ );
      }

      /** \brief Constructor using a vector of corner coordinates and the dimension
       *  \note the geometry type is guessed from the number of vertices, thus this will only work up to dim 3
       */
      template< class CoordVector >
      BasicGeometry ( const CoordVector &coords )
      {
        GeometryType type;
        type.makeFromVertices( mydim, coords.size() );
        mapping_ = MappingFactory::construct( type.id(), coords, mappingStorage_ );
      }

      /** \brief obtain a geometry for a subentity
       *
       *  Assume that we have a geometry for some entity d-dimensional E.
       *  This method can provide a geometry for the i-th subentity of E
       *  (of codimension d - mydimension).
       *
       *  \note This method can be more efficient than just building up the
       *        geometry for the subentity. For example, the subgeometry
       *        automatically inherits affinity.
       *
       *  \param[in]  father  geometry of entity \em E
       *  \param[in]  i       number of the subentity (in generic numbering)
       */
      template< int fatherdim >
      BasicGeometry ( const BasicGeometry< fatherdim, Traits > &father, int i )
      {
        const unsigned int codim = fatherdim - mydim;
        mapping_ = father.mapping_->template trace< codim >( i, mappingStorage_ );
      }

      /** \brief Copy constructor */
      BasicGeometry ( const BasicGeometry &other )
        : mapping_( other.mapping_ ? other.mapping_->clone( mappingStorage_ ) : nullptr )
      {}

      /** \brief Destructor */
      ~BasicGeometry ()
      {
        if( mapping_ )
          mapping_->~Mapping();
      }

      /** \brief Assignment from other BasicGeometry */
      const BasicGeometry &operator= ( const BasicGeometry &other )
      {
        if( mapping_ )
          mapping_->~Mapping();
        mapping_ = (other.mapping_) ? other.mapping_->clone( mappingStorage_ ) : nullptr;
        return *this;
      }

      /** \brief bool cast
       *
       *  Like a pointer, a BasicGeometry casts to <b>true</b> if and only if
       *  it is properly initialized.
       *  If a geometry casts to <b>false</b>, none of the interface methods
       *  may be called.
       */
      operator bool () const
      {
        return bool( mapping_ );
      }

      /** \brief Return the topological type of this geometry */
      GeometryType type () const
      {
        return mapping_->type();
      }

      /** \brief Return the number of corners */
      int corners () const
      {
        return mapping_->numCorners();
      }

      /** \brief Return the world coordinates of the i-th corner */
      GlobalCoordinate corner ( const int i ) const
      {
        return mapping_->corner( i );
      }

      /** \brief Map local to global coordinates */
      GlobalCoordinate global ( const LocalCoordinate &local ) const
      {
        return mapping_->global( local );
      }

      /** \brief Map global to local coordinates */
      LocalCoordinate local ( const GlobalCoordinate &global ) const
      {
        return mapping_->local( global );
      }

      /** \brief return center of element */
      GlobalCoordinate center () const
      {
        return mapping_->center();
      }

      /** \brief Return true if this is an affine geometry */
      bool affine () const
      {
        return mapping_->affine();
      }

      /** \brief Return the factor \$|det F|\$ that appears in the integral transformation formula */
      ctype integrationElement ( const LocalCoordinate &local ) const
      {
        return mapping_->integrationElement( local );
      }

      /** \brief Return the volume of the element */
      ctype volume () const
      {
        return mapping_->volume();
      }

      /** \brief Compute the transpose of the Jacobian matrix of the
       *         transformation from the reference element into the world
       *         space
       */
      const JacobianTransposed &jacobianTransposed ( const LocalCoordinate &local ) const
      {
        return mapping_->jacobianTransposed( local );
      }

      /** \brief Compute the transpose of the inverse Jacobian matrix of the transformation
          from the reference element into the world space */
      const JacobianInverseTransposed &jacobianInverseTransposed ( const LocalCoordinate &local ) const
      {
        return mapping_->jacobianInverseTransposed( local );
      }

    private:

      /** \brief Always points to mappingStorage_, but has the correct type */
      Mapping* mapping_;

      /** \brief A chunk of raw memory storing the actual object
       *
       * We don't know its type, but we don't want to do classical
       * dynamic polymorphism, because heap allocation is expensive.
       */
      char mappingStorage_[ MappingFactory::maxMappingSize ];
    };



    // Geometry
    // --------

    /** \class   Geometry
     *  \ingroup GenericGeometry
     *  \brief   generic implementation of a DUNE (global) geometry
     *
     *  Geometry inherits all its features from BasicGeometry. It only adds
     *  GlobalGeometryTraits< Grid > as Traits parameter to the template
     *  parameter list.
     *
     * \tparam mydim Dimension of the entity
     * \tparam cdim Dimension of the coordinate space
     * \tparam Grid The grid this geometry will be used in
     */
    template< int mydim, int cdim, class Grid >
    class Geometry
      : public BasicGeometry< mydim, GlobalGeometryTraits< Grid > >
    {
      typedef BasicGeometry< mydim, GlobalGeometryTraits< Grid > > Base;

    protected:
      typedef typename Base::Mapping Mapping;

    public:

      Geometry ()
      {}

      /** \brief Copy constructor from another geometry */
      template< class Geo >
      explicit Geometry ( const Geo &geo )
        : Base( geo.type(), geo, geo.affine() )
      {}

      /** \brief Constructor with a GeometryType and a set of coordinates */
      template< class CoordVector >
      Geometry ( const GeometryType &type, const CoordVector &coords )
        : Base( type, coords )
      {}

      /** \todo Please doc me! */
      template< int fatherdim >
      Geometry ( const Geometry< fatherdim, cdim, Grid > &father, int i )
        : Base( father, i )
      {}
    };



    // LocalGeometry
    // -------------

    /** \class   LocalGeometry
     *  \ingroup GenericGeometry
     *  \brief   generic implementation of a DUNE (local) geometry
     *
     *  LocalGeometry inherits all its features from BasicGeometry. It only adds
     *  LocalGeometryTraits< Grid > as Traits parameter to the template
     *  parameter list.
     *
     * \tparam mydim Dimension of the entity
     * \tparam cdim Dimension of the coordinate space
     * \tparam Grid The grid this geometry will be used in
     */
    template< int mydim, int cdim, class Grid >
    class LocalGeometry
      : public BasicGeometry< mydim, LocalGeometryTraits< Grid > >
    {
      typedef BasicGeometry< mydim, LocalGeometryTraits< Grid > > Base;

    protected:
      typedef typename Base::Mapping Mapping;

    public:
      /** \brief Copy constructor from another geometry */
      template< class Geo >
      explicit LocalGeometry ( const Geo &geo )
        : Base( geo.type(), geo, geo.affine() )
      {}

      /** \brief Constructor with a GeometryType and a set of coordinates */
      template< class CoordVector >
      LocalGeometry ( const GeometryType &type, const CoordVector &coords )
        : Base( type, coords )
      {}

      /** \todo Please doc me! */
      template< int fatherdim >
      LocalGeometry ( const Geometry< fatherdim, cdim, Grid > &father, int i )
        : Base( father, i )
      {}
    };

  }

}

#endif // #ifndef DUNE_GENERICGEOMETRY_GEOMETRY_HH
