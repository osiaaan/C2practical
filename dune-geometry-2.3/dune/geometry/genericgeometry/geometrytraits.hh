// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_GEOMETRYTRAITS_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_GEOMETRYTRAITS_HH

#include "../type.hh"
#include "matrixhelper.hh"
#include "cornermapping.hh"

namespace Dune
{
  namespace GenericGeometry
  {

    // DuneCoordTraits
    // ---------------

    template< class ct >
    struct DuneCoordTraits
    {
      typedef ct ctype;

      template< int dim >
      struct Vector
      {
        typedef FieldVector< ctype, dim > type;
      };

      template< int rows, int cols >
      struct Matrix
      {
        typedef FieldMatrix< ctype, rows, cols > type;
      };

      // This limit is, e.g., used in the termination criterion of the Newton
      // scheme within the generic implementation of the method local
      static const ctype epsilon ()
      {
        return 1e-6;
      }
    };



    // MappingTraits
    // -------------
    /** \ingroup GenericGeometry
     *  \brief Default mapping traits using Dune::FieldVector and
     *  Dune::FieldMatrix
     */
    template< class CT, unsigned int dim, unsigned int dimW >
    struct MappingTraits
    {
      typedef CT CoordTraits;

      static const unsigned int dimension = dim;
      static const unsigned int dimWorld = dimW;

      typedef typename CoordTraits :: ctype FieldType;
      typedef typename CoordTraits :: template Vector< dimension > :: type LocalCoordinate;
      typedef typename CoordTraits :: template Vector< dimWorld > :: type GlobalCoordinate;

      typedef typename CoordTraits :: template Matrix< dimWorld, dimension > :: type
      JacobianType;
      typedef typename CoordTraits :: template Matrix< dimension, dimWorld > :: type
      JacobianTransposedType;

      typedef GenericGeometry :: MatrixHelper< CoordTraits > MatrixHelper;
    };



    /** \brief If not affine only volume is cached (based on intElCompute)
     * otherwise all quantities can be cached
     */
    enum EvaluationType
    {
      //! assign if method called using barycenter
      ComputeOnDemand,
      //! assign in constructor using barycenter
      PreCompute
    };



    // DefaultGeometryTraits
    // ---------------------

    /** \ingroup GenericGeometry
     *  \brief   default settings for BasicGeometry
     *
     *  The class BasicGeometry requires a template argument <em>Traits</em>.
     *  These traits specify which reference mapping shall be used by the
     *  geometry and tweaks some performance settings.
     *
     *  This default implementation serves two purposed. Firstly, it documents
     *  the expected parameters. Secondly, the user of BasicGeometry can
     *  derive his traits class from DefaultGeometryTraits. Then, only the
     *  non-default settings have to be specified. Moreover, deriving from
     *  DefaultGeometryTraits makes the user code more robust to changes in
     *  the generic geometries.
     *
     *  \note DefaultGeometryTraits can directly be used for the
     *        <em>Traits</em> argument of BasicGeometry.
     *
     * \tparam ctype Type used for coordinate coefficients
     * \tparam dimG Dimension of the grid (note: This template parameter exists
     * for backward-compatibility only.  It is not used anywhere.)
     * \tparam dimW Dimension of the range space of this geometry
     * \tparam alwaysAffine Set to true if geometry is always affine (enables a few optimizations)
     */
    template< class ctype, int dimG, int dimW, bool alwaysAffine = false >
    struct DefaultGeometryTraits
    {
      //! types needed in matrix-vector operations
      typedef DuneCoordTraits< ctype > CoordTraits;

      //! dimension of the world
      static const int dimWorld = dimW;

      /** \brief will there be only one geometry type for a dimension?
       *
       *  If multiple geometry types are requested for a dimension, all methods
       *  of the geometry implementation are virtual (but no other branching for
       *  the geometry type is used).
       *
       *  If there is only a single geometry type for a certain dimension,
       *  <em>hasSingleGeometryType::v</em> can be set to true.
       *  In this case, virtual methods are not necessary and the geometries
       *  are a little faster.
       *
       *  If <em>hasSingleGeometryType::v</em> is set to true, an additional
       *  parameter <em>topologyId</em> is required.
       *  Here's an example:
       *  \code
       *  static const unsigned int topologyId = SimplexTopology< dim >::type::id;
       *  \endcode
       */
      template< int dim >
      struct hasSingleGeometryType
      {
        static const bool v = false;
        static const unsigned int topologyId = ~0u;
      };

      /** \brief specifies the reference mapping to be used
       *
       *  \tparam  Topology  type of topology for which the mapping
       *                     implementation is specified
       *
       *  This structure contains a single tydedef <em>type</em> specifying
       *  the implementation of the reference mapping. Basically, it looks like
       *  \code
       *  typedef CornerMapping< ... > type;
       *  \endcode
       */
      template< class Topology >
      struct Mapping
      {
        typedef CoordStorage< CoordTraits, Topology, dimWorld > CornerStorage;
        typedef CornerMapping< CoordTraits, Topology, dimWorld, CornerStorage, alwaysAffine > type;
      };

      /** \brief specifies how constant values are to be cached
       *
       *  This structure contains 3 parameters of type
       *  GenericGeometry::EvaluationType:
       *  - evaluateJacobianTransposed
       *  - evaluateJacobianInverseTransposed
       *  - evaluateIntegrationElement
       *  .
       *  These parameters control how eagerly these evaluations shall be
       *  performed in the case of an affine mapping.
       */
      struct Caching
      {
        static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
        static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
        static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
      };

      /** \brief type of additional user data to be stored in each mapping
       *
       *  Each HybridMapping and NonHybridMapping stores a user data structure
       *  of this type.
       */
      struct UserData {};
    };



    /** \struct  GlobalGeometryTraits
     *  \ingroup GenericGeometry
     *  \brief   grid specific information required by GenericGeometry::Geometry
     *
     *  Every implementation of a DUNE Geometry is required to have the same
     *  template parameter list:
     *  \code
     *  template< int mydim, int cdim, class Grid >
     *  \endcode
     *  Consequently, there is no direct way to pass compile time static
     *  information to a unified implementation such as the generic geometries.
     *  The structure GeometryTraits realizes an indirect way to do this.
     *
     *  For every grid implementation using the generic geometries, this
     *  structure must be specialized. The following default implementation
     *  can be used (via subclassing) to provide the necessary information. It
     *  contains exactly the fields that are necessary:
     *  \code
     *  template< class ctype, int dimW >
     *  struct DefaultGeometryTraits
     *  {
     *    typedef DuneCoordTraits< ctype > CoordTraits;
     *
     *    static const int dimWorld = dimW;
     *
     *    //   hybrid   [ true if Codim 0 is hybrid ]
     *    static const bool hybrid = true;
     *
     *    template< class Topology >
     *    struct Mapping
     *    {
     *      typedef MappingTraits< CoordTraits, Topology :: dimension, dimWorld > Traits;
     *      typedef CoordPointerStorage< Topology, typename Traits :: GlobalCoordinate >
     *        CornerStorage;
     *      typedef CornerMapping< Topology, Traits, CornerStorage > type;
     *    };
     *
     *    struct Caching
     *    {
     *      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
     *    };
     *  };
     *  \endcode
     *
     *  This implementation specifies the information used by
     *  GenericGeometry::Geometry.
     *
     *  \tparam  Grid  type of the grid, this traits class applies to
     */
    template< class Grid >
    struct GlobalGeometryTraits;

    template< class Grid >
    struct GlobalGeometryTraits< const Grid >
      : public GlobalGeometryTraits< Grid >
    {};



    /** \struct  LocalGeometryTraits
     *  \ingroup GenericGeometry
     *  \brief   grid specific information required by GenericGeometry::LocalGeometry
     *
     *  Every implementation of a DUNE Geometry is required to have the same
     *  template parameter list:
     *  \code
     *  template< int mydim, int cdim, class Grid >
     *  \endcode
     *  Consequently, there is no direct way to pass compile time static
     *  information to a unified implementation such as the generic geometries.
     *  The structure GeometryTraits realizes an indirect way to do this.
     *
     *  For every grid implementation using the generic geometries, this
     *  structure must be specialized. The following default implementation
     *  can be used (via subclassing) to provide the necessary information. It
     *  contains exactly the fields that are necessary:
     *  \code
     *  template< class ctype, int dimW >
     *  struct DefaultGeometryTraits
     *  {
     *    typedef DuneCoordTraits< ctype > CoordTraits;
     *
     *    static const int dimWorld = dimW;
     *
     *    //   hybrid   [ true if Codim 0 is hybrid ]
     *    static const bool hybrid = true;
     *
     *    template< class Topology >
     *    struct Mapping
     *    {
     *      typedef MappingTraits< CoordTraits, Topology :: dimension, dimWorld > Traits;
     *      typedef CoordPointerStorage< Topology, typename Traits :: GlobalCoordinate >
     *        CornerStorage;
     *      typedef CornerMapping< Topology, Traits, CornerStorage > type;
     *    };
     *
     *    struct Caching
     *    {
     *      static const EvaluationType evaluateJacobianTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateJacobianInverseTransposed = ComputeOnDemand;
     *      static const EvaluationType evaluateIntegrationElement = ComputeOnDemand;
     *    };
     *  };
     *  \endcode
     *
     *  This implementation specifies the information used by
     *  GenericGeometry::LocalGeometry.
     *
     *  \tparam  Grid  type of the grid, this traits class applies to
     */
    template< class Grid >
    struct LocalGeometryTraits;

    template< class Grid >
    struct LocalGeometryTraits< const Grid >
      : public LocalGeometryTraits< Grid >
    {};
  }

}

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_GEOMETRYTRAITS_HH
