// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_GENERICGEOMETRY_CORNERMAPPING_HH
#define DUNE_GEOMETRY_GENERICGEOMETRY_CORNERMAPPING_HH

#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/genericgeometry/referenceelements.hh>
#include <dune/geometry/genericgeometry/matrixhelper.hh>

namespace Dune
{

  namespace GenericGeometry
  {

    // External Forward Declarations
    // -----------------------------

    template< class CT, unsigned int dim, unsigned int dimW >
    struct MappingTraits;



    // GenericCornerMapping
    // --------------------

    template< class Topology, class Traits, bool affine, unsigned int offset = 0 >
    class GenericCornerMapping;

    template< class Traits, bool affine, unsigned int offset >
    class GenericCornerMapping < Point, Traits, affine, offset >
    {
      typedef Point Topology;

    public:
      static const unsigned int dim = Topology::dimension;
      static const unsigned int dimW = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;
      typedef typename Traits::JacobianTransposedType JacobianTransposedType;

      static const bool alwaysAffine = true;

      template< class CoordStorage >
      static const GlobalCoordinate &origin ( const CoordStorage &coords )
      {
        dune_static_assert( CoordStorage::size, "Invalid offset." );
        return coords[ offset ];
      }

      template< class CoordStorage >
      static void phi_set ( const CoordStorage &coords,
                            const LocalCoordinate &x,
                            const FieldType &factor,
                            GlobalCoordinate &p )
      {
        const GlobalCoordinate &y = origin( coords );
        for( unsigned int i = 0; i < dimW; ++i )
          p[ i ] = factor * y[ i ];
      }

      template< class CoordStorage >
      static void phi_add ( const CoordStorage &coords,
                            const LocalCoordinate &x,
                            const FieldType &factor,
                            GlobalCoordinate &p )
      {
        const GlobalCoordinate &y = origin( coords );
        for( unsigned int i = 0; i < dimW; ++i )
          p[ i ] += factor * y[ i ];
      }

      template< class CoordStorage >
      static bool Dphi_set ( const CoordStorage &coords,
                             const LocalCoordinate &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        return true;
      }

      template< class CoordStorage >
      static bool Dphi_add ( const CoordStorage &coords,
                             const LocalCoordinate &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        return true;
      }
    };


    template< class BaseTopology, class Traits, bool affine, unsigned int offset >
    class GenericCornerMapping< Prism< BaseTopology >, Traits, affine, offset >
    {
      typedef Prism< BaseTopology > Topology;

      typedef GenericCornerMapping< BaseTopology, Traits, affine, offset >
      BottomMapping;
      typedef GenericCornerMapping
      < BaseTopology, Traits, affine, offset + BaseTopology::numCorners >
      TopMapping;

    public:
      static const unsigned int dim = Topology::dimension;
      static const unsigned int dimW = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;
      typedef typename Traits::JacobianTransposedType JacobianTransposedType;

      static const bool alwaysAffine = ((dim < 2) || affine);

      template< class CoordStorage >
      static const GlobalCoordinate &origin ( const CoordStorage &coords )
      {
        return BottomMapping::origin( coords );
      }

      template< class CoordStorage >
      static void phi_set ( const CoordStorage &coords,
                            const LocalCoordinate &x,
                            const FieldType &factor,
                            GlobalCoordinate &p )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = FieldType( 1 ) - xn;
        BottomMapping::phi_set( coords, x, factor * cxn, p );
        TopMapping::phi_add( coords, x, factor * xn, p );
      }

      template< class CoordStorage >
      static void phi_add ( const CoordStorage &coords,
                            const LocalCoordinate &x,
                            const FieldType &factor,
                            GlobalCoordinate &p )
      {
        const FieldType xn = x[ dim-1 ];
        const FieldType cxn = FieldType( 1 ) - xn;
        BottomMapping::phi_add( coords, x, factor * cxn, p );
        TopMapping::phi_add( coords, x, factor * xn, p );
      }

      template< class CoordStorage >
      static bool Dphi_set ( const CoordStorage &coords,
                             const LocalCoordinate &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        const FieldType xn = x[ dim-1 ];
        bool isAffine = true;
        if( alwaysAffine )
        {
          const FieldType cxn = FieldType( 1 ) - xn;
          BottomMapping::Dphi_set( coords, x, factor * cxn, J );
          TopMapping::Dphi_add( coords, x, factor * xn, J );
        }
        else
        {
          JacobianTransposedType Jtop;
          isAffine &= BottomMapping::Dphi_set( coords, x, factor, J );
          isAffine &= TopMapping::Dphi_set( coords, x, factor, Jtop );

          FieldType norm = FieldType( 0 );
          for( unsigned int i = 0; i < dim-1; ++i )
          {
            Jtop[ i ] -= J[ i ];
            norm += Jtop[ i ].two_norm2();
            J[ i ].axpy( xn, Jtop[ i ] );
          }
          isAffine &= (norm < 1e-12);
        }
        BottomMapping::phi_set( coords, x, -factor, J[ dim-1 ] );
        TopMapping::phi_add( coords, x, factor, J[ dim-1 ] );
        return isAffine;
      }

      template< class CoordStorage >
      static bool Dphi_add ( const CoordStorage &coords,
                             const LocalCoordinate &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        const FieldType xn = x[ dim-1 ];
        bool isAffine = true;
        if( alwaysAffine )
        {
          const FieldType cxn = FieldType( 1 ) - xn;
          BottomMapping::Dphi_add( coords, x, factor * cxn, J );
          TopMapping::Dphi_add( coords, x, factor * xn, J );
        }
        else
        {
          JacobianTransposedType Jbottom, Jtop;
          isAffine &= BottomMapping::Dphi_set( coords, x, FieldType( 1 ), Jbottom );
          isAffine &= TopMapping::Dphi_set( coords, x, FieldType( 1 ), Jtop );

          FieldType norm = FieldType( 0 );
          for( unsigned int i = 0; i < dim-1; ++i )
          {
            Jtop[ i ] -= Jbottom[ i ];
            norm += Jtop[ i ].two_norm2();
            J[ i ].axpy( factor, Jbottom[ i ] );
            J[ i ].axpy( factor*xn, Jtop[ i ] );
          }
          isAffine &= (norm < 1e-12);
        }
        BottomMapping::phi_add( coords, x, -factor, J[ dim-1 ] );
        TopMapping::phi_add( coords, x, factor, J[ dim-1 ] );
        return isAffine;
      }
    };


    template< class BaseTopology, class Traits, bool affine, unsigned int offset >
    class GenericCornerMapping < Pyramid< BaseTopology >, Traits, affine, offset >
    {
      typedef Pyramid< BaseTopology > Topology;

      typedef GenericCornerMapping< BaseTopology, Traits, affine, offset >
      BottomMapping;
      typedef GenericCornerMapping
      < Point, Traits, affine, offset + BaseTopology::numCorners >
      TopMapping;

    public:
      static const unsigned int dim = Topology::dimension;
      static const unsigned int dimW = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;
      typedef typename Traits::JacobianTransposedType JacobianTransposedType;

      static const bool alwaysAffine = (BottomMapping::alwaysAffine || affine);

      template< class CoordStorage >
      static const GlobalCoordinate &origin ( const CoordStorage &coords )
      {
        return BottomMapping::origin( coords );
      }

      template< class CoordStorage >
      static void phi_set ( const CoordStorage &coords,
                            const LocalCoordinate &x,
                            const FieldType &factor,
                            GlobalCoordinate &p )
      {
        const FieldType xn = x[ dim-1 ];
        if( alwaysAffine )
        {
          const GlobalCoordinate &top = TopMapping::origin( coords );
          const GlobalCoordinate &bottom = BottomMapping::origin( coords );

          BottomMapping::phi_set( coords, x, factor, p );
          for( unsigned int i = 0; i < dimW; ++i )
            p[ i ] += (factor * xn) * (top[ i ] - bottom[ i ]);
        }
        else
        {
          TopMapping::phi_set( coords, x, factor * xn, p );
          const FieldType cxn = FieldType( 1 ) - xn;
          if( cxn > 1e-12 )
          {
            const FieldType icxn = FieldType( 1 ) / cxn;
            LocalCoordinate xb;
            for( unsigned int i = 0; i < dim-1; ++i )
              xb[ i ] = icxn * x[ i ];

            BottomMapping::phi_add( coords, xb, factor * cxn, p );
          }
        }
      }

      template< class CoordStorage >
      static void phi_add ( const CoordStorage &coords,
                            const LocalCoordinate &x,
                            const FieldType &factor,
                            GlobalCoordinate &p )
      {
        const FieldType xn = x[ dim-1 ];
        if( alwaysAffine )
        {
          const GlobalCoordinate &top = TopMapping::origin( coords );
          const GlobalCoordinate &bottom = BottomMapping::origin( coords );

          BottomMapping::phi_add( coords, x, factor, p );
          for( unsigned int i = 0; i < dimW; ++i )
            p[ i ] += (factor * xn) * (top[ i ] - bottom[ i ]);
        }
        else
        {
          TopMapping::phi_add( coords, x, factor * xn, p );
          const FieldType cxn = FieldType( 1 ) - xn;
          if( cxn > 1e-12 )
          {
            const FieldType icxn = FieldType( 1 ) / cxn;
            LocalCoordinate xb;
            for( unsigned int i = 0; i < dim-1; ++i )
              xb[ i ] = icxn * x[ i ];

            BottomMapping::phi_add( coords, xb, factor * cxn, p );
          }
        }
      }

      template< class CoordStorage >
      static bool Dphi_set ( const CoordStorage &coords,
                             const LocalCoordinate &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        GlobalCoordinate &q = J[ dim-1 ];
        bool isAffine;
        if( alwaysAffine )
        {
          const GlobalCoordinate &top = TopMapping::origin( coords );
          const GlobalCoordinate &bottom = BottomMapping::origin( coords );

          isAffine = BottomMapping::Dphi_set( coords, x, factor, J );
          for( unsigned int i = 0; i < dimW; ++i )
            q[ i ] = factor * (top[ i ] - bottom[ i ]);
        }
        else
        {
          const FieldType xn = x[ dim-1 ];
          const FieldType cxn = FieldType( 1 ) - xn;
          const FieldType icxn = FieldType( 1 ) / cxn;
          LocalCoordinate xb;
          for( unsigned int i = 0; i < dim-1; ++i )
            xb[ i ] = icxn * x[ i ];
          isAffine = BottomMapping::Dphi_set( coords, xb, factor, J );

          TopMapping::phi_set( coords, x, factor, q );
          BottomMapping::phi_add( coords, xb, -factor, q );
          xb *= factor;
          for( unsigned int j = 0; j < dim-1; ++j )
          {
            for( unsigned int i = 0; i < dimW; ++i )
              q[ i ] += J[ j ][ i ] * xb[ j ];
          }
        }
        return isAffine;
      }

      template< class CoordStorage >
      static bool Dphi_add ( const CoordStorage &coords,
                             const LocalCoordinate &x,
                             const FieldType &factor,
                             JacobianTransposedType &J )
      {
        GlobalCoordinate &q = J[ dim-1 ];
        bool isAffine;
        if( alwaysAffine )
        {
          const GlobalCoordinate &top = TopMapping::origin( coords );
          const GlobalCoordinate &bottom = BottomMapping::origin( coords );

          isAffine = BottomMapping::Dphi_add( coords, x, factor, J );
          for( unsigned int i = 0; i < dimW; ++i )
            q[ i ] += factor * (top[ i ] - bottom[ i ]);
        }
        else
        {
          const FieldType xn = x[ dim-1 ];
          const FieldType cxn = FieldType( 1 ) - xn;
          const FieldType icxn = FieldType( 1 ) / cxn;
          LocalCoordinate xb;
          for( unsigned int i = 0; i < dim-1; ++i )
            xb[ i ] = icxn * x[ i ];
          isAffine = BottomMapping::Dphi_add( coords, xb, factor, J );

          TopMapping::phi_add( coords, x, factor, q );
          BottomMapping::phi_add( coords, xb, -factor, q );
          xb *= factor;
          for( unsigned int j = 0; j < dim-1; ++j )
          {
            for( unsigned int i = 0; i < dimW; ++i )
              q[ i ] += J[ j ][ i ] * xb[ j ];
          }
        }
        return isAffine;
      }
    };



    // SubMappingCoords
    // ----------------

    template< class Mapping, unsigned int codim >
    class SubMappingCoords
    {
      typedef typename Mapping::GlobalCoordinate GlobalCoordinate;
      typedef typename Mapping::ReferenceElement ReferenceElement;

      static const unsigned int dimension = ReferenceElement::dimension;

      const Mapping &mapping_;
      const unsigned int i_;

    public:
      SubMappingCoords ( const Mapping &mapping, unsigned int i )
        : mapping_( mapping ), i_( i )
      {}

      const GlobalCoordinate &operator[] ( unsigned int j ) const
      {
        const unsigned int k
          = ReferenceElement::template subNumbering< codim, dimension - codim >( i_, j );
        return mapping_.corner( k );
      }
    };



    // CoordStorage
    // ------------

    /** \class CoordStorage
     *  \ingroup GenericGeometry
     *  \brief
     */
    template< class CoordTraits, class Topology, unsigned int dimW >
    class CoordStorage
    {
      typedef CoordStorage< CoordTraits, Topology, dimW > This;

    public:
      static const unsigned int size = Topology::numCorners;

      static const unsigned int dimWorld = dimW;

      typedef typename CoordTraits::template Vector< dimWorld >::type GlobalCoordinate;

      template< class SubTopology >
      struct SubStorage
      {
        typedef CoordStorage< CoordTraits, SubTopology, dimWorld > type;
      };

      template< class CoordVector >
      explicit CoordStorage ( const CoordVector &coords )
      {
        for( unsigned int i = 0; i < size; ++i )
          coords_[ i ] = coords[ i ];
      }

      const GlobalCoordinate &operator[] ( unsigned int i ) const
      {
        return coords_[ i ];
      }

    private:
      GlobalCoordinate coords_[ size ];
    };



    // CoordPointerStorage
    // -------------------

    /** \class CoordPointerStorage
     *  \ingroup GenericGeometry
     *  \brief
     */
    template< class CoordTraits, class Topology, unsigned int dimW >
    class CoordPointerStorage
    {
      typedef CoordPointerStorage< CoordTraits, Topology, dimW > This;

    public:
      static const unsigned int size = Topology::numCorners;

      static const unsigned int dimWorld = dimW;

      typedef typename CoordTraits::template Vector< dimWorld >::type GlobalCoordinate;

      template< class SubTopology >
      struct SubStorage
      {
        typedef CoordPointerStorage< CoordTraits, SubTopology, dimWorld > type;
      };

      template< class CoordVector >
      explicit CoordPointerStorage ( const CoordVector &coords )
      {
        for( unsigned int i = 0; i < size; ++i )
          coords_[ i ] = &(coords[ i ]);
      }

      const GlobalCoordinate &operator[] ( unsigned int i ) const
      {
        return *(coords_[ i ]);
      }

    private:
      const GlobalCoordinate *coords_[ size ];
    };



    // CornerMapping
    // -------------

    /** \class CornerMapping
     *  \ingroup GenericGeometry
     *  \brief implementation of GenericGeometry::Mapping for first order
     *  lagrange type reference mappings.
     */
    template< class CoordTraits, class Topo, unsigned int dimW,
        class CStorage = CoordStorage< CoordTraits, Topo, dimW >,
        bool affine = false >
    class CornerMapping
    {
      typedef CornerMapping< CoordTraits, Topo, dimW, CStorage, affine > This;

    public:
      typedef Topo Topology;
      typedef CStorage CornerStorage;
      typedef MappingTraits< CoordTraits, Topology::dimension, dimW > Traits;

      static const unsigned int dimension = Traits::dimension;
      static const unsigned int dimWorld = Traits::dimWorld;

      typedef typename Traits::FieldType FieldType;
      typedef typename Traits::LocalCoordinate LocalCoordinate;
      typedef typename Traits::GlobalCoordinate GlobalCoordinate;
      typedef typename Traits::JacobianType JacobianType;
      typedef typename Traits::JacobianTransposedType JacobianTransposedType;

      typedef GenericGeometry::ReferenceElement< Topology, FieldType > ReferenceElement;

      template< unsigned int codim, unsigned int i >
      struct SubTopology
      {
        typedef typename GenericGeometry::SubTopology< Topo, codim, i >::type Topology;
        typedef typename CStorage::template SubStorage< Topology >::type CornerStorage;
        typedef CornerMapping< CoordTraits, Topology, dimWorld, CornerStorage, affine > Trace;
      };

    private:
      typedef GenericGeometry::GenericCornerMapping< Topology, Traits, affine > GenericMapping;

    public:
      static const bool alwaysAffine = GenericMapping::alwaysAffine;

      template< class CoordVector >
      explicit CornerMapping ( const CoordVector &coords )
        : coords_( coords )
      {}

      const GlobalCoordinate &corner ( int i ) const
      {
        return coords_[ i ];
      }

      void global ( const LocalCoordinate &x, GlobalCoordinate &y ) const
      {
        GenericMapping::phi_set( coords_, x, FieldType( 1 ), y );
      }

      bool jacobianTransposed ( const LocalCoordinate &x,
                                JacobianTransposedType &JT ) const
      {
        return GenericMapping::Dphi_set( coords_, x, FieldType( 1 ), JT );
      }

      template< unsigned int codim, unsigned int i >
      typename SubTopology< codim, i >::Trace trace () const
      {
        typedef typename SubTopology< codim, i >::Trace Trace;
        typedef SubMappingCoords< This, codim > CoordVector;
        return Trace( CoordVector( *this, i ) );
      }

    protected:
      CornerStorage coords_;
    };

  } // namespace GenericGeometry

} // namespace Dune

#endif // #ifndef DUNE_GEOMETRY_GENERICGEOMETRY_CORNERMAPPING_HH
