// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GEOMETRY_REFINEMENT_PYRAMIDTRIANGULATION_CC
#define DUNE_GEOMETRY_REFINEMENT_PYRAMIDTRIANGULATION_CC

#include <dune/common/array.hh>
#include <dune/common/fvector.hh>
#include <dune/common/typetraits.hh>

#include <dune/geometry/genericgeometry/geometrytraits.hh>
#include <dune/geometry/genericgeometry/topologytypes.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/geometry/type.hh>

#include "base.cc"
#include "simplex.cc"

namespace Dune
{
  namespace RefinementImp
  {
    /*!
     * \brief This namespace contains the \ref Refinement
     * implementation for triangulating pyramids
     * (GeometryType::pyramid -> GeometryType::simplex)
     */
    namespace PyramidTriangulation
    {
      // ////////////
      //
      //  Utilities
      //

      using Simplex::getPermutation;
      using Simplex::referenceToKuhn;

      // ////////////////////////////////////
      //
      //  Refine a pyramid with simplices
      //

      // forward declaration of the iterator base
      template<int dimension, class CoordType, int codimension>
      class RefinementIteratorSpecial;

      /*
       * The permutations 0 and 1 of the Kuhn-decomposition of a cube into simplices form a pyramid.
       * The resulting pyramid is not oriented the same as the reference pyramid and so the Kuhn-coordinates
       * have to be transformed using the method below.
       */
      template<int dimension, class CoordType> FieldVector<CoordType, dimension>
      transformCoordinate(/*! Point to transform */ FieldVector<CoordType, dimension> point)
      {
        FieldVector<CoordType, dimension> transform;
        transform[0]=1-point[0];
        transform[1]=1-point[1];
        transform[2]=point[2];
        return transform;
      }

      /*!
       * \brief Implementation of the refinement of a pyramid into simplices.
       *
       * Note that the virtual vertices of two intersecting simplices might have copies, i.e.
       * by running over all vertices using the VertexIterator you might run over some twice.
       */
      template<int dimension_, class CoordType>
      class RefinementImp
      {
      public:
        enum {dimension = dimension_};

        typedef CoordType ctype;
        enum {dimensionworld = dimension};

        template<int codimension>
        struct Codim;
        typedef typename Codim<dimension>::SubEntityIterator VertexIterator;
        typedef FieldVector<CoordType, dimension> CoordVector;
        typedef typename Codim<0>::SubEntityIterator ElementIterator;
        typedef FieldVector<int, dimension+1> IndexVector;

        static int nVertices(int level);
        static VertexIterator vBegin(int level);
        static VertexIterator vEnd(int level);

        static int nElements(int level);
        static ElementIterator eBegin(int level);
        static ElementIterator eEnd(int level);

      private:
        friend class RefinementIteratorSpecial<dimension, CoordType, 0>;
        friend class RefinementIteratorSpecial<dimension, CoordType, dimension>;

        typedef Simplex::RefinementImp<dimension, CoordType> BackendRefinement;

        enum { nKuhnSimplices = 2 };
      };

      template<int dimension, class CoordType>
      template<int codimension>
      struct RefinementImp<dimension, CoordType>::Codim
      {
        class SubEntityIterator;
        typedef Dune::MultiLinearGeometry<CoordType,dimension-codimension,dimension> Geometry;
      };

      template<int dimension, class CoordType>
      int
      RefinementImp<dimension, CoordType>::
      nVertices(int level)
      {
        return BackendRefinement::nVertices(level) * nKuhnSimplices;
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vBegin(int level)
      {
        return VertexIterator(level);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::VertexIterator
      RefinementImp<dimension, CoordType>::
      vEnd(int level)
      {
        return VertexIterator(level, true);
      }

      template<int dimension, class CoordType>
      int
      RefinementImp<dimension, CoordType>::
      nElements(int level)
      {
        return BackendRefinement::nElements(level) * nKuhnSimplices;
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eBegin(int level)
      {
        return ElementIterator(level);
      }

      template<int dimension, class CoordType>
      typename RefinementImp<dimension, CoordType>::ElementIterator
      RefinementImp<dimension, CoordType>::
      eEnd(int level)
      {
        return ElementIterator(level, true);
      }

      // //////////////
      //
      // The iterator
      //

      // vertices
      template<int dimension, class CoordType>
      class RefinementIteratorSpecial<dimension, CoordType, dimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::CoordVector CoordVector;
        typedef typename Refinement::template Codim<dimension>::Geometry Geometry;

        RefinementIteratorSpecial(int level, bool end = false);

        void increment();

        CoordVector coords() const;

        Geometry geometry() const;

        int index() const;
      protected:
        typedef typename Refinement::BackendRefinement BackendRefinement;
        typedef typename BackendRefinement::template Codim<dimension>::SubEntityIterator BackendIterator;
        enum { nKuhnSimplices = 2 };

        int level;

        int kuhnIndex;
        BackendIterator backend;
        const BackendIterator backendEnd;
      };

      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      RefinementIteratorSpecial(int level_, bool end)
        : level(level_), kuhnIndex(0),
          backend(BackendRefinement::vBegin(level)),
          backendEnd(BackendRefinement::vEnd(level))
      {
        if (end)
          kuhnIndex = nKuhnSimplices;
      }

      template<int dimension, class CoordType>
      void
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      increment()
      {
        ++backend;
        if(backend == backendEnd)
        {
          backend = BackendRefinement::vBegin(level);
          ++kuhnIndex;
        }
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, dimension>::CoordVector
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      coords() const
      {
        return transformCoordinate(referenceToKuhn(backend.coords(),
                                                   getPermutation<dimension>(kuhnIndex)));
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, dimension>::Geometry
      RefinementIteratorSpecial<dimension, CoordType, dimension>::geometry () const
      {
        std::vector<CoordVector> corners(1);
        corners[0] = referenceToKuhn(backend.coords(), getPermutation<dimension>(kuhnIndex));
        return Geometry(GeometryType(0), corners);
      }

      template<int dimension, class CoordType>
      int
      RefinementIteratorSpecial<dimension, CoordType, dimension>::
      index() const
      {
        return kuhnIndex*BackendRefinement::nVertices(level) + backend.index();
      }

      // elements
      template<int dimension, class CoordType>
      class RefinementIteratorSpecial<dimension, CoordType, 0>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef typename Refinement::IndexVector IndexVector;
        typedef typename Refinement::CoordVector CoordVector;
        typedef typename Refinement::template Codim<0>::Geometry Geometry;

        RefinementIteratorSpecial(int level, bool end = false);

        void increment();

        IndexVector vertexIndices() const;
        int index() const;
        CoordVector coords() const;

        Geometry geometry() const;

      private:
        CoordVector global(const CoordVector &local) const;

      protected:
        typedef typename Refinement::BackendRefinement BackendRefinement;
        typedef typename BackendRefinement::template Codim<0>::SubEntityIterator BackendIterator;
        enum { nKuhnSimplices = 2};

        int level;

        int kuhnIndex;
        BackendIterator backend;
        const BackendIterator backendEnd;
      };

      template<int dimension, class CoordType>
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      RefinementIteratorSpecial(int level_, bool end)
        : level(level_), kuhnIndex(0),
          backend(BackendRefinement::eBegin(level)),
          backendEnd(BackendRefinement::eEnd(level))
      {
        if (end)
          kuhnIndex = nKuhnSimplices;
      }

      template<int dimension, class CoordType>
      void
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      increment()
      {
        ++backend;
        if (backend == backendEnd)
        {
          backend = BackendRefinement::eBegin(level);
          ++kuhnIndex;
        }
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::IndexVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      vertexIndices() const
      {
        IndexVector indices = backend.vertexIndices();

        int base = kuhnIndex * BackendRefinement::nVertices(level);
        indices += base;

        return indices;
      }

      template<int dimension, class CoordType>
      int
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      index() const
      {
        return kuhnIndex*BackendRefinement::nElements(level) + backend.index();
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::CoordVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      coords() const
      {
        return global(backend.coords());
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::Geometry
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      geometry() const
      {
        const typename BackendIterator::Geometry &
        bgeo = backend.geometry();
        std::vector<CoordVector> corners(dimension+1);
        for(int i = 0; i <= dimension; ++i)
          corners[i] = global(bgeo.corner(i));

        return Geometry(bgeo.type(), corners);
      }

      template<int dimension, class CoordType>
      typename RefinementIteratorSpecial<dimension, CoordType, 0>::
      CoordVector
      RefinementIteratorSpecial<dimension, CoordType, 0>::
      global(const CoordVector &local) const
      {
        return transformCoordinate(referenceToKuhn(local,
                                                   getPermutation<dimension>(kuhnIndex)));
      }

      // common
      template<int dimension, class CoordType>
      template<int codimension>
      class RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator
        : public ForwardIteratorFacade<typename RefinementImp<dimension, CoordType>::template Codim<codimension>::SubEntityIterator, int>,
          public RefinementIteratorSpecial<dimension, CoordType, codimension>
      {
      public:
        typedef RefinementImp<dimension, CoordType> Refinement;
        typedef SubEntityIterator This;

        SubEntityIterator(int level, bool end = false);

        bool equals(const This &other) const;
      protected:
        using RefinementIteratorSpecial<dimension, CoordType, codimension>::kuhnIndex;
        using RefinementIteratorSpecial<dimension, CoordType, codimension>::backend;
      };

#ifndef DOXYGEN
      template<int dimension, class CoordType>
      template<int codimension>
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      SubEntityIterator(int level, bool end)
        : RefinementIteratorSpecial<dimension, CoordType, codimension>(level, end)
      {}

      template<int dimension, class CoordType>
      template<int codimension>
      bool
      RefinementImp<dimension, CoordType>::Codim<codimension>::SubEntityIterator::
      equals(const This &other) const
      {
        return kuhnIndex == other.kuhnIndex && backend == other.backend;
      }
#endif

    } // namespace PyramidTriangulation
  } // namespace RefinementImp

  namespace RefinementImp
  {
    // ///////////////////////
    //
    // The refinement traits
    //
#ifndef DOXYGEN
    template<unsigned topologyId, class CoordType, unsigned coerceToId>
    struct Traits<
        topologyId, CoordType, coerceToId, 3,
        typename enable_if<
            (GenericGeometry::PyramidTopology<3>::type::id >> 1) ==
            (topologyId >> 1) &&
            (GenericGeometry::SimplexTopology<3>::type::id >> 1) ==
            (coerceToId >> 1)
            >::type>
    {
      typedef PyramidTriangulation::RefinementImp<3, CoordType> Imp;
    };
#endif

  } // namespace RefinementImp
} // namespace Dune

#endif // DUNE_GEOMETRY_REFINEMENT_PYRAMIDTRIANGULATION_CC
