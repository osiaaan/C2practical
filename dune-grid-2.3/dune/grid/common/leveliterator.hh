// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRID_LEVELITERATOR_HH
#define DUNE_GRID_LEVELITERATOR_HH

#include <cstddef>
#include <iterator>

#include <dune/grid/common/entityiterator.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune
{

  /**********************************************************************/
  /** @brief Enables iteration over all entities
          of a given codimension and level of a grid.
          See also the documentation of Dune::EntityPointer.

      \note The LevelIterator interface is deprecated. Use the EntityIterator
            interface instead.

     @ingroup GIEntityPointer
   */
  template<int codim, PartitionIteratorType pitype, class GridImp,
      template<int,PartitionIteratorType,class> class LevelIteratorImp>
  class LevelIterator
    : public EntityIterator< codim, GridImp, LevelIteratorImp< codim, pitype, GridImp > >
  {
    typedef EntityIterator< codim, GridImp, LevelIteratorImp< codim, pitype, GridImp > > Base;

  public:
    /**
       @brief Preincrement operator.

       @note Forwarded to LevelIteratorImp.increment()
     */
    LevelIterator& operator++()
    {
      ++static_cast< Base & >( *this );
      return *this;
    }

  };

}

namespace std {

  template
  < int codim, Dune::PartitionIteratorType pitype, class GridImp,
      template<int,Dune::PartitionIteratorType,class> class LevelIteratorImp>
  struct iterator_traits<Dune::LevelIterator<codim, pitype, GridImp,
          LevelIteratorImp> > {
    typedef ptrdiff_t difference_type;
    typedef const typename Dune::LevelIterator<codim, pitype, GridImp,
        LevelIteratorImp>::Entity value_type;
    typedef value_type* pointer;
    typedef value_type& reference;
    typedef forward_iterator_tag iterator_category;
  };

} // namespace std

#endif // DUNE_GRID_LEVELITERATOR_HH
