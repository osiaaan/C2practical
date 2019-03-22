/******************************************************************************

  ALUGrid - a library providing a mesh manager supporting simplicial
  and hexahedral meshes, local grid adaptivity for use in parallel
  computations including dynamic load balancing.

  Copyright (C) 1998 - 2002 Bernhard Schupp
  Copyright (C) 1998 - 2002 Mario Ohlberger
  Copyright (C) 2004 - 2012 Robert Kloefkorn
  Copyright (C) 2005 - 2012 Andreas Dedner
  Copyright (C) 2010 - 2012 Martin Nolte

  The DUNE ALUGrid module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The ALUGrid library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

******************************************************************************/

/******************************************************************************

  ALUGrid - a library providing a mesh manager supporting simplicial
  and hexahedral meshes, local grid adaptivity for use in parallel
  computations including dynamic load balancing.

  Copyright (C) 1998 - 2002 Bernhard Schupp
  Copyright (C) 1998 - 2002 Mario Ohlberger
  Copyright (C) 2004 - 2012 Robert Kloefkorn
  Copyright (C) 2005 - 2012 Andreas Dedner
  Copyright (C) 2010 - 2012 Martin Nolte

  The DUNE ALUGrid module is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License as
  published by the Free Software Foundation; either version 2 of
  the License, or (at your option) any later version.

  The ALUGrid library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

******************************************************************************/

// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_MACROGRIDVIEW_HH
#define DUNE_MACROGRIDVIEW_HH

#include <dune/common/typetraits.hh>
#include <dune/common/exceptions.hh>

#include <dune/grid/common/capabilities.hh>
#include <dune/grid/common/gridview.hh>

namespace Dune
{

  template< class GridImp, PartitionIteratorType pitype >
  class MacroGridView;

  template< class GridImp, PartitionIteratorType pitype >
  struct MacroGridViewTraits
  {
    typedef MacroGridView< GridImp, pitype > GridViewImp;

    /** \brief type of the grid */
    typedef typename remove_const<GridImp>::type Grid;

    /** \brief type of the index set */
    typedef typename Grid :: Traits :: LevelIndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Grid :: Traits :: LevelIntersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Grid :: Traits :: LevelIntersectionIterator
      IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Grid :: Traits :: CollectiveCommunication CollectiveCommunication;

    template< int cd >
    struct Codim
    {
      typedef typename Grid :: Traits
        :: template Codim< cd > :: template Partition< pitype > :: LevelIterator
        Iterator;

      typedef typename Grid :: Traits :: template Codim< cd > :: Entity Entity;
      typedef typename Grid :: Traits :: template Codim< cd > :: EntityPointer
        EntityPointer;

      typedef typename Grid :: template Codim< cd > :: Geometry Geometry;
      typedef typename Grid :: template Codim< cd > :: LocalGeometry
        LocalGeometry;

      /** \brief Define types needed to iterate over entities of a given partition type */
      template< PartitionIteratorType pit >
      struct Partition
      {
        /** \brief iterator over a given codim and partition type */
        typedef typename Grid :: template Codim< cd >
          :: template Partition< pit > :: LevelIterator
          Iterator;
      };
    };

    enum { conforming = Capabilities :: isLevelwiseConforming< Grid > :: v };
  };


  template< class GridImp, PartitionIteratorType pitype >
  class MacroGridView 
  {
    typedef MacroGridView< GridImp, pitype > ThisType;

  public:
    typedef MacroGridViewTraits< GridImp, pitype > Traits;

    /** \brief type of the grid */
    typedef typename Traits::Grid Grid;

    /** \brief type of the index set */
    typedef typename Traits :: IndexSet IndexSet;

    /** \brief type of the intersection */
    typedef typename Traits :: Intersection Intersection;

    /** \brief type of the intersection iterator */
    typedef typename Traits :: IntersectionIterator IntersectionIterator;

    /** \brief type of the collective communication */
    typedef typename Traits :: CollectiveCommunication CollectiveCommunication;

    /** \brief Codim Structure */
    template< int cd >
    struct Codim : public Traits :: template Codim<cd> {};
 
    enum { conforming = Traits :: conforming };

    MacroGridView ( const Grid &grid )
    : grid_( &grid ),
      level_( 0 )
    {}

    /** \brief obtain a const reference to the underlying hierarchic grid */
    const Grid &grid () const
    {
      assert( grid_ );
      return *grid_;
    }

    /** \brief obtain the index set */
    const IndexSet &indexSet () const
    {
      return grid().levelIndexSet( level_ );
    }
    
    /** \brief obtain number of entities in a given codimension */
    int size ( int codim ) const
    {
      return grid().size( level_, codim );
    }

    /** \brief obtain number of entities with a given geometry type */
    int size ( const GeometryType &type ) const
    {
      return grid().size( level_, type );
    }

    /** \brief obtain begin iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator begin () const
    {
      return grid().template lbegin< cd, pitype >( level_ );
    }

    /** \brief obtain begin iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator begin () const
    {
      return grid().template lbegin< cd, pit >( level_ );
    }

    /** \brief obtain end iterator for this view */
    template< int cd >
    typename Codim< cd > :: Iterator end () const
    {
      return grid().template lend< cd, pitype >( level_ );
    }

    /** \brief obtain end iterator for this view */
    template< int cd, PartitionIteratorType pit >
    typename Codim< cd > :: template Partition< pit > :: Iterator end () const
    {
      return grid().template lend< cd, pit >( level_ );
    }

    /** \brief obtain begin intersection iterator with respect to this view */
    IntersectionIterator
    ibegin ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.impl().ilevelbegin();
    }

    /** \brief obtain end intersection iterator with respect to this view */
    IntersectionIterator
    iend ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.impl().ilevelend();
    }

    /** \brief obtain collective communication object */
    const CollectiveCommunication &comm () const
    {
      return grid().comm();
    }

    /** \brief Return size of the overlap region for a given codim on the grid view.  */
    int overlapSize(int codim) const
    {
      return grid().overlapSize(level_, codim);
    }

    /** \brief Return size of the ghost region for a given codim on the grid view.  */
    int ghostSize(int codim) const 
    {
      return grid().ghostSize(level_, codim);
    }

    /** communicate data on this view */
    template< class DataHandleImp, class DataType >
    void communicate ( CommDataHandleIF< DataHandleImp, DataType > &data,
                       InterfaceType iftype,
                       CommunicationDirection dir ) const
    {
      return grid().communicate( data, iftype, dir, level_ );
    }

    //** extra methods for load balancing */
    int
    master ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.impl().master();
    }
    int
    macroId ( const typename Codim< 0 > :: Entity &entity ) const
    {
      return entity.impl().macroId();
    }
    int
    weight ( const IntersectionIterator &intersectionIterator ) const // should perhaps be intersection but that class a default is used...
    {
      return intersectionIterator.impl().weight();
    }

  private:
    const Grid *grid_;
    const int level_;
  };
}

#endif
