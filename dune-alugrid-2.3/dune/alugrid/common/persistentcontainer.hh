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

#ifndef DUNE_ALU_PERSISTENTCONTAINER_HH
#define DUNE_ALU_PERSISTENTCONTAINER_HH

#include <dune/grid/utility/persistentcontainer.hh>
#include <dune/grid/utility/persistentcontainervector.hh>

#include <dune/alugrid/grid.hh>

namespace Dune
{

  // ALUGridPersistentContainer
  // --------------------------

  template< class G, class T >
  class ALUGridPersistentContainer
  : public PersistentContainerVector< G, typename G::HierarchicIndexSet, std::vector< T > >
  {
    typedef PersistentContainerVector< G, typename G::HierarchicIndexSet, std::vector< T > > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    /** \brief \deprecated typedef of class Grid */
    typedef typename Base::Grid GridType;

    ALUGridPersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
    : Base( grid.hierarchicIndexSet(), codim, value )
    {}
  };


  // PersistentContainer for ALUGrid
  // -------------------------------

  template< int dim, int dimworld, ALUGridElementType eltype, ALUGridRefinementType refinementtype, class Comm, class T >
  class PersistentContainer< ALUGrid< dim, dimworld, eltype, refinementtype, Comm >, T >
  : public ALUGridPersistentContainer< ALUGrid< dim, dimworld, eltype, refinementtype, Comm >, T >
  {
    typedef ALUGridPersistentContainer< ALUGrid< dim, dimworld, eltype, refinementtype, Comm >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    /** \brief \deprecated typedef of class Grid */
    typedef typename Base::Grid GridType;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
    : Base( grid, codim, value )
    {}
  };

  template< int dim, int dimworld, ALU2DSPACE ElementType elType, class T >
  class PersistentContainer< ALU2dGrid< dim, dimworld, elType >, T >
  : public ALUGridPersistentContainer< ALU2dGrid< dim, dimworld, elType >, T >
  {
    typedef ALUGridPersistentContainer< ALU2dGrid< dim, dimworld, elType >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    /** \brief \deprecated typedef of class Grid */
    typedef typename Base::Grid GridType;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
    : Base( grid, codim, value )
    {}
  };

  template< ALU3dGridElementType elType, class Comm, class T >
  class PersistentContainer< ALU3dGrid< elType, Comm >, T >
  : public ALUGridPersistentContainer< ALU3dGrid< elType, Comm >, T >
  {
    typedef ALUGridPersistentContainer< ALU3dGrid< elType, Comm >, T > Base;

  public:
    typedef typename Base::Grid Grid;
    typedef typename Base::Value Value;

    /** \brief \deprecated typedef of class Grid */
    typedef typename Base::Grid GridType;

    PersistentContainer ( const Grid &grid, int codim, const Value &value = Value() )
    : Base( grid, codim, value )
    {}
  };

} // namespace Dune

#endif // #ifndef DUNE_ALU_PERSISTENTCONTAINER_HH
