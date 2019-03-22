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

#ifndef DUNE_ALU2DGRIDDATAHANDLE_HH
#define DUNE_ALU2DGRIDDATAHANDLE_HH

#include <iostream>

#include <dune/grid/common/adaptcallback.hh>

#include <dune/alugrid/2d/alu2dinclude.hh>

namespace ALU2DGrid
{

  // AdaptRestrictProlong2dImpl
  // --------------------------

  template< class Grid, class AdaptDataHandle >
  class AdaptRestrictProlong2dImpl
  : public AdaptRestrictProlong2d< Grid::dimensionworld,(Grid::elementType == ALU2DSPACE triangle ? 3 : 4) >
  {
    typedef Dune::MakeableInterfaceObject< typename Grid::template Codim< 0 >::Entity > EntityType;
    typedef typename EntityType::ImplementationType RealEntityType;
    typedef typename Dune::ALU2dImplTraits< Grid::dimensionworld, Grid::elementType >::HElementType HElementType;

    Grid & grid_;
    EntityType & reFather_;
    EntityType & reSon_;
    RealEntityType & realFather_;
    RealEntityType & realSon_;

    AdaptDataHandle &rp_;

    int maxlevel_;

  public:
    //! Constructor
    AdaptRestrictProlong2dImpl ( Grid &grid,
                                 EntityType &f, RealEntityType &rf,
                                 EntityType &s, RealEntityType &rs,
                                 AdaptDataHandle &rp )
     : grid_(grid)
      , reFather_(f)
      , reSon_(s)
      , realFather_(rf) 
      , realSon_(rs) 
      , rp_(rp) 
      , maxlevel_(-1) 
    {}

    virtual ~AdaptRestrictProlong2dImpl () 
    {}

    //! restrict data , elem is always the father 
    int preCoarsening ( HElementType &father )
    {
      maxlevel_ = std::max( maxlevel_, father.level() );
      //father.resetRefinedTag();
      realFather_.setElement( father );
      rp_.preCoarsening( reFather_ );

      return 0;
    }

    //! prolong data, elem is the father  
    int postRefinement ( HElementType &father )
    {
      maxlevel_ = std::max( maxlevel_, father.level()+1 );
      //father.resetRefinedTag();
      realFather_.setElement( father );
      rp_.postRefinement( reFather_ );

      return 0;
    }

    int maxLevel () const { return maxlevel_; }
  };

} // namespace ALU2DGrid

#endif // #ifndef DUNE_ALU2DGRIDDATAHANDLE_HH
