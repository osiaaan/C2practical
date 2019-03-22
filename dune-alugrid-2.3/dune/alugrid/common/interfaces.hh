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

#ifndef DUNE_ALUGRID_INTERFACES_HH
#define DUNE_ALUGRID_INTERFACES_HH

#include <dune/common/typetraits.hh>

/** @file
  @author Robert Kloefkorn
  @brief Provides a Interfaces for detection of specific behavior
*/

namespace Dune {

  //! Tagging interface to indicate that Grid provides typedef ObjectStreamType
  struct HasObjectStream {};

  //! Helper template (implicit specialisation if GridImp exports an object
  //! stream
  template <bool hasStream, class GridImp, class DefaultImp>
  struct GridObjectStreamOrDefaultHelper {
    typedef typename GridImp::InStreamType     InStreamType;
    typedef typename GridImp::OutStreamType    OutStreamType;
  };
  
  //! Helper template (explicit specialisation if GridImp doesn't export an
  //! object stream -> DefaultImplementation is exported)
  template <class GridImp, class DefaultImp>
  struct GridObjectStreamOrDefaultHelper<false, GridImp, DefaultImp> {
    typedef DefaultImp InStreamType;
    typedef DefaultImp OutStreamType; 
  };

  //! Template to choose right Object stream type for a given class
  template <class GridImp, class DefaultImp>
  struct GridObjectStreamOrDefault 
  {
    typedef GridObjectStreamOrDefaultHelper<
                Conversion<GridImp, HasObjectStream>::exists, 
                GridImp, 
                DefaultImp> GridObjectStreamTraits; 
    
    typedef typename GridObjectStreamTraits :: InStreamType   InStreamType;  //  read  stream
    typedef typename GridObjectStreamTraits :: OutStreamType  OutStreamType; //  write stream 
  };

  //! Tagging interface to indicate that class is of Type DofManager 
  struct IsDofManager {};

  //! Tagging interface to indicate that Grid has HierarchicIndexSet  
  struct HasHierarchicIndexSet {};

} // end namespace Dune
#endif
