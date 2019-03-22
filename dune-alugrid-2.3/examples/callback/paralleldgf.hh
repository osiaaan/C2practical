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

#ifndef ALUGRID_PARALLEL_DGF_HH
#define ALUGRID_PARALLEL_DGF_HH

#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/alugrid/dgf.hh>

#include <dune/alugrid/common/structuredgridfactory.hh>

namespace Dune 
{
  template <class Grid> 
  class CreateParallelGrid ;  

  template < int dim, int dimworld, ALUGridElementType eltype, 
             ALUGridRefinementType refineType, class Comm > 
  class CreateParallelGrid< ALUGrid< dim, dimworld, eltype, refineType, Comm > >
  {
    typedef ALUGrid< dim, dimworld, eltype, refineType, Comm > Grid ;

  public:  
    static GridPtr< Grid > create( const std::string& filename ) 
    {
      // this only works for Cartesian grid using DGF's IntervalBlock
      if( eltype == simplex ) 
        return GridPtr< Grid >( filename );

#if ! HAVE_ALUGRID
      typedef StructuredGridFactory< Grid > SGF;
      return SGF :: createCubeGrid( filename );
#else 
      return GridPtr< Grid > (filename);
#endif // if ! HAVE_ALUGRID 
    }
  };
}

#endif
