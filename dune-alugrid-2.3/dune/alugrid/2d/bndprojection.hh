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

#ifndef DUNE_ALU2D_BNDPROJECTION_HH
#define DUNE_ALU2D_BNDPROJECTION_HH

#include <dune/alugrid/common/bndprojection.hh>

#include <dune/alugrid/2d/alu2dinclude.hh>

namespace Dune
{

  template< class Grid >
  class ALU2dGridBoundaryProjection
  : public ALU2DSPACE VtxProjection< Grid::dimensionworld,(Grid::elementType == ALU2DSPACE triangle ? 3 : 4) >
  {
    typedef ALU2DSPACE VtxProjection< Grid::dimensionworld,(Grid::elementType == ALU2DSPACE triangle ? 3 : 4) > Base;

  public:
    enum { ncoord = Base::ncoord };

    typedef typename Base::hbndel_t hbndel_t;
    typedef typename Base::helement_t helement_t;

    typedef typename Grid::DuneBoundaryProjectionType DuneBoundaryProjectionType;

    typedef typename DuneBoundaryProjectionType::CoordinateType CoordinateType;

    explicit ALU2dGridBoundaryProjection ( const Grid &grid )
    : grid_( grid )
    {}

    int operator() ( const hbndel_t *hbndel, const double local, double (&global)[ ncoord ] ) const
    {
      return callProjection( grid_.boundaryProjection( hbndel->segmentIndex() ), global );
    }

    int operator() ( const helement_t *helement, const double (&local)[ 2 ], double (&global)[ ncoord ] ) const
    {
      return callProjection( grid_.globalProjection(), global );
    }

  private:
    static int callProjection ( const DuneBoundaryProjectionType *prj, double (&global)[ ncoord ] )
    {
      if( prj ) 
      {
        CoordinateType x, y;
        for( int i = 0; i < ncoord; ++i )
          x[ i ] = global[ i ];
        y = (*prj)( x );
        for( int i = 0; i < ncoord; ++i )
          global[ i ] = y[ i ];
      }
      return 1;
    }

    const Grid &grid_;
  };

} // end namespace Dune 

#endif // #ifndef DUNE_ALU2D_BNDPROJECTION_HH
