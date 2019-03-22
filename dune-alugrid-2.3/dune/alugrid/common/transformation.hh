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

#ifndef DUNE_ALUGRID_TRANSFORMATION_HH
#define DUNE_ALUGRID_TRANSFORMATION_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

namespace Dune
{

  template< class ctype, int dimw >
  struct ALUGridTransformation
  {
    static const int dimension = dimw;

    typedef FieldVector< ctype, dimension > WorldVector;
    typedef FieldMatrix< ctype, dimension, dimension > WorldMatrix;

    ALUGridTransformation ( const WorldMatrix &matrix, const WorldVector &shift )
    : matrix_( matrix ),
      shift_( shift )
    {}

    WorldVector evaluate ( const WorldVector &x ) const
    {
      WorldVector y = shift_;
      matrix_.umv( x, y );
      return y;
    }

    WorldVector evaluateInverse ( const WorldVector &y ) const
    {
      // Note: We assume the matrix to be orthogonal, here
      WorldVector ys = y - shift_;
      WorldVector x;
      matrix_.mtv( ys, x );
      return x;
    }

  private:
    WorldMatrix matrix_;
    WorldVector shift_;
  };

}

#endif // #ifndef DUNE_ALUGRID_TRANSFORMATION_HH
