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

// (c) Robert Kloefkorn 2010
#include <config.h>

#include "ghost_info.h"

namespace ALUGrid
{

  template<int points>
  void MacroGhostInfoStorage<points>::
  inlineGhostElement(ObjectStream & os ) const
  {
   // local face number 
    os.writeObject( _fce );

    // global vertex number of the hexas vertices  
    for(int i=0; i<noVx; ++i) os.writeObject( _vx[i] );
    
    // global vertex numbers of the face not existing on this partition  
    for(int i=0; i<noFaceVx; ++i) 
    {
      os.writeObject( _vxface[i] );
      os.writeObject( _p[i][0] ); 
      os.writeObject( _p[i][1] ); 
      os.writeObject( _p[i][2] ); 
    }
  }

  template<int points>
  void MacroGhostInfoStorage<points>::
  readData(ObjectStream & os ) 
  {
    // read local face number
    os.readObject ( _fce );

    // read vertices of element
    for(int i=0; i<noVx; ++i)
    {
      os.readObject ( _vx[i] );
    }

#ifdef ALUGRIDDEBUG 
    for( int i=0; i<noVx; ++i )
    {
      for( int j=0; j<noVx; ++j)
      {
        if( j == i ) continue;
        alugrid_assert ( _vx[ i ] != _vx[ j ] );
      }
    }
#endif

    // read vertices of face an coordinates 
    for(int i=0; i<noFaceVx; ++i)
    {
      os.readObject ( _vxface[i] );
      alucoord_t (&pr) [3] = _p[i];

      os.readObject (pr[0]);
      os.readObject (pr[1]);
      os.readObject (pr[2]);
    }

    alugrid_assert ( _fce != invalidFace );
  }

  MacroGhostInfoHexa::
  MacroGhostInfoHexa(const Gitter::Geometric::hexa_GEO * hexa, 
                     const int fce)
  {
    alugrid_assert ( points == this->nop() );
    int oppFace = Gitter::Geometric::hexa_GEO::oppositeFace[fce];
    for(int vx=0; vx<points; vx++)
    {
      const Gitter::Geometric::VertexGeo * vertex = hexa->myvertex(oppFace,vx);
      this->_vxface[vx] = vertex->ident();
      const alucoord_t (&p) [3] = vertex->Point();
      this->_p[vx][0] = p[0];
      this->_p[vx][1] = p[1];
      this->_p[vx][2] = p[2];
    }
    
    for(int i=0; i<noVx; i++) 
    {
      this->_vx[i] = hexa->myvertex(i)->ident();
    }
    this->_fce = fce;
  }

  MacroGhostInfoTetra::MacroGhostInfoTetra ( const Gitter::Geometric::tetra_GEO *tetra, const int fce )
  {
    alugrid_assert ( points == this->nop() );
    const Gitter::Geometric::VertexGeo * vertex = tetra->myvertex( fce );
    alugrid_assert ( vertex );
    for(int vx=0; vx<points; ++vx)
    {
      this->_vxface[vx] = vertex->ident();
      const alucoord_t (&p) [3] = vertex->Point();
      this->_p[vx][0] = p[0];
      this->_p[vx][1] = p[1];
      this->_p[vx][2] = p[2];
    }
    
    for( int i = 0; i < noVx; ++i )
      this->_vx[i] = tetra->myvertex(i)->ident();

    this->_fce = tetra->orientation() ? -fce-1 : fce;
  }



  // Template Instantiation
  // ----------------------

  template class MacroGhostInfoStorage< 1 >;
  template class MacroGhostInfoStorage< 4 >;
  class MacroGhostInfoHexa; 
  class MacroGhostInfoTetra; 

} // namespace ALUGrid
