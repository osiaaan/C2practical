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

#ifndef ALUGRID_SRC_2d_VMMAP_H
#define ALUGRID_SRC_2d_VMMAP_H

#include <map>
#include <vector>

namespace ALU2DGrid
{

  template< int N, int NV >
  class Multivertexadapter
  {
  public:
    typedef Vertex < N > vertex_t;
    typedef Element < N, NV > element_t;

    typedef Macro < element_t > macroelement_t;

  private:
    typedef struct value
    {
      void * a;
      void * d;
      int b;
      int c;

      value(vertex_t * x = 0, int y = 0) : a(x), b(y), c(0) {}
     ~value() { }
    } val_t;

    typedef std::map< std::vector< vertex_t * >, val_t > map_t;

    std::vector< map_t > edmaps;
    std::vector< map_t > f4maps;

    Multivertexadapter(const Multivertexadapter &);
    Multivertexadapter & operator=(const Multivertexadapter &);

  public :
    Multivertexadapter();
   ~Multivertexadapter() {}

    void refresh(Listwalk < macroelement_t > & );

    vertex_t * find( vertex_t *, vertex_t *, int );

    vertex_t * find( vertex_t *, vertex_t *, vertex_t *, vertex_t *, int );

    void insert( vertex_t *, vertex_t *, vertex_t *, int );

    void insert( vertex_t *, vertex_t *, vertex_t *, vertex_t *, vertex_t *, int );
  };

} // namespace ALU2DGrid

#endif // #ifndef ALUGRID_SRC_2d_VMMAP_H
