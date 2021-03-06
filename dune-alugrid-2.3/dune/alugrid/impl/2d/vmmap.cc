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

#include <config.h>

#include "grid.h"
#include "handle.h"
#include "vmmap.h"

namespace ALU2DGrid
{

  template < int N, int NV >
  Multivertexadapter < N, NV >::Multivertexadapter()
  {
    edmaps.reserve(30);
    f4maps.reserve(30);
    edmaps.push_back(map_t ()); 
    f4maps.push_back(map_t());
  }

  template < int N, int NV >
  void Multivertexadapter < N, NV >::refresh( Listwalk < macroelement_t > & walk) {

    while(edmaps.size() > 0) edmaps.pop_back();

    edmaps.push_back(map_t());

    while(f4maps.size() > 0) f4maps.pop_back();

    f4maps.push_back(map_t());

    
    for( walk.first(); ! walk.done(); walk.next()) {

      for(int e = 0; e < walk.getitem()->numfaces(); e ++) {

        vertex_t * edg [2];

        walk.getitem()->edge_vtx(e, edg);

        std::vector< vertex_t * > v;

        v.push_back(edg[0]);

        v.push_back(edg[1]);

        sort(v.begin(), v.end());

        edmaps[0][v].b ++;

      }

    }

  }

  template < int N, int NV >
  typename Multivertexadapter < N, NV >::vertex_t *
  Multivertexadapter < N, NV >::find( vertex_t * a, vertex_t * b, int l) {

    std::vector< vertex_t * > e;

    e.push_back(a);

    e.push_back(b);

    sort(e.begin(), e.end());

    if(! l < edmaps.size()) edmaps.push_back(map_t ());

    map_t & map = edmaps[l];

    typename map_t::iterator edge = map.find(e);

    return edge == map.end() ? 0 : (vertex_t *) (*edge).second.a;

  }

  template < int N, int NV >
  typename Multivertexadapter < N, NV >::vertex_t *
  Multivertexadapter < N, NV >::find(vertex_t *a, vertex_t *b, vertex_t *c, vertex_t *d, int l) {

    std::vector< vertex_t * > v;

    v.push_back(a); 

    v.push_back(b);

    v.push_back(c); 

    v.push_back(d);

    sort(v.begin(), v.end());

    if(! l < f4maps.size()) f4maps.push_back(map_t ());

    map_t & map = f4maps[l];

    typename map_t::iterator face = map.find(v);

    if(face == map.end()) return 0;

    else {

      vertex_t * hit = (vertex_t *) (*face).second.a;

      map.erase(face);

      return hit;

    }

  }

  template < int N, int NV >
  void Multivertexadapter < N, NV >::insert( vertex_t * a, vertex_t * b, 

          vertex_t * ev, int l) {

    std::vector< vertex_t * > e;

    e.push_back(a);

    e.push_back(b);

    sort(e.begin(), e.end());  

    edmaps[l][e] = val_t(ev,0);

  }

  template< int N, int NV >
  void Multivertexadapter< N, NV >
    ::insert ( vertex_t *a, vertex_t *b, vertex_t *c, vertex_t *d, vertex_t *cv, int l )
  {
    std::vector< vertex_t * > v;
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    v.push_back(d);
    std::sort( v.begin(), v.end() );
    f4maps[l][v] = val_t(cv,0);
  }



  // Template Instantiation
  // ----------------------

  template class Multivertexadapter < 2,3 >;
  template class Multivertexadapter < 3,3 >;
  template class Multivertexadapter < 2,4 >;
  template class Multivertexadapter < 3,4 >;

} // namespace ALU2DGrid
