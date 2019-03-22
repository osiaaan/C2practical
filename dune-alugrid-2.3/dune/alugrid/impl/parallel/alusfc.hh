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

#ifndef ALUGRID_SFC_H_INCLUDED
#define ALUGRID_SFC_H_INCLUDED

#include <cmath>
#include <vector>

#include "mpAccess.h"

namespace ALUGridMETIS
{
  // Bresenham line drawing algorithm to partition the space filling curve
  template< class vertexmap_t, class connect_t, class vec_t >
  bool CALL_spaceFillingCurve(const ALUGrid::MpAccessGlobal& mpa, // communicator
                              const int numProcs,                 // number of partitions 
                              vertexmap_t& vertexMap,             // the space filling curve
                              connect_t&   connect,               // connectivity set
                              vec_t& graphSizes,                  // graph sizes to be communicated
                              const bool keepMapEntries )         // true if vertex entries should no be deleted 
  {
    // my rank 
    const int me = mpa.myrank();

    // clear connectivity set 
    connect.clear();

    typedef typename vertexmap_t :: iterator   iterator ;
    const iterator vertexEnd = vertexMap.end();
    long int sum = 0 ;
    // compute sum at first 
    for( iterator it = vertexMap.begin(); it != vertexEnd; ++ it ) 
    {
      sum += (*it).first.weight();
    }

    // clear map only when storeLinkageInVertices is not enabled 
    // since the vertices are still needed in that situation 
    const bool clearMap = ! ALUGrid :: Gitter :: storeLinkageInVertices && ! keepMapEntries ;

    const bool graphSizeCalculation = graphSizes.size() > 0 ;
    const int sizeOfVertexData = ALUGrid::LoadBalancer::GraphVertex::sizeOfData ;

    int destination = 0;
    long int d = -sum ;
    for( iterator it = vertexMap.begin(); it != vertexEnd; ++ it ) 
    {
      // increase destination if neccessary 
      if( d >= sum )
      {
        ++destination;
        d -= 2 * sum;
      }

      // get current rank
      const int source = (*it).second ;

      // set new rank information 
      (*it).second = destination ;
      // add weight 
      d += (2 * numProcs) * ((*it).first.weight());

      // add communication sizes of graph
      if( graphSizeCalculation ) 
        graphSizes[ destination ] += sizeOfVertexData ;

      // if the element currently belongs to me
      // then check the new destination 
      if( source == me && destination != me )
      {
        // insert into linkage set as send rank  
        connect.insert( ALUGrid::MpAccessLocal::sendRank( destination ) );
      }
      else if( source != me ) 
      {
        if( clearMap ) 
        {
          // mark element for delete 
          (*it).second = -1 ;
        }

        if( destination == me )
        {
          // insert into linkage set (receive ranks have negative numbers), see MpAccessLocal 
          connect.insert( ALUGrid::MpAccessLocal::recvRank( source ) );
        }
      }
    }

    if( clearMap )
    {
      // erase elements that are not further needed to save memory 
      for (iterator it = vertexMap.begin (); it != vertexEnd; )
      {
        // if element does neither belong to me not will belong to me, erase it 
        if( (*it).second < 0 )
        {
          vertexMap.erase( it++ );
        }
        else
          ++ it;
      }
    }

    alugrid_assert ( destination < numProcs );
    // return true if partitioning is ok, should never be false 
    return (destination < numProcs);
  } // end of simple sfc splitting without edges 

} // namespace ALUGridMETIS

#endif // #ifndef ALUGRID_SFC_H_INCLUDED
