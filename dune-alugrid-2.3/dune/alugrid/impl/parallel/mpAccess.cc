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

#include <algorithm>
#include <iostream>

#include "mpAccess.h"

namespace ALUGrid
{

  void MpAccessLocal::printLinkage ( std::ostream &out ) const
  {
    const char* tag = ( symmetricLinkage() ) ? "linkage" : "sendLinkage" ;
    out << "  MpAccessLocal::printLinkage() " << myrank () << " (" << tag << ") -> ";
    typedef linkage_t :: const_iterator const_iterator;
    for (const_iterator i = _sendLinkage.begin (); i != _sendLinkage.end(); ++i )
      out << (*i).first << " ";
    out << std::endl;
    if( ! symmetricLinkage() )
    {
      out << "  MpAccessLocal::printLinkage() " << myrank () << " (recvLinkage) -> ";
      for (const_iterator i = _recvLinkage.begin (); i != _recvLinkage.end(); ++i )
        out << (*i).first << " ";
      out << std::endl;
    }
  }

  void MpAccessLocal::computeDestinations( const linkage_t& linkage, vector_t& dest )  
  {
    typedef linkage_t::const_iterator const_iterator ;
    dest.resize( linkage.size() );
    const const_iterator linkageEnd = linkage.end ();
    for( const_iterator i = linkage.begin (); i != linkageEnd; ++i )
    {
      dest[ (*i).second ] = (*i).first;
    }
  }

  int MpAccessLocal::insertRequestNonSymmetric( const std::set< int >& req )
  {
    const int me = myrank ();

    {
      typedef std::map< int, int >::iterator iterator ;
      typedef std::set< int >::const_iterator const_iterator;

      const iterator sendEnd = _sendLinkage.end ();
      const iterator recvEnd = _recvLinkage.end ();
      const const_iterator reqEnd = req.end ();
      int sendLink = 0 ;
      int recvLink = 0 ;
      for (const_iterator i = req.begin (); i != reqEnd; ++i )
      {
        const int val = (*i);
        // negative ranks mean receive ranks 
        if( val < 0 ) 
        {
          const int rank = recvRank( val );
          // if rank was not inserted, insert with current link number 
          if( rank != me && (_recvLinkage.find ( rank ) == recvEnd ) )
          {
            _recvLinkage.insert( std::make_pair( rank, recvLink++) );
          }
        }
        else 
        {
          const int rank = sendRank( val );
          // if rank was not inserted, insert with current link number 
          if( rank != me && (_sendLinkage.find ( rank ) == sendEnd ) )
          {
            _sendLinkage.insert( std::make_pair( rank, sendLink++) );
          }
        }
      }
    }

    // compute send destinations
    computeDestinations( _sendLinkage, _sendDest );
  
    // compute send destinations
    computeDestinations( _recvLinkage, _recvDest );

    // set pointers to recv stuff 
    _currentRecvLinkage = & _recvLinkage ;
    _currentRecvDest    = & _recvDest ;

    return sendLinks();
  }

  int MpAccessLocal::insertRequestSymmetric( const std::set< int >& req )
  {
    const int me = myrank ();

    {
      typedef std::map< int, int >::iterator iterator ;
      typedef std::set< int >::const_iterator const_iterator;

      const iterator linkageEnd = _sendLinkage.end ();
      const const_iterator reqEnd = req.end ();
      int link = 0 ;
      for (const_iterator i = req.begin (); i != reqEnd; ++i )
      {
        const int val = (*i);
        const int rank = (val < 0 ) ? recvRank( val ) : sendRank( val );
        // if rank was not inserted, insert with current link number 
        if( rank != me && (_sendLinkage.find ( rank ) == linkageEnd ) )
        {
          _sendLinkage.insert( std::make_pair( rank, link++) );
        }
      }
    }

    linkage_t().swap( _recvLinkage );
    // since linkage is symmetric, use sendLinkage for both 
    _currentRecvLinkage = & _sendLinkage;

    // compute send destinations
    computeDestinations( _sendLinkage, _sendDest );
  
    vector_t().swap( _recvDest );
    // since linkage is symmetric, use sendDest for both 
    _currentRecvDest = & _sendDest ;

    return sendLinks();
  }

  // insertRequestSymmetric needs a global communication 
  // this method is used to build the pattern for the loadBalancing 
  // where we don't know who is sending whom something 
  int MpAccessLocal::insertRequestSymmetricGlobalComm ( const std::set< int >& req )
  {
    const int me = myrank ();

    std::vector< int > out;
    out.reserve( req.size() );

    typedef linkage_t::iterator iterator ;
    typedef std::set< int >::const_iterator const_iterator;

    {
      const iterator linkageEnd = _sendLinkage.end ();
      const const_iterator reqEnd = req.end ();
      for (const_iterator i = req.begin (); i != reqEnd; ++i )
      {
        const int val = (*i);
        const int rank = (val < 0 ) ? recvRank( val ) : sendRank( val );
        if( rank != me && (_sendLinkage.find (rank) == linkageEnd ) )
          out.push_back ( rank );
      }
    }

    const iterator linkageEnd = _sendLinkage.end ();
    std::vector< std::vector< int > > in = gcollect (out);
    { 
      for (std::vector< int >::const_iterator i = out.begin (); i != out.end (); ++i )
      {
        if (_sendLinkage.find (*i) == linkageEnd )
        {
          int n = _sendLinkage.size ();
          _sendLinkage [*i] = n;
        }
      }
    }

    int cnt = 0;
    for( int i = 0; i < psize(); ++i )
    {
      if( std::find (in[ i ].begin(), in[ i ].end(), me) != in[ i ].end()  )
      {
        alugrid_assert (i != me);
        if (_sendLinkage.find (i) == linkageEnd ) 
        {
          int n = _sendLinkage.size ();
          _sendLinkage [i] = n;
          cnt ++; 
        }
      }
    }

    linkage_t().swap( _recvLinkage );
    // since linkage is symmetric, use sendLinkage for both 
    _currentRecvLinkage = & _sendLinkage;

    // compute send destinations
    computeDestinations( _sendLinkage, _sendDest );
  
    vector_t().swap( _recvDest );
    // since linkage is symmetric, use sendDest for both 
    _currentRecvDest = & _sendDest ;

    return sendLinks();
  }

} // namespace ALUGrid
