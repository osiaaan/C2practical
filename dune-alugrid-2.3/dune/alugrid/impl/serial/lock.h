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

// (c) bernhard schupp, 1997 - 1998
#ifndef LOCK_H_INCLUDED
#define LOCK_H_INCLUDED

#include <dune/alugrid/common/alugrid_assert.hh>
#include <cstring>
#include <iostream>

namespace ALUGrid
{

  class FSLock
  {
    char * _fname;

  public:
    FSLock ( const char * = "" );
   ~FSLock ();
  };

  inline FSLock::FSLock ( const char * name )
  : _fname( 0 )
  {
    _fname = new char[ std::strlen( name ) + 100 ];
    alugrid_assert ( _fname );
    sprintf (_fname, "%s.lock", name) ;
    FILE *fp = std::fopen( _fname, "w" );
    if( !fp )
    {
      delete[] _fname;
      _fname = 0;
      std::cerr << "WARNING (ignored): Could not create lock file." << std::endl;
    }
    else
    {
      // only test in debug mode 
#ifdef ALUGRIDDEBUG 
      int test = 
#endif
      std::fclose( fp );
      alugrid_assert (test == 0) ;
    }
  }

  inline FSLock::~FSLock ()
  {
    if( _fname )
    {
      int test = std::remove( _fname );
      if( test != 0 )
        std::cerr << "WARNING (ignored): Could not remove lock file." << std::endl;
      delete[] _fname;
      _fname = 0;
    }
  }

} // namespace ALUGrid

#endif // #ifndef LOCK_H_INCLUDED
