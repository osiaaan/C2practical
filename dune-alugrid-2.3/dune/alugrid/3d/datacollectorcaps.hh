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

#ifndef DUNE_ALUGRID_3D_DATACOLLECTORCAPS_HH
#define DUNE_ALUGRID_3D_DATACOLLECTORCAPS_HH

namespace ALUGrid
{

  namespace DataCollectorCaps
  {

    template< class DataCollector >
    class HasUserDefinedPartitioning
    {
      typedef char Small;
      struct Big { char dummy[2]; };

      template< class T, T > struct TypeCheck;

      typedef bool (DataCollector::*Method)() const;

      template< class T >
      static Small test ( TypeCheck< Method, &T::userDefinedPartitioning > * );
      template< class T >
      static Big test ( ... );

    public:
      static const bool v = (sizeof( test< DataCollector >( 0 ) ) == sizeof( Small ));
    };

    template< class DataCollector >
    class HasUserDefinedLoadWeights
    {
      typedef char Small;
      struct Big { char dummy[2]; };

      template< class T, T > struct TypeCheck;

      typedef bool (DataCollector::*Method)() const;

      template< class T >
      static Small test ( TypeCheck< Method, &T::userDefinedLoadWeights > * );
      template< class T >
      static Big test ( ... );

    public:
      static const bool v = (sizeof( test< DataCollector >( 0 ) ) == sizeof( Small ));
    };

  } // namespace DataCollectorCaps

} // namespace Dune

#endif // #ifndef DUNE_ALUGRID_3D_DATACOLLECTORCAPS_HH
