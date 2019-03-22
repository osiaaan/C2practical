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

// (c) Robert Kloefkorn 2004 - 2013
#ifndef ITERATOR_STI_H_INCLUDED
#define ITERATOR_STI_H_INCLUDED

namespace ALUGrid
{

  ////////////////////////////////////////////////////////////////////
  //
  // Schnittstelle des Iterationsobjekts vgl. Gamma, Helm, Johnson &
  // Vlissides: Design Patterns; Addison Wesley 
  // Die Schnittstellenbeschreibung wird sowohl polymorph als auch
  // in den verschiedenen Schablonen f"ur einfache Iterationsobjekte
  // s.a. Datei 'walk.h' verwendet.
  //
  ////////////////////////////////////////////////////////////////////
  template < class A > class IteratorSTI 
  {
  protected:  
    IteratorSTI () {}
  public :
    typedef A val_t;
    virtual ~IteratorSTI () {}
    virtual void first () = 0;
    virtual void next () = 0;
    virtual int done () const = 0;
    virtual int size () = 0;
    virtual val_t & item () const = 0;
    virtual IteratorSTI < A > * clone () const = 0;
  };

} // namespace ALUGrid

#endif // #ifndef ITERATOR_STI_H_INCLUDED
