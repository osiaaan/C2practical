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

#ifndef ALUGRID_SRC_SERIAL_REFCOUNT_HH
#define ALUGRID_SRC_SERIAL_REFCOUNT_HH

namespace ALUGrid
{

  // Einfacher Referenzenz"ahler mit cast-around-const
  // feature, der zum Z"ahlen der Referenzen auf Fl"achen
  // Kanten und Knoten verwendet wird. Vorteil: Objekte,
  // die einen Z"ahler dieser Klasse enthalten, werden
  // durch Inkrementierung bzw. Dekrementierung des Z"ahlers
  // nicht ver"andert (k"onnen also auch 'const' sein).

  class Refcount
  {
#ifdef ALUGRIDDEBUG
#ifdef DEBUG_ALUGRID
    // Der Globale Z"ahler soll helfen, nicht gel"oschte
    // Gitterobjekte oder Iteratorobjekte zu erkennen.
    // (Wird aber nur in den DEBUG-Versionen angelegt.) 
    //
    // Refcounting only turned on, if NDEBUG is not defined and
    // DEBUG_ALUGRID is defined 
    class Globalcount
    {
      int _c;
    public:
      Globalcount ();
      ~Globalcount ();
      void operator++ ( int ) const;
      void operator-- ( int ) const;
    };

    static Globalcount _g;
#endif // #ifdef DEBUG_ALUGRID
#endif // #ifdef ALUGRIDDEBUG

    mutable unsigned char _c;

  public:
    void reset () { _c = 0; }
    bool positive () const { return _c > 0; }
    Refcount ();
    ~Refcount ();
    int operator++ ( int ) const;
    int operator++ () const;
    int operator-- ( int) const;
    int operator-- () const;
    bool operator! () const;
    operator int () const;
  };

  class IteratorRefcount 
#ifdef ALUGRID_ITERATORS_WITH_MYALLOC
  : public Refcount 
#endif // #ifdef ALUGRID_ITERATORS_WITH_MYALLOC
  {
  public:
#ifndef ALUGRID_ITERATORS_WITH_MYALLOC
    void reset () { }
    bool positive () const { return false; }
    int operator++ ( int ) const { return 0; }
    int operator++ () const { return 0; }
    int operator-- ( int ) const { return 0; }
    int operator-- () const { return 0; }  
    bool operator! () const { return false; }
    operator int () const { return 0; }
#endif // #ifndef ALUGRID_ITERATORS_WITH_MYALLOC
  };

} // namespace ALUGrid

#endif // #ifndef ALUGRID_SRC_SERIAL_REFCOUNT_HH
