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

// (c) bernhard schupp 1998
// modifications for dune 
// (c) robert kloefkorn 2007 
#ifndef MYALLOC_H_INCLUDED
#define MYALLOC_H_INCLUDED

#include <cstddef>
#include <string>

// if defined the memory allocation from dlmalloc is used
#define ALUGRID_USES_DLMALLOC 
// if defined, standard C++ new and delete are used
//#define DONT_USE_ALUGRID_ALLOC 

namespace ALUGrid
{
  class ALUGridException 
  {
  protected:
    ALUGridException () {}
  public:  
    virtual ~ALUGridException () {}
    virtual std::string what () const = 0 ;
  };

#ifndef DONT_USE_ALUGRID_ALLOC 

  class MyAlloc
  {
    // max number of storable items per stack 
    static const std::size_t MAX_HOLD_ADD;
    // overestimation factor 
    static const double MAX_HOLD_MULT ;

    // true if initialized 
    static bool _initialized ;

    // if true objects are not free, only pushed to stack 
    static bool _freeAllowed ;
    
  public :
    static const bool ALUGridUsesDLMalloc = 
#ifdef ALUGRID_USES_DLMALLOC
      true ;
#else 
      false ;
#endif

      class Initializer 
      {
        // initializer versucht, die statischen Objekte der Speicherverwaltung
        // vor allem anderen zu initialisieren, damit keine Fehler auftreten,
        // falls statische Objekte irgendwo Instanzen mit MyAlloc als Basis-
        // klasse angelegen.
        public :
          Initializer () ;
         ~Initializer () ;
      } ;
      
      class OutOfMemoryException : public ALUGridException 
      { 
      public:
        virtual std::string what () const { return "OutOfMemoryException"; }
      };
      friend class Initializer;

      // if called, freeing objects is allowed again 
      static void unlockFree(void *);

      // if called free of objects is not allowed 
      static void lockFree (void *); 

      // try to make free memory available for the system 
      static void clearFreeMemory () ;

      // return size of allocated memory 
      static size_t allocatedMemory () ;


    protected :
      MyAlloc () {}
     ~MyAlloc () {}

    public :

      // new version of operator new 
      void * operator new (size_t) throw (OutOfMemoryException) ;
      // corresponding version of operator delete 
      void operator delete (void *,size_t) ;
  } ;

  static MyAlloc :: Initializer allocatorInitializer ;

#else //  #ifndef DONT_USE_ALUGRID_ALLOC 

  // dummy class 
  class MyAlloc 
  {
  public:  
    // this is false here anyway 
    static const bool ALUGridUsesDLMalloc = false ;

    // if called, freeing objects is allowed again 
    inline static void unlockFree(void *) {}

    // if called free of objects is not allowed 
    inline static void lockFree (void *) {} 

    // try to make free memory available for the system 
    inline static void clearFreeMemory () {}

    // return size of allocated memory 
    inline static size_t allocatedMemory () { return 0; }

  };

#endif // #else //  #ifndef DONT_USE_ALUGRID_ALLOC 

} // namespace ALUGrid

#endif // #ifndef MYALLOC_H_INCLUDED
