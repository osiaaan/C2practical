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

#ifndef DUNE_ALU3DGRIDMEMORY_HH
#define DUNE_ALU3DGRIDMEMORY_HH

#include <dune/alugrid/common/alugrid_assert.hh>
#include <cstdlib>
#include <vector>

namespace ALUGrid
{

  template< class T, int length >
  class ALUGridFiniteStack;

}

namespace Dune
{

  //! organize the memory management for entitys used by the NeighborIterator
  template <class Object>
  class ALUMemoryProvider
  {
    enum { maxStackObjects = 256 };
    typedef ::ALUGrid::ALUGridFiniteStack< Object *, maxStackObjects > StackType;

    StackType objStack_;

    typedef ALUMemoryProvider < Object > MyType;

    StackType &objStack () { return objStack_; }

  public:
    typedef Object ObjectType;

    //!default constructor 
    ALUMemoryProvider() {}

    //! do not copy pointers  
    ALUMemoryProvider(const ALUMemoryProvider<Object> & org)
      : objStack_( org.objStack_ )
    {}

    //! call deleteEntity 
    ~ALUMemoryProvider ();

    //! i.e. return pointer to Entity
    template <class FactoryType>
    ObjectType * getObject(const FactoryType &factory, int level);

    //! i.e. return pointer to Entity
    template <class FactoryType, class EntityImp>
    inline ObjectType * getEntityObject(const FactoryType& factory, int level , EntityImp * fakePtr ) 
    {
      if( objStack().empty() )
      {
        return ( new ObjectType(EntityImp(factory,level) )); 
      }
      else
      {
        return stackObject();
      }
    }

    //! return object, if created default constructor is used 
    ObjectType * getEmptyObject ();

    //! i.e. return pointer to Entity
    ObjectType * getObjectCopy(const ObjectType & org);

    //! free, move element to stack, returns NULL 
    void freeObject (ObjectType * obj);

  protected:
    inline ObjectType * stackObject() 
    {
      alugrid_assert ( ! objStack().empty() );
      // finite stack does also return object on pop
      return objStack().pop();
    }

  };


  //************************************************************************
  //
  //  ALUMemoryProvider implementation
  //
  //************************************************************************
  template <class Object> template <class FactoryType>
  inline typename ALUMemoryProvider<Object>::ObjectType * 
  ALUMemoryProvider<Object>::getObject
  (const FactoryType &factory, int level )
  {
    if( objStack().empty() )
    {
      return ( new Object (factory, level) ); 
    }
    else
    {
      return stackObject();
    }
  }

  template <class Object>
  inline typename ALUMemoryProvider<Object>::ObjectType * 
  ALUMemoryProvider<Object>::getObjectCopy
  (const ObjectType & org )
  {
    if( objStack().empty() )
    {
      return ( new Object (org) ); 
    }
    else
    {
      return stackObject();
    }
  }

  template <class Object>
  inline typename ALUMemoryProvider<Object>::ObjectType * 
  ALUMemoryProvider<Object>::getEmptyObject () 
  {
    if( objStack().empty() )
    {
      return new Object () ; 
    }
    else
    {
      return stackObject();
    }
  }

  template <class Object>
  inline ALUMemoryProvider<Object>::~ALUMemoryProvider()
  {
    StackType& objStk = objStack_;
    while ( ! objStk.empty() )
    {
      ObjectType * obj = objStk.pop();
      delete obj;
    }
  }

  template <class Object>
  inline void ALUMemoryProvider<Object>::freeObject(Object * obj)
  {
    StackType& stk = objStack();
    if( stk.full() ) 
      delete obj;
    else 
      stk.push( obj );
  }

#undef USE_FINITE_STACK

} // namespace Dune

#endif // #ifndef DUNE_ALU3DGRIDMEMORY_HH
