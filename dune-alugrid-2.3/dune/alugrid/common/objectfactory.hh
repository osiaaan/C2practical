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

#ifndef DUNE_ALUGRIDOBJECTFACTORY_HH
#define DUNE_ALUGRIDOBJECTFACTORY_HH

#include <dune/alugrid/common/memory.hh>

#if defined USE_PTHREADS || defined _OPENMP  
#define USE_SMP_PARALLEL
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

#if HAVE_DUNE_FEM 
#include <dune/fem/misc/threads/threadmanager.hh>
#endif

namespace Dune 
{
  template <class InterfaceType>
  struct MakeableInterfaceObject ;

  template <class GridImp>
  class ALUGridObjectFactory
  {
    template <class OF, int codim>
    class ALUGridEntityFactory;

    /////////////////////////////////////////////////////
    //
    //  partial specialization of method getNewEntity
    //
    /////////////////////////////////////////////////////
    template <class GridObjectFactory>
    class ALUGridEntityFactory<GridObjectFactory,0>
    {
    public:  
      enum { codim = 0 };
      typedef typename GridImp :: template Codim<codim> :: Entity Entity;
      typedef MakeableInterfaceObject<Entity> EntityObject;
      typedef typename EntityObject :: ImplementationType EntityImp;
      
      inline static EntityObject * 
      getNewEntity (const GridObjectFactory& factory, int level)
      {
        return factory.entityProvider_.getEntityObject( factory, level, (EntityImp *) 0);
      }

      inline static void freeEntity(const GridObjectFactory& factory, EntityObject * e )
      {
        factory.entityProvider_.freeObject( e );
      }
    };

    template <class GridObjectFactory>
    class ALUGridEntityFactory<GridObjectFactory,1>
    {
    public:  
      enum { codim = 1 };
      typedef typename GridImp :: template Codim<codim> :: Entity Entity;
      typedef MakeableInterfaceObject<Entity> EntityObject;
      typedef typename EntityObject :: ImplementationType EntityImp;
      
      inline static EntityObject * 
      getNewEntity (const GridObjectFactory& factory, int level)
      {
        return factory.faceProvider_.getEntityObject( factory, level, (EntityImp *) 0);
      }

      inline static void freeEntity(const GridObjectFactory& factory, EntityObject * e )
      {
        factory.faceProvider_.freeObject( e );
      }
    };

    template <class GridObjectFactory>
    class ALUGridEntityFactory<GridObjectFactory,2>
    {
    public:  
      enum { codim = 2 };
      typedef typename GridImp :: template Codim<codim> :: Entity Entity;
      typedef MakeableInterfaceObject<Entity> EntityObject;
      typedef typename EntityObject :: ImplementationType EntityImp;
      
      inline static EntityObject * 
      getNewEntity (const GridObjectFactory& factory, int level)
      {
        return factory.edgeProvider_.getEntityObject( factory, level, (EntityImp *) 0);
      }

      inline static void freeEntity(const GridObjectFactory& factory, EntityObject * e )
      {
        factory.edgeProvider_.freeObject( e );
      }
    };
    
    template <class GridObjectFactory>
    class ALUGridEntityFactory<GridObjectFactory,3>
    {
    public:  
      enum { codim = 3 };
      typedef typename GridImp :: template Codim<codim> :: Entity Entity;
      typedef MakeableInterfaceObject<Entity> EntityObject;
      typedef typename EntityObject :: ImplementationType EntityImp;
      
      inline static EntityObject * 
      getNewEntity (const GridObjectFactory& factory, int level)
      {
        return factory.vertexProvider_.getEntityObject( factory, level, (EntityImp *) 0);
      }

      inline static void freeEntity(const GridObjectFactory& factory, EntityObject * e )
      {
        factory.vertexProvider_.freeObject( e );
      }
    }; // end of ALUGridEntityFactory

    enum { vxCodim = GridImp :: dimension };
  public:  
    typedef GridImp GridType;
    typedef ALUGridObjectFactory FactoryType;

    typedef MakeableInterfaceObject<typename GridType :: Traits::template Codim<0>::Entity> EntityObject;
    typedef MakeableInterfaceObject<typename GridType :: Traits::template Codim<1>::Entity> FaceObject;
    typedef MakeableInterfaceObject<typename GridType :: Traits::template Codim<2>::Entity> EdgeObject;
    typedef MakeableInterfaceObject<typename GridType :: Traits::template Codim< vxCodim >::Entity> VertexObject;

    typedef typename GridType :: LeafIntersectionIteratorImp  LeafIntersectionIteratorImp ;
    typedef typename GridType :: LevelIntersectionIteratorImp LevelIntersectionIteratorImp ;

    // declare friendship 
    friend class ALUGridEntityFactory<FactoryType,0>;
    friend class ALUGridEntityFactory<FactoryType,1>;
    friend class ALUGridEntityFactory<FactoryType,2>;
    friend class ALUGridEntityFactory<FactoryType,3>;

  protected:  
    typedef ALUMemoryProvider< EntityObject > EntityProvider;
    typedef ALUMemoryProvider< FaceObject >   FaceProvider;
    typedef ALUMemoryProvider< EdgeObject >   EdgeProvider;
    typedef ALUMemoryProvider< VertexObject > VertexProvider;

    mutable EntityProvider   entityProvider_;
    mutable FaceProvider     faceProvider_;
    mutable EdgeProvider     edgeProvider_;
    mutable VertexProvider   vertexProvider_;
  
    typedef ALUMemoryProvider< LeafIntersectionIteratorImp > LeafIntersectionIteratorProviderType;
    typedef ALUMemoryProvider< LevelIntersectionIteratorImp >   LevelIntersectionIteratorProviderType;

    mutable LeafIntersectionIteratorProviderType leafInterItProvider_;
    mutable LevelIntersectionIteratorProviderType levelInterItProvider_;

    const GridType& grid_ ;

#ifdef USE_SMP_PARALLEL
  public:  
#endif
    ALUGridObjectFactory( const ALUGridObjectFactory& other ) : grid_( other.grid_ ) {}

  public:
    const  GridType& grid() const { return grid_; }

    ALUGridObjectFactory( const GridType& grid ) : grid_( grid ) {}

    template <int codim> 
    inline MakeableInterfaceObject<typename GridType :: Traits::template Codim<codim>::Entity> * 
    getNewEntity ( int level = -1 ) const
    {
      return ALUGridEntityFactory<FactoryType,codim>::getNewEntity( *this, level);
    }

    template <int codim>
    inline void freeEntity (MakeableInterfaceObject<typename GridType :: Traits::template Codim<codim>::Entity> * en) const
    {
      ALUGridEntityFactory<FactoryType,codim>::freeEntity(*this, en);
    }
  
    LeafIntersectionIteratorImp& getIntersection( const int wLevel, const LeafIntersectionIteratorImp* ) const
    {
      return * (leafInterItProvider_.getObject( *this, wLevel )); 
    }

    LevelIntersectionIteratorImp& getIntersection(const int wLevel, const LevelIntersectionIteratorImp* ) const 
    {
      return * (levelInterItProvider_.getObject( *this, wLevel )); 
    }

    //! free intersection 
    void freeIntersection(LeafIntersectionIteratorImp  & it) const { leafInterItProvider_.freeObject( &it ); }
    void freeIntersection(LevelIntersectionIteratorImp & it) const { levelInterItProvider_.freeObject( &it ); }

    // return thread number 
    static inline int threadNumber()
    {
#ifdef _OPENMP
      return omp_get_thread_num();
#elif HAVE_DUNE_FEM 
      return Fem :: ThreadManager :: thread() ;
#else
      return 0;
#endif
    }

    // return maximal possible number of threads 
    static inline int maxThreads() {
#ifdef _OPENMP
      return omp_get_max_threads();
#elif HAVE_DUNE_FEM 
      return Fem :: ThreadManager :: maxThreads() ;
#else
      return 1;
#endif
    }
  }; /// end class ALUGridObjectFactory

}  // end namespace Dune 
#endif
