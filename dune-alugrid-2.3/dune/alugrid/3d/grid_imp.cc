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

#ifndef DUNE_ALUGRID_GRID_IMP_CC
#define DUNE_ALUGRID_GRID_IMP_CC

#if COMPILE_ALUGRID_INLINE == 0
#include <config.h>
#endif

// Dune includes
#include <dune/common/stdstreams.hh>

// Local includes
#include "entity.hh"
#include "iterator.hh"
#include "datahandle.hh"

#include "grid.hh"

#if COMPILE_ALUGRID_INLINE
#define alu_inline inline 
#else
#define alu_inline
#endif

namespace Dune
{

  template< class Comm >
  template< class GridType > 
  alu_inline 
  void ALU3dGridVertexList< Comm >:: 
  setupVxList(const GridType & grid, int level) 
  {
    // iterates over grid elements of given level and adds all vertices to
    // given list  
    enum { codim = 3 };

    VertexListType & vxList = vertexList_;

    unsigned int vxsize = grid.hierarchicIndexSet().size(codim); 
    if( vxList.size() < vxsize ) vxList.reserve(vxsize); 
    std::vector<int> visited_(vxsize);
    
    for(unsigned int i=0; i<vxsize; i++) 
    {
      visited_[i] = 0;
    }

    vxList.resize(0);
    
    const ALU3dGridElementType elType = GridType:: elementType;

    typedef ALU3DSPACE ALU3dGridLevelIteratorWrapper< 0, Dune::All_Partition, Comm > ElementLevelIteratorType;
    typedef typename ElementLevelIteratorType :: val_t val_t;
    
    typedef ALU3dImplTraits< elType, Comm > ImplTraits;
    typedef typename ImplTraits::IMPLElementType IMPLElementType;
    typedef typename ImplTraits::VertexType VertexType;

    enum { nVx = ElementTopologyMapping < elType > :: numVertices };
    
    ElementLevelIteratorType it ( grid, level, grid.nlinks() );

    int count = 0;
    for( it.first(); !it.done() ; it.next())
    {
      val_t & item = it.item();

      IMPLElementType * elem = 0; 
      if( item.first )
        elem = static_cast<IMPLElementType *> (item.first); 
      else if( item.second )
        elem = static_cast<IMPLElementType *> (item.second->getGhost().first); 

      alugrid_assert ( elem );

      for(int i=0; i<nVx; ++i)
      {
        VertexType * vx = elem->myvertex(i); 
        alugrid_assert ( vx );

        // insert only interior and border vertices 
        if( vx->isGhost() ) continue;

        const int idx = vx->getIndex(); 
        if(visited_[idx] == 0)
        {
          vxList.push_back(vx);
          ++count;
        }
        visited_[idx] = 1;
      }
    }
    alugrid_assert ( count == (int) vxList.size());;
    up2Date_ = true;    
  }


  template< class Comm >
  template< class GridType >
  alu_inline 
  void ALU3dGridLeafVertexList< Comm >::
  setupVxList(const GridType & grid) 
  {
    // iterates over grid elements of given level and adds all vertices to
    // given list  
    enum { codim = 3 };

    VertexListType & vxList = vertexList_;
    size_t vxsize = grid.hierarchicIndexSet().size(codim);
    if( vxList.capacity() < vxsize) vxList.reserve(vxsize);
    vxList.resize(vxsize);

    for(size_t i=0; i<vxsize; ++i) 
    {   
      ItemType & vx = vxList[i];
      vx.first  = 0;
      vx.second = -1;
    }

    const ALU3dGridElementType elType = GridType:: elementType;

    typedef ALU3DSPACE ALU3dGridLeafIteratorWrapper< 0, Dune::All_Partition, Comm > ElementIteratorType;
    typedef typename ElementIteratorType :: val_t val_t;
    
    typedef ALU3dImplTraits< elType, Comm > ImplTraits;
    typedef typename ImplTraits::IMPLElementType IMPLElementType;
    typedef typename ImplTraits::VertexType VertexType; 

    enum { nVx = ElementTopologyMapping < elType > :: numVertices };
    
    ElementIteratorType it ( grid, grid.maxLevel() , grid.nlinks() );

#ifdef ALUGRIDDEBUG
    int count = 0;
#endif
    for( it.first(); !it.done() ; it.next())
    {
      val_t & item = it.item();

      IMPLElementType * elem = 0; 
      if( item.first )
        elem = static_cast<IMPLElementType *> (item.first); 
      else if( item.second )
        elem = static_cast<IMPLElementType *> (item.second->getGhost().first); 
      
      alugrid_assert ( elem );
      int level = elem->level();

      for(int i=0; i<nVx; ++i)
      {
        VertexType * vx = elem->myvertex(i); 
        alugrid_assert ( vx );

        // insert only interior and border vertices 
        if( vx->isGhost() ) continue;

        const int idx = vx->getIndex(); 
        ItemType & vxpair = vxList[idx];
        if( vxpair.first == 0 )
        {
          vxpair.first  = vx; 
          vxpair.second = level; 
#ifdef ALUGRIDDEBUG
          ++ count;
#endif
        }
        // always store max level of vertex as grdi definition says
        else 
        {
          // set the max level for each vertex, see Grid definition
          if (vxpair.second < level) vxpair.second = level;
        }
      }
    }
    
    //std::cout << count << "c | s " << vxList.size() << "\n";
    // make sure that the found number of vertices equals to stored ones 
    //alugrid_assert ( count == (int)vxList.size() );
    up2Date_ = true;    
  }



  // ALU3dGrid
  // ---------

  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  const ALU3dGrid< elType, Comm > &
  ALU3dGrid< elType, Comm >::operator= ( const ALU3dGrid< elType, Comm > &other )
  {
    DUNE_THROW(GridError,"Do not use assignment operator of ALU3dGrid! \n");
    return (*this);
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  ALU3dGrid< elType, Comm >::~ALU3dGrid ()
  {
    delete communications_;
    delete vertexProjection_; 
    //delete bndPrj_;
    if( bndVec_ )
    {
      const size_t bndSize = bndVec_->size();
      for(size_t i=0; i<bndSize; ++i)
      {
        delete (*bndVec_)[i];
      } 
      delete bndVec_; bndVec_ = 0;
    }

    for(unsigned int i=0; i<levelIndexVec_.size(); i++) delete levelIndexVec_[i];
    delete globalIdSet_; globalIdSet_ = 0;
    delete leafIndexSet_; leafIndexSet_ = 0;
    delete sizeCache_; sizeCache_ = 0;

    /*
    if(myGrid().container().iterators_attached())
    {
      dwarn << "WRANING: There still exists instances of iterators giving access to this grid which is about to be removed! in: " << __FILE__ << " line: " << __LINE__ << std::endl;
    }
    */
    delete mygrid_; mygrid_ = 0;
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  int ALU3dGrid< elType, Comm >::size ( int level, int codim ) const
  {
    // if we dont have this level return 0 
    if( (level > maxlevel_) || (level < 0) ) return 0;

    alugrid_assert ( codim >= 0);
    alugrid_assert ( codim < dimension+1 );

    alugrid_assert ( sizeCache_ );
    return sizeCache_->size(level,codim);
  }


  template< ALU3dGridElementType elType, class Comm >
  size_t ALU3dGrid< elType, Comm >::numBoundarySegments () const
  {
    return myGrid().numMacroBndSegments();
  }


  // --size 
  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  int ALU3dGrid< elType, Comm >::size ( int level, GeometryType type ) const 
  {
    if(elType == tetra && !type.isSimplex()) return 0; 
    if(elType == hexa  && !type.isCube   ()) return 0; 
    return size( level, dimension - type.dim() );
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  int ALU3dGrid< elType, Comm >::size ( int codim ) const
  {
    alugrid_assert ( codim >= 0 );
    alugrid_assert ( codim <= dimension );

    alugrid_assert ( sizeCache_ );
    return sizeCache_->size(codim);
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  int ALU3dGrid< elType, Comm >::size ( GeometryType type ) const
  {
    if(elType == tetra && !type.isSimplex()) return 0; 
    if(elType == hexa  && !type.isCube   ()) return 0; 
    return size( dimension - type.dim() );
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  int ALU3dGrid< elType, Comm >::ghostSize ( int codim ) const
  {
    return ( ghostCellsEnabled() && codim == 0 ) ? 1 : 0 ;
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  int ALU3dGrid< elType, Comm >::ghostSize ( int level, int codim ) const
  {
    return ghostSize( codim );
  }

  // calc all necessary things that might have changed 
  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::setMaxLevel ( int mxl )
  {
    maxlevel_ = std::max(maxlevel_,mxl);
  }

  
  // calc all necessary things that might have changed 
  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::updateStatus()
  {
    calcMaxLevel();
    calcExtras();
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::calcMaxLevel ()
  {
    // old fashioned way 
    int testMaxLevel = 0;
    typedef ALU3DSPACE ALU3dGridLeafIteratorWrapper< 0, All_Partition, Comm > IteratorType;
    IteratorType w (*this, maxLevel(), nlinks() );
    
    typedef typename IteratorType :: val_t val_t ; 
    typedef typename ALU3dImplTraits< elType, Comm >::IMPLElementType IMPLElementType;  
    
    for (w.first () ; ! w.done () ; w.next ())
    {
      val_t & item = w.item();

      IMPLElementType * elem = 0; 
      if( item.first )
        elem = static_cast<IMPLElementType *> (item.first); 
      else if( item.second )
        elem = static_cast<IMPLElementType *> (item.second->getGhost().first); 

      alugrid_assert ( elem );
      
      int level = elem->level();
      if(level > testMaxLevel) testMaxLevel = level;
    }
    maxlevel_ = comm().max( testMaxLevel );
    alugrid_assert ( maxlevel_ == comm().max( maxlevel_ ));
  }


  // --calcExtras
  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::calcExtras ()
  {
    // make sure maxLevel is the same on all processes ????
    //alugrid_assert ( maxlevel_ == comm().max( maxlevel_ ));

    if(sizeCache_) delete sizeCache_;
    sizeCache_ = new SizeCacheType (*this);

    // unset up2date before recalculating the index sets, 
    // becasue they will use this feature
    leafVertexList_.unsetUp2Date();
    for(size_t i=0; i<MAXL; ++i) 
    { 
      vertexList_[i].unsetUp2Date();
      levelEdgeList_[i].unsetUp2Date();
    }

    if( comm().size() > 1 )
    {
      for( int i = 0; i < dimension; ++i )
      {
        ghostLeafList_[i].unsetUp2Date();
        for(size_t l=0; l<MAXL; ++l) ghostLevelList_[i][l].unsetUp2Date();
      }
    }

    // update all index set that are already in use 
    for(size_t i=0; i<levelIndexVec_.size(); ++i)
    {
      if(levelIndexVec_[i]) 
        (*(levelIndexVec_[i])).calcNewIndex( this->template lbegin<0>( i ),
                                             this->template lend<0>( i ) );
    }

    if(leafIndexSet_) 
      leafIndexSet_->calcNewIndex( this->template leafbegin<0>(), this->template leafend<0>() );

    // build global ID set new (to be revised)
    if( globalIdSet_ ) globalIdSet_->updateIdSet();

    coarsenMarked_ = 0;
    refineMarked_  = 0;
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  const typename ALU3dGrid< elType, Comm >::Traits::LeafIndexSet &
  ALU3dGrid< elType, Comm >::leafIndexSet () const
  {
    if(!leafIndexSet_) leafIndexSet_ = new LeafIndexSetImp ( *this,
                                                             this->template leafbegin<0>(),
                                                             this->template leafend<0>() );
    return *leafIndexSet_;
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  const typename ALU3dGrid< elType, Comm >::Traits::LevelIndexSet &
  ALU3dGrid< elType, Comm >::levelIndexSet ( int level ) const 
  {
    // check if level fits in vector  
    alugrid_assert ( level >= 0 );
    alugrid_assert ( level < (int) levelIndexVec_.size() );

    if( levelIndexVec_[level] == 0 )
      levelIndexVec_[level] = 
        new LevelIndexSetImp ( *this, 
                               this->template lbegin<0> (level), 
                               this->template lend<0> (level), 
                               level );
    return *(levelIndexVec_[level]);
  }


  // global refine 
  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::globalRefine ( int refCount )
  {
    alugrid_assert ( (refCount + maxLevel()) < MAXL ); 
    
    for( int count = refCount; count > 0; --count )
    {
      const LeafIteratorType end = leafend();
      for( LeafIteratorType it = leafbegin(); it != end; ++it )
        mark( 1, *it );
      const bool refined = adapt();
      if( refined )
        postAdapt();
    }
  }

  // preprocess grid
  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  bool ALU3dGrid< elType, Comm >::preAdapt()
  {
    return (coarsenMarked_ > 0);
  }


  // adapt grid
  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  bool ALU3dGrid< elType, Comm >::adapt ()
  {
    bool ref = false;

    if( lockPostAdapt_ == true ) 
    {
      DUNE_THROW(InvalidStateException,"Make sure that postAdapt is called after adapt was called and returned true!");
    }

    bool mightCoarse = preAdapt();
    // if prallel run, then adapt also global id set
    if(globalIdSet_)
    {
      //std::cout << "Start adapt with globalIdSet prolong \n";
      int defaultChunk = newElementsChunk_;
      int actChunk     = refineEstimate_ * refineMarked_; 
  
      // guess how many new elements we get 
      int newElements = std::max( actChunk , defaultChunk ); 
      
      globalIdSet_->setChunkSize( newElements );
      ref = myGrid().duneAdapt(*globalIdSet_); // adapt grid 
    }
    else 
    {
      ref = myGrid().adaptWithoutLoadBalancing(); 
    }

    // in parallel this is different 
    if( this->comm().size() == 1 ) 
    {
      ref = ref && refineMarked_ > 0;
    }

    if(ref || mightCoarse)
    {
      // calcs maxlevel and other extras 
      updateStatus();

      // notify that postAdapt must be called 
      lockPostAdapt_ = true;
    }
    return ref;
  }

  // post process grid  
  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::postAdapt ()
  {
    { 
      // old fashioned way 
      typedef ALU3DSPACE ALU3dGridLeafIteratorWrapper< 0, All_Partition, Comm > IteratorType;
      IteratorType w (*this, maxLevel(), nlinks() );
      
      typedef typename IteratorType::val_t val_t;
      typedef typename ALU3dImplTraits< elType, Comm >::IMPLElementType IMPLElementType;  
      
      for (w.first () ; ! w.done () ; w.next ())
      {
        val_t & item = w.item();

        alugrid_assert ( item.first || item.second );
        IMPLElementType * elem = 0; 
        if( item.first )
          elem = static_cast<IMPLElementType *> (item.first); 
        else if( item.second )
        {
          elem = static_cast<IMPLElementType *>( item.second->getGhost().first );
          alugrid_assert ( elem );
        }
        if (elem->hasBeenRefined())
        {
          elem->resetRefinedTag();
          // on bisected grids its possible that not only leaf elements where added so
          // we have to move up the hierarchy to make sure that the refined tag on parents are also removed
          while (elem->up())
          {
            elem = static_cast<IMPLElementType *>(elem->up());
            elem->resetRefinedTag();
          }
        }
      }
    }

    // make that postAdapt has been called
    lockPostAdapt_ = false;
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  bool ALU3dGrid< elType, Comm >
    ::writeGrid_Ascii( const std::string filename, alu3d_ctype time , bool scientific ) const
  {     
    // write ascii only works for serial grids at the moment 
    if ( this->comm().size() > 1 )
    {
      DUNE_THROW(GridError,"ALU3dGrid::writeGrid_Ascii not implemented for parallel grids!");
    }
    
    std::ofstream file ( filename.c_str() );
    if(file)
    {     
      typedef typename ALU3dImplTraits< elType, Comm >::HElementType    HElementType;
      typedef typename ALU3dImplTraits< elType, Comm >::VertexType      VertexType;
      typedef typename ALU3dImplTraits< elType, Comm >::BNDFaceType     BNDFaceType;
      typedef typename ALU3dImplTraits< elType, Comm >::IMPLElementType IMPLElementType;
      typedef typename ALU3dImplTraits< elType, Comm >::HasFaceType     HasFaceType;
      typedef typename ALU3dImplTraits< elType, Comm >::GEOVertexType   GEOVertexType;
          
      ALU3DSPACE LeafIterator < HElementType > leafElements( myGrid() );
      file << "!" << elType2Name( elType ) << " Elements = " << leafElements->size() << std::endl;

      ALU3DSPACE LeafIterator < VertexType > leafVertices( myGrid() );
      { 
        file << std::endl;

        if( scientific ) 
        {
          // use scientific mode 
          file << std::scientific;
        }
      
        // write coordinates of the vertices 
        int vxsize = leafVertices->size();
        file << vxsize << std::endl;

        typedef double ShortVecType[3];
        ShortVecType * vxvec = new ShortVecType [vxsize];
        alugrid_assert ( vxvec );
        
        for( leafVertices->first(); !leafVertices->done() ; leafVertices->next() )
        {
          const GEOVertexType & vertex = 
            static_cast<GEOVertexType &> (leafVertices->item()); 
          const double (&p)[3] = vertex.Point();
          int vxidx = vertex.getIndex();
          double (&v)[3] = vxvec[vxidx];
          for(int i=0; i<3; ++i) v[i] = p[i];
        }

        for(int i=0; i<vxsize; ++i)
        {
          file << vxvec[i][0] << " " << vxvec[i][1] << " " << vxvec[i][2] << std::endl;
        }

        delete [] vxvec;
      }

      file << std::endl;
      // write element vertices 
      {
        const int novx = (elType == tetra) ? 4 : 8;
        file << leafElements->size() << std::endl;
        for( leafElements->first(); !leafElements->done() ; leafElements->next() )
        {
          IMPLElementType & item = static_cast<IMPLElementType &> (leafElements->item());
          for(int i=0; i<novx; ++i)
          {
            const int vxnum = item.myvertex(i)->getIndex();
            file << vxnum << " ";
          }
          file << std::endl;
        }
      }

      // write boundary faces 
      {
        file << std::endl;
        const int nofaces  = (elType == tetra) ? 4 : 6;
        int bndfaces = 0;
        for( leafElements->first(); !leafElements->done() ; leafElements->next() )
        {
          IMPLElementType & item = static_cast<IMPLElementType &> (leafElements->item());
          for(int i=0; i<nofaces; ++i)
          {
            std::pair < HasFaceType * , int > nbpair = item.myneighbour(i);
            if(nbpair.first->isboundary())
            {
              ++bndfaces;
            }
          }
        }
        file << bndfaces << std::endl;
      }
      // write boundary faces 
      {
        const int bndvxnum = (elType == tetra) ? 3 : 4;
        const int nofaces  = (elType == tetra) ? 4 : 6;
        for( leafElements->first(); !leafElements->done() ; leafElements->next() )
        {
          IMPLElementType & item = static_cast<IMPLElementType &> (leafElements->item());
          for(int i=0; i<nofaces; ++i)
          {
            std::pair < HasFaceType * , int > nbpair = item.myneighbour(i);
            if(nbpair.first->isboundary())
            {
              BNDFaceType * face = static_cast<BNDFaceType *> (nbpair.first);
              file << -(face->bndtype()) << " " << bndvxnum << " ";
              for(int j=0; j<bndvxnum; j++)
              {
                int vxnum = face->myvertex(0,j)->getIndex();
                file << vxnum << " ";
              }
              file << std::endl;
            }
          }
        }
      }

      {
        file << std::endl;

        // write coordinates of the vertices 
        int vxnum = 0;
        for( leafVertices->first(); !leafVertices->done() ; leafVertices->next() )
        {
          file << vxnum << " -1" << std::endl;
          ++vxnum;
        }
      }
    }
    return true;
  }


  template< ALU3dGridElementType elType, class Comm >
  template <GrapeIOFileFormatType ftype>
  alu_inline 
  bool ALU3dGrid< elType, Comm >
    ::writeGrid ( const std::string filename, alu3d_ctype time ) const
  {
    switch(ftype)
    {
      case xdr  : return writeGrid_Xdr(filename,time);
      case ascii: return writeGrid_Ascii(filename,time);
      default: derr << "Wrong file type in writeGrid method~ \n";
    }
    return false;
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  bool ALU3dGrid< elType, Comm >
    ::writeGrid_Xdr ( const std::string filename, alu3d_ctype time ) const
  {
    myGrid().duneBackup( filename.c_str() );
    return true;
  }


  template< ALU3dGridElementType elType, class Comm >
  inline bool ALU3dGrid< elType, Comm >
    ::writeMacroGrid ( const std::string path, const std::string name ) const
  {
    std::stringstream filename;
    filename << path << "/" << name << "." << comm().rank();

    std::ofstream macro( filename.str().c_str() );

    if( macro ) 
    {
      // dump distributed macro grid as ascii files 
      myGrid().container().backupCMode( macro );
    }
    else 
      std::cerr << "WARNING: couldn't open file `" <<  filename.str() << "' for writing!" << std::endl;

    return true;
  }


  template< ALU3dGridElementType elType, class Comm >
  template <GrapeIOFileFormatType ftype>
  alu_inline 
  bool ALU3dGrid< elType, Comm >
    ::readGrid ( const std::string filename, alu3d_ctype & time )
  {
    {
      std::string mName(filename);
      mName += ".macro";
      const char * macroName = mName.c_str();
    
      { //check if file exists 
        std::ifstream check ( macroName );
        if( !check 
            // only abort on rank 0
            // on all other ranks this can be empty 
            && comm().rank() == 0      
            ) 
          DUNE_THROW(GridError,"cannot read file " << macroName << "\n");
        check.close();
      }
    
      // if grid exists delete first 
      if( mygrid_ ) delete mygrid_;
      mygrid_ = createALUGrid( macroName );
    }


    alugrid_assert (mygrid_ != 0);

    // check for element type 
    this->checkMacroGrid ();
    
    myGrid().duneRestore(filename.c_str());

    // don't restore time 
    time = 0;
  
    // calculate new maxlevel 
    // calculate indices 
    updateStatus();

    // reset refinement markers 
    postAdapt();

    return true;
  }

  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::backup( std::ostream& stream ) const
  {
    // backup grid to given stream 
    myGrid().duneBackup( stream );
  }

  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::restore( std::istream& stream )
  {
    // if grid exists delete first 
    if( mygrid_ ) delete mygrid_;
    // create new grid from stream 
    mygrid_ = createALUGrid( stream );

    // check for grid 
    if( ! mygrid_ ) 
    {
      DUNE_THROW(InvalidStateException,"ALUGrid::restore failed");
    }

    // check for element type 
    this->checkMacroGrid ();
    
    // restore hierarchy from given stream 
    myGrid().duneRestore( stream );

    // calculate new maxlevel 
    // calculate indices 
    updateStatus();

    // reset refinement markers 
    postAdapt();
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::checkMacroGridFile ( const std::string filename )
  {
    if(filename == "") return;
    
    std::ifstream file(filename.c_str());
    if(!file) 
    {
      std::cerr << "Couldn't open file '" << filename <<"' !" << std::endl;
      DUNE_THROW(IOError,"Couldn't open file '" << filename <<"' !");
    }
    
    const std::string aluid((elType == tetra) ? "!Tetrahedra" : "!Hexahedra");
    const std::string oldAluId((elType == tetra) ? "!Tetraeder" : "!Hexaeder");
    std::string idline;
    std::getline(file,idline);
    std::stringstream idstream(idline);
    std::string id;
    idstream >> id;

    if(id == aluid ) 
    {
      return; 
    }
    else if ( id == oldAluId ) 
    {
      derr << "\nKeyword '" << oldAluId << "' is deprecated! Change it to '" << aluid << "' in file '" << filename<<"'! \n";
      return ; 
    }
    else 
    {
      std::cerr << "Delivered file '"<<filename<<"' does not contain keyword '" 
        << aluid << "'. Found id '" <<id<< "'. Check the macro grid file! Bye." << std::endl;
      DUNE_THROW(IOError,"Wrong file format! ");
    }
  }


  template< ALU3dGridElementType elType, class Comm >
  alu_inline 
  void ALU3dGrid< elType, Comm >::checkMacroGrid ()
  {
    typedef typename ALU3dImplTraits< elType, Comm >::HElementType HElementType;
    typedef ALU3DSPACE PureElementLeafIterator< HElementType > IteratorType;
    IteratorType w( this->myGrid()  ); 
    for (w->first () ; ! w->done () ; w->next ())
    {
      ALU3dGridElementType type = (ALU3dGridElementType) w->item().type();
      if( type != elType )
      {
        derr << "\nERROR: " << elType2Name(elType) << " Grid tries to read a ";
        derr << elType2Name(type) << " macro grid file! \n\n";
        alugrid_assert (type == elType);
        DUNE_THROW(GridError,"\nERROR: " << elType2Name(elType) << " Grid tries to read a " << elType2Name(type) << " macro grid file! ");
      }
    }
  }


  alu_inline 
  const char * elType2Name( ALU3dGridElementType elType )
  {
    switch( elType )
    {
      case tetra  : return "Tetrahedra";
      case hexa   : return "Hexahedra";
      case mixed  : return "Mixed";
      default     : return "Error";
    }
  }


#if COMPILE_ALUGRID_LIB 
  // Instantiation
  template class ALU3dGrid< hexa, ALUGridNoComm >;
  template class ALU3dGrid< tetra, ALUGridNoComm >;

  template bool ALU3dGrid< tetra, ALUGridNoComm >::readGrid< xdr > (const std::string, alu3d_ctype & );  
  template bool ALU3dGrid< tetra, ALUGridNoComm >::readGrid< ascii > (const std::string, alu3d_ctype & );  
  template bool ALU3dGrid< tetra, ALUGridNoComm >::writeGrid< xdr > (const std::string, alu3d_ctype ) const ;  
  template bool ALU3dGrid< tetra, ALUGridNoComm >::writeGrid< ascii > (const std::string, alu3d_ctype ) const ;  

  template bool ALU3dGrid< hexa, ALUGridNoComm >::readGrid< xdr > (const std::string, alu3d_ctype & );  
  template bool ALU3dGrid< hexa, ALUGridNoComm >::readGrid< ascii > (const std::string, alu3d_ctype & );  
  template bool ALU3dGrid< hexa, ALUGridNoComm >::writeGrid< xdr > (const std::string, alu3d_ctype ) const ;  
  template bool ALU3dGrid< hexa, ALUGridNoComm >::writeGrid< ascii > (const std::string, alu3d_ctype ) const ;

  template class ALU3dGrid< hexa, ALUGridMPIComm >;
  template class ALU3dGrid< tetra, ALUGridMPIComm >;

  template bool ALU3dGrid< tetra, ALUGridMPIComm >::readGrid< xdr > (const std::string, alu3d_ctype & );  
  template bool ALU3dGrid< tetra, ALUGridMPIComm >::readGrid< ascii > (const std::string, alu3d_ctype & );  
  template bool ALU3dGrid< tetra, ALUGridMPIComm >::writeGrid< xdr > (const std::string, alu3d_ctype ) const ;  
  template bool ALU3dGrid< tetra, ALUGridMPIComm >::writeGrid< ascii > (const std::string, alu3d_ctype ) const ;  

  template bool ALU3dGrid< hexa, ALUGridMPIComm >::readGrid< xdr > (const std::string, alu3d_ctype & );  
  template bool ALU3dGrid< hexa, ALUGridMPIComm >::readGrid< ascii > (const std::string, alu3d_ctype & );  
  template bool ALU3dGrid< hexa, ALUGridMPIComm >::writeGrid< xdr > (const std::string, alu3d_ctype ) const ;  
  template bool ALU3dGrid< hexa, ALUGridMPIComm >::writeGrid< ascii > (const std::string, alu3d_ctype ) const ;
#endif // #if COMPILE_ALUGRID_LIB 

} // end namespace Dune

#endif // end DUNE_ALUGRID_GRID_IMP_CC
