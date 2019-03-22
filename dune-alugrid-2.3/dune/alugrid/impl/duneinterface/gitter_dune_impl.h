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

#ifndef GITTER_DUNE_IMPL_H_INCLUDED
#define GITTER_DUNE_IMPL_H_INCLUDED

#include <iomanip>
#include <sstream>

#include "../serial/gitter_impl.h"
#include "../serial/lock.h"

namespace ALUGrid
{

  template < class A > class PureElementAccessIterator : public AccessIterator <A> 
  {
    public :
      Refcount ref;
      virtual IteratorSTI < A > * pureElementIterator (const A *) const = 0;
      virtual IteratorSTI < A > * pureElementIterator (const IteratorSTI < A > *) const = 0;
    public :

    // the only difference to AccessIterator is, that in the Constructors of this class 
    // the method pureElementIterator is called instead of iterator, which gives an iterator 
    // that doesn't iterator over periodic elements 

    class Handle
    : public AccessIterator< A >::Handle
    {
        // type of handle 
        typedef typename PureElementAccessIterator < A >::Handle ThisType;
      public :
        Handle ( AccessIterator< A > & );
        Handle ( const ThisType & );
        Handle (); 
    };

    protected :
      PureElementAccessIterator () {}
      virtual ~PureElementAccessIterator () { alugrid_assert (!ref); }
  };

  template < class A > class PureElementLeafIterator; 

  class GitterDuneBasis :  public virtual GitterBasis 
  {
    enum IndexType { no_index = 0 , hierarchic_index = 1, leaf_index = 3 };

  protected:
    // adaptation callback handler 
    AdaptRestrictProlongType * _arp;

    // call preCoarsening and postRefinement of arp
    virtual int preCoarsening  (Gitter::helement_STI &);
    virtual int postRefinement (Gitter::helement_STI &);
    
    // call preCoarsening and postRefinement of arp
    virtual int preCoarsening  (Gitter::hbndseg_STI &);
    virtual int postRefinement (Gitter::hbndseg_STI &);

    virtual void setAdaptRestrictProlongOp ( AdaptRestrictProlongType & arp );
    virtual void removeAdaptRestrictProlongOp ();

    // maxlevel of the grid 
    int maxlevel_;

    friend class PureElementLeafIterator < Gitter::helement_STI >;
    // return leafIterator using pureElement Iterators 
    virtual IteratorSTI < Gitter::helement_STI > * leafIterator (const Gitter::helement_STI *) = 0;
    virtual IteratorSTI < Gitter::helement_STI > * leafIterator (const IteratorSTI < Gitter::helement_STI > *) = 0;

  public:
    GitterDuneBasis() : _arp(0), maxlevel_(0) {}
    
    template <class ostream_t>
    void backupIndices  (ostream_t & out);

    template <class istream_t>
    void restoreIndices (istream_t & in );

    // write status of grid  
    void duneBackup (const char*); 

    // write status of grid  
    void duneBackup ( std::ostream &out );

    // read status of grid 
    void duneRestore (const char*);
   
    // read status of grid 
    void duneRestore ( std::istream &in );
   
    // done call notify and loadBalancer  
    bool duneAdapt (AdaptRestrictProlongType & arp);
  };

  class GitterDuneImpl : public GitterBasisImpl , public GitterDuneBasis 
  {
    // return LeafIterator which only iterates over elements 
    virtual IteratorSTI < Gitter::helement_STI > * leafIterator (const Gitter::helement_STI *);
    virtual IteratorSTI < Gitter::helement_STI > * leafIterator (const IteratorSTI < Gitter::helement_STI > *);

    friend class PureElementLeafIterator < Gitter::helement_STI >;
  public:

    //! constructor creating grid from std::istream 
    GitterDuneImpl ( std::istream &in, ProjectVertex *ppv = 0 )
    : GitterBasisImpl ( in, ppv )
    {}
    
    //! constructor creating grid from ObjectStream 
    inline GitterDuneImpl ( ObjectStream& in, ProjectVertex* ppv = 0 ) 
      : GitterBasisImpl ( in, ppv ) 
    {}
    
    //! constructor creating grid from macro grid file 
    inline GitterDuneImpl (const char *filename, ProjectVertex* ppv = 0 ) 
      : GitterBasisImpl (filename, ppv ) 
    {}
    
    //! constructor creating empty grid 
    inline GitterDuneImpl () 
      : GitterBasisImpl () 
    {}

    // compress memory of given grid and return new object (holding equivalent information)
    static GitterDuneImpl* compress( GitterDuneImpl* grd ) 
    {
      // only do the backup-restore thing if dlmalloc is enabled
      if( MyAlloc :: ALUGridUsesDLMalloc )
      {
        // backup stream 
        std::stringstream backup;
        // backup grid 
        grd->duneBackup( backup );
        delete grd; grd = 0;
        // free allocated memory (only works if all grids are deleted at this point)
        MyAlloc::clearFreeMemory ();
        // restore saved grid 
        grd = new GitterDuneImpl( backup );
        alugrid_assert ( grd );
        grd->duneRestore( backup );
      }
      return grd;
    }

  };


  // this LeafIterator only iterates over elements, i.e. tetra,hexa 
  template < class A > class PureElementLeafIterator : public MyAlloc {
    GitterDuneBasis * _grd;
    IteratorSTI < A > * _w;
    const A * _a;
    void * operator new (size_t) { return 0; }
    void operator delete (void *) {}
    inline PureElementLeafIterator ();
    public :
      inline PureElementLeafIterator (GitterDuneBasis &);
      inline PureElementLeafIterator (const PureElementLeafIterator < A > & );
      inline ~PureElementLeafIterator ();
      inline IteratorSTI < A > * operator -> () const;
      inline IteratorSTI < A > & operator * () const;
  };

  //
  //    #    #    #  #          #    #    #  ######
  //    #    ##   #  #          #    ##   #  #
  //    #    # #  #  #          #    # #  #  #####
  //    #    #  # #  #          #    #  # #  #
  //    #    #   ##  #          #    #   ##  #
  //    #    #    #  ######     #    #    #  ######
  //

  template <class ostream_t> 
  inline void GitterDuneBasis::backupIndices (ostream_t & out)
  {
    // get byte order of stream 
    out.put( RestoreInfo::systemByteOrder() );

    // backup indices, our index type is hierarchic_index 
    unsigned char indices = hierarchic_index;
    out.put( indices );

    enum { numOfIndexManager = Gitter::Geometric::BuilderIF:: numOfIndexManager };
    // store max indices 
    for(int i=0; i< numOfIndexManager; ++i)
      this->indexManager(i).backupIndexSet(out);

    { // backup index of elements 
      AccessIterator <helement_STI>::Handle ew (container ());
      for (ew.first (); ! ew.done (); ew.next ()) ew.item ().backupIndex (out);
    }

    // TODO: backup face and edge indices 
    
    {
      // backup index of vertices 
      LeafIterator < vertex_STI > w ( *this );
      for( w->first(); ! w->done(); w->next () ) w->item().backupIndex(out);
    }

    return;
  }

  template <class istream_t>
  inline void GitterDuneBasis ::restoreIndices (istream_t & in) 
  {
    // get byte order of stream 
    char byteOrder = in.get();  

    unsigned char indices = no_index;
    indices = in.get();

    // set VERBOSE to 20 and you have the indices value printed
#ifdef ALUGRIDDEBUG
    if( debugOption( 20 ) )
      std::cout << "INFO: GitterDuneBasis::restoreIndices.indices = " << (int)indices << std::endl;
#endif // #ifdef ALUGRIDDEBUG

    typedef Gitter::Geometric::BuilderIF  BuilderIF;
    enum { numOfIndexManager = BuilderIF::numOfIndexManager };
    
    // restore dune indices (see backUpIndices method)
    if(indices == hierarchic_index) 
    {
      // create vector, default all internal types for 
      // elements to vertices 
      RestoreInfo restoreInfo( byteOrder );
      
      for(int i=0; i< numOfIndexManager; ++i)
        this->indexManager(i).restoreIndexSet( in, restoreInfo );

      // will fail if numbering was changed 
      // and one forgot to apply changes here 
      alugrid_assert ( BuilderIF ::IM_Vertices+1 == 4 );

      // resize and reset 
      for(size_t i=0; i<restoreInfo.size(); ++i)
      {
        restoreInfo( i ).resize( this->indexManager(i).getMaxIndex(), true );
      }

      // restore index of elements 
      // mark all visited items as not a hole 
      {
        AccessIterator < helement_STI >:: Handle ew(container());
        for ( ew.first(); !ew.done(); ew.next()) ew.item().restoreIndex (in, restoreInfo);
      }
      // restore index of vertices
      // mark all visited items as not a hole 
      {
        LeafIterator < vertex_STI > w ( *this );
        for( w->first(); ! w->done(); w->next () ) w->item().restoreIndex(in, restoreInfo );
      }
      
      // reconstruct holes 
      {
        IndexManagerType& elementManager = this->indexManager(BuilderIF::IM_Elements);
        elementManager.generateHoles( restoreInfo( BuilderIF::IM_Elements ) );
      }
      
      // TODO indices for faces and edges 
      
      {
        IndexManagerType& vertexManager = this->indexManager(BuilderIF::IM_Vertices);
        vertexManager.generateHoles( restoreInfo( BuilderIF ::IM_Vertices ) );
      }
      return;
    }

    if( indices == leaf_index ) // convert indices to leafindices 
    {
      int idx = 0;
      PureElementLeafIterator < helement_STI > ew(*this);
      for ( ew->first(); !ew->done(); ew->next()) 
      {
        ew->item().setIndex( idx );
        ++idx;
      }
      this->indexManager( 0 ).setMaxIndex( idx );
#ifdef ALUGRIDDEBUG
      if( debugOption( 20 ) )
        std::cout << "INFO: GitterDuneBasis::restoreIndices created new leaf indices with size " << idx << "." << std::endl;
#endif // #ifdef ALUGRIDDEBUG
    }
    else
      std::cerr << "WARNING (ignored): indices (id = " << indices << ") not read in GitterDuneBasis::restoreIndices." << std::endl;
  }

  // wird von Dune verwendet 
  inline void GitterDuneBasis::duneBackup ( std::ostream &out )
  {
    // backup macro grid 
    container ().backupCMode ( out );
    // backup hierarchy 
    Gitter::backup ( out );
    // backup indices 
    backupIndices ( out );
  }

  // wird von Dune verwendet 
  inline void GitterDuneBasis::duneRestore ( std::istream &in )
  {
    // macro grid is created during grid creation
    // restore hierarchy 
    Gitter::restore (in);
    // restore indices 
    restoreIndices (in);
  }


  template < class A > inline PureElementAccessIterator < A >::
  Handle::Handle (AccessIterator < A > & f) 
   : AccessIterator < A >::Handle () 
  {
    this->removeObj();
    
    this->_fac = &f;
    this->_fac->ref ++;

    alugrid_assert ( this->_w == 0 );
    // this is the difference to the normal AccessIterator, we insert
    // pureElementIterator, all other things are the same   
    this->_w = this->_fac->iterator(this->_a);
    return;
  } 

  template < class A > inline PureElementAccessIterator < A >::Handle::
  Handle (const ThisType& p) 
    : AccessIterator < A >::Handle (p)
  {
  } 

  template < class A > inline PureElementAccessIterator < A >::Handle::Handle () 
    : AccessIterator < A >::Handle () {} 

  template < class A > PureElementLeafIterator < A >::PureElementLeafIterator () : _grd (0), _w (0) {
    return;
  }


  // new LEafIterator which only iterates over elements 
  template < class A > inline PureElementLeafIterator < A >::
  PureElementLeafIterator (GitterDuneBasis & g) : _grd (&g), _w (0) , _a(0) {
    _w = _grd->leafIterator (_a);
    return;
  }

  template < class A > inline PureElementLeafIterator < A >::PureElementLeafIterator (const PureElementLeafIterator < A > & x) : _grd (x._grd), _w (0) {
    _w = _grd->leafIterator (x._w);
    return;
  }
    
  template < class A > inline PureElementLeafIterator < A >::~PureElementLeafIterator () {
    if(_w) delete _w;
    return;
  }
    
  template < class A > inline IteratorSTI < A > * PureElementLeafIterator < A >::operator -> () const {
    return _w;
  }
    
  template < class A > inline IteratorSTI < A > & PureElementLeafIterator < A >::operator * () const {
    return * _w;
  }

  namespace {
  std::string ZeroPadNumber(int num)
  {
      std::ostringstream ss;
      ss << std::setw( 7 ) << std::setfill( '0' ) << num;
      return ss.str();
  }
  }

} // namespace ALUGrid

#endif // #ifndef GITTER_DUNE_IMPL_H_INCLUDED
