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

#ifndef LOADBALANCE_HH
#define LOADBALANCE_HH

#include <dune/grid/common/gridenums.hh>
#include <dune/grid/common/datahandleif.hh>

#if HAVE_ZOLTAN 
#include <zoltan.h>
#endif

template< class Grid >
class SimpleLoadBalanceHandle
: public Dune::LoadBalanceHandleIF< SimpleLoadBalanceHandle<Grid> >
{
  typedef SimpleLoadBalanceHandle This;
  typedef Dune::LoadBalanceHandleIF< This > Base;

public:
  static const int dimension = Grid :: dimension;

  template< int codim >
  struct Codim
  {
    typedef typename Grid :: Traits :: template Codim< codim > :: Entity Entity;
    typedef typename Grid :: Traits :: template Codim< codim > :: EntityPointer
      EntityPointer;
  };
  typedef typename Codim< 0 > :: Entity Element;

private:
  const Grid &grid_;

public:
  SimpleLoadBalanceHandle ( const Grid &grid )
  : grid_( grid ),
    angle_( 0 )
  {}

  bool userDefinedPartitioning () const
  {
    return true;
  }
  // return true if user defined load balancing weights are provided
  bool userDefinedLoadWeights () const
  {
    return false;
  }

  // returns true if user defined partitioning needs to be readjusted 
  bool repartition () 
  { 
    angle_ += 2.*M_PI/50.;
    return true;
  }
  // return load weight of given element 
  int loadWeight( const Element &element ) const 
  { 
    return 1;
  }
  // return destination (i.e. rank) where the given element should be moved to 
  // this needs the methods userDefinedPartitioning to return true
  int destination( const Element &element ) const 
  { 
    typename Element::Geometry::GlobalCoordinate w = element.geometry().center();
    if (w[0]*w[0]+w[1]*w[1] > 0.2 && this->grid_.comm().size()>0)
    {
      double phi=arg(std::complex<double>(w[0],w[1]));
      if (w[1]<0) phi+=2.*M_PI;
      phi += angle_;
      phi *= double(this->grid_.comm().size()-1)/(2.*M_PI);
      int p = int(phi) % (this->grid_.comm().size()-1);
      return p+1;
    }
    else
      return 0;
  }
private:
  double angle_;
};

#if HAVE_ZOLTAN 
template< class Grid >
class ZoltanLoadBalanceHandle
: public Dune::LoadBalanceHandleIF< ZoltanLoadBalanceHandle<Grid> >
{
  typedef ZoltanLoadBalanceHandle This;
  typedef Dune::LoadBalanceHandleIF< This > Base;

private:
  typedef typename Grid::GlobalIdSet GlobalIdSet;
  typedef typename GlobalIdSet::IdType GIdType;
  static const int dimension = Grid :: dimension;
  static const int NUM_GID_ENTRIES = 1;
  template< int codim >
  struct Codim
  {
    typedef typename Grid :: Traits :: template Codim< codim > :: Entity Entity;
    typedef typename Grid :: Traits :: template Codim< codim > :: EntityPointer EntityPointer;
  };

  struct ZoltanPartitioning{
    int changes; // 1 if partitioning was changed, 0 otherwise 
    int numGidEntries;  // Number of integers used for a global ID 
    int numLidEntries;  // Number of integers used for a global ID 
    int numExport;      // Number of vertices I must send to other processes
    int numImport;      // Number of vertices I must send to other processes
    unsigned int *importLocalGids;  // Global IDs of the vertices I must send 
    unsigned int *importGlobalGids;  // Global IDs of the vertices I must send 
    unsigned int *exportLocalGids;  // Global IDs of the vertices I must send 
    unsigned int *exportGlobalGids;  // Global IDs of the vertices I must send 
    int *importProcs;    // Process to which I send each of the vertices 
    int *exportProcs;    // Process to which I send each of the vertices 
    int *importToPart;
    int *exportToPart;
  };
  struct FixedElements {
    int fixed_entities;
    std::vector<int> fixed_GID;
    std::vector<int> fixed_Process;
    FixedElements() : fixed_GID(0), fixed_Process(0) {}
  };
  struct HGraphData {      /* Zoltan will partition vertices, while minimizing edge cuts */
    int numMyVertices;  /* number of vertices that I own initially */
    ZOLTAN_ID_TYPE *vtxGID;        /* global ID of these vertices */
    int numMyHEdges;    /* number of my hyperedges */
    int numAllNbors;    /* number of vertices in my hyperedges */
    ZOLTAN_ID_TYPE *edgeGID;       /* global ID of each of my hyperedges */
    int *nborIndex;     /* index into nborGID array of edge's vertices */
    ZOLTAN_ID_TYPE *nborGID;       /* Vertices of edge edgeGID[i] begin at nborGID[nborIndex[i]] */
    FixedElements fixed_elmts;
    HGraphData() : vtxGID(0), edgeGID(0), nborIndex(0), nborGID(0) {}
    ~HGraphData() { freeMemory();}
    void freeMemory() 
    {
      free(nborGID); 
      free(nborIndex);
      free(edgeGID);
      free(vtxGID);
    }
  };
public:
  typedef typename Codim< 0 > :: Entity Element;
  ZoltanLoadBalanceHandle ( const Grid &grid);
  ~ZoltanLoadBalanceHandle();

  bool userDefinedPartitioning () const
  {
    return true;
  }
  // return true if user defined load balancing weights are provided
  bool userDefinedLoadWeights () const
  {
    return false;
  }

  // returns true if user defined partitioning needs to be readjusted 
  bool repartition ()
  { 
    if (!first_)
    {
      Zoltan_LB_Free_Part(&(new_partitioning_.importGlobalGids), 
                   &(new_partitioning_.importLocalGids), 
                   &(new_partitioning_.importProcs), 
                   &(new_partitioning_.importToPart) );
      Zoltan_LB_Free_Part(&(new_partitioning_.exportGlobalGids), 
                   &(new_partitioning_.exportLocalGids), 
                   &(new_partitioning_.exportProcs), 
                   &(new_partitioning_.exportToPart) );
    }
    generateHypergraph();
    /******************************************************************
    ** Zoltan can now partition the vertices of hypergraph.
    ** In this simple example, we assume the number of partitions is
    ** equal to the number of processes.  Process rank 0 will own
    ** partition 0, process rank 1 will own partition 1, and so on.
    ******************************************************************/
    Zoltan_LB_Partition(zz_, // input (all remaining fields are output)
          &new_partitioning_.changes,        // 1 if partitioning was changed, 0 otherwise 
          &new_partitioning_.numGidEntries,  // Number of integers used for a global ID 
          &new_partitioning_.numLidEntries,  // Number of integers used for a local ID 
          &new_partitioning_.numImport,      // Number of vertices to be sent to me 
          &new_partitioning_.importGlobalGids,  // Global IDs of vertices to be sent to me 
          &new_partitioning_.importLocalGids,   // Local IDs of vertices to be sent to me 
          &new_partitioning_.importProcs,    // Process rank for source of each incoming vertex 
          &new_partitioning_.importToPart,   // New partition for each incoming vertex 
          &new_partitioning_.numExport,      // Number of vertices I must send to other processes
          &new_partitioning_.exportGlobalGids,  // Global IDs of the vertices I must send 
          &new_partitioning_.exportLocalGids,   // Local IDs of the vertices I must send 
          &new_partitioning_.exportProcs,    // Process to which I send each of the vertices 
          &new_partitioning_.exportToPart);  // Partition to which each vertex will belong 
    first_ = false;
    return (new_partitioning_.changes == 1);
  }
  // return load weight of given element 
  int loadWeight( const Element &element ) const 
  { 
    return -1; // not used
  }
  // return destination (i.e. rank) where the given element should be moved to 
  // this needs the methods userDefinedPartitioning to return true
  int destination( const Element &element ) const 
  { 
	  std::vector<int> elementGID(NUM_GID_ENTRIES);
    // GIdType id = globalIdSet_.id(element);
    // id.getKey().extractKey(elementGID);
	  elementGID[0] = grid_.macroView().macroId(element); //   element.impl().macroID();

    // add one to the GIDs, so that they match the ones from Zoltan
    transform(elementGID.begin(), elementGID.end(), elementGID.begin(), bind2nd(std::plus<int>(), 1));

    int p = int(grid_.comm().rank());

    for (int i = 0; i<new_partitioning_.numExport; ++i)
    {
      if (std::equal(elementGID.begin(),elementGID.end(), &new_partitioning_.exportGlobalGids[i*new_partitioning_.numGidEntries]) )
      {
        p = new_partitioning_.exportProcs[i];
        break;
      }
    }
    return p;
  }
private:
  void generateHypergraph();

  // ZOLTAN query functions
  static int get_number_of_vertices(void *data, int *ierr);
  static void get_vertex_list(void *data, int sizeGID, int sizeLID,
              ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
              int wgt_dim, float *obj_wgts, int *ierr);
  static void get_hypergraph_size(void *data, int *num_lists, int *num_nonzeroes,
                                  int *format, int *ierr);
  static void get_hypergraph(void *data, int sizeGID, int num_edges, int num_nonzeroes,
                             int format, ZOLTAN_ID_PTR edgeGID, int *vtxPtr,
                             ZOLTAN_ID_PTR vtxGID, int *ierr);
  static int get_num_fixed_obj(void *data, int *ierr);
  static void get_fixed_obj_list(void *data, int num_fixed_obj,
                                 int num_gid_entries, ZOLTAN_ID_PTR fixed_gids, int *fixed_part, int *ierr);

  const Grid &grid_;
  const GlobalIdSet &globalIdSet_;

  Zoltan_Struct *zz_;
  HGraphData hg_;
  ZoltanPartitioning new_partitioning_;
  bool first_;
};
#endif // if HAVE_ZOLTAN

#include "loadbalance_inline.hh"

#endif // #ifndef LOADBALNCE_HH
