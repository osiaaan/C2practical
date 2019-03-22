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

#include <config.h>

#include <fstream>

#include "../serial/gatherscatter.hh"
#include "gitter_dune_pll_impl.h"

namespace ALUGrid
{

  IteratorSTI < Gitter::helement_STI > * GitterDunePll::
  leafIterator (const helement_STI *)
  {
    return new Insert < PureElementAccessIterator < Gitter::helement_STI >::Handle,
      TreeIterator < Gitter::helement_STI, is_leaf < Gitter::helement_STI> > > (container ());
  }

  IteratorSTI < Gitter::helement_STI > * GitterDunePll ::
  leafIterator (const IteratorSTI < helement_STI > * p)
  {
    return new Insert < PureElementAccessIterator < Gitter::helement_STI >::Handle,
      TreeIterator < Gitter::helement_STI, is_leaf < Gitter::helement_STI> > >
      (*(const Insert < PureElementAccessIterator < Gitter::helement_STI >::Handle,
         TreeIterator < Gitter::helement_STI, is_leaf < Gitter::helement_STI> > > *) p);
  }

  bool GitterDunePll::duneNotifyNewGrid ()
  {
    LoadBalancer::DataBase db;
    return checkPartitioning( db, (GatherScatter*) 0 );
  }

  void GitterDunePll::duneNotifyMacroGridChanges ()
  {
    GitterPll::notifyMacroGridChanges ();
    rebuildGhostCells ();
  }

  // done call notify and loadBalancer  
  bool GitterDunePll::duneAdapt ( AdaptRestrictProlongType &arp )
  {
    this->setAdaptRestrictProlongOp(arp);
    const bool refined = this->adaptWithoutLoadBalancing();
    this->removeAdaptRestrictProlongOp ();
    return refined;
  }

  bool GitterDunePll::duneLoadBalance () 
  {
    return loadBalancerGridChangesNotify ( ( GatherScatterType* ) 0 );
  }

  // returns true if grid was repartitioned 
  bool GitterDunePll::duneLoadBalance (GatherScatterType & gs, AdaptRestrictProlongType & arp) 
  {
    // set restriction/prolongation operator 
    this->setAdaptRestrictProlongOp(arp);

    // do load balancing 
    const bool repartion = loadBalancerGridChangesNotify ( &gs );

    // remove restriction/prolongation operator 
    this->removeAdaptRestrictProlongOp ();
    return repartion;
  }

  std::pair< IteratorSTI < GitterPll::vertex_STI > *, IteratorSTI < GitterPll::vertex_STI > *> 
  GitterDunePll::borderIteratorTT (const vertex_STI * v, int link )
  {
    // return default vertex iterator 
    return this->iteratorTT(v, link);
  }

  std::pair< IteratorSTI < GitterPll::hedge_STI > *, IteratorSTI < GitterPll::hedge_STI > *> 
  GitterDunePll::borderIteratorTT (const hedge_STI * e, int link )
  {
    // return edge iterator over all edges 
    is_def_true< hedge_STI > * s = 0;
    return this->createEdgeIteratorTT(s, link);
  }

  std::pair< IteratorSTI < GitterPll::hface_STI > *, IteratorSTI < GitterPll::hface_STI > *> 
  GitterDunePll::borderIteratorTT (const hface_STI * f, int link )
  {
    // return face iterator over all faces 
    is_def_true< hface_STI > rule;
    return this->createFaceIteratorTT( rule , link);
  }

  std::pair< IteratorSTI < GitterPll::hface_STI > *, IteratorSTI < GitterPll::hface_STI > *> 
  GitterDunePll::leafBorderIteratorTT (const hface_STI * f, int link )
  {
    // return face iterator over all faces that are 
    // leaf faces in the DUNE context
    is_leaf_entity < hface_STI > rule;
    return this->createFaceIteratorTT( rule , link);
  }

  std::pair< IteratorSTI < GitterPll::hface_STI > *, IteratorSTI < GitterPll::hface_STI > *> 
  GitterDunePll::levelBorderIteratorTT (const hface_STI * f, int link , int level)
  {
    // return face iterator over faces with given level 
    any_has_level < hface_STI > rule(level);
    return this->createFaceIteratorTT( rule, link);
  }

  template <class ObjectStreamType, class HItemType> 
  void GitterDunePll::sendSlaves (
      ObjectStreamType & sendBuff, 
      HItemType * fakeItem ,
      GatherScatterType & dataHandle, const int link )
  {
    // temporary buffer 
    SmallObjectStream osTmp; 

    std::pair< IteratorSTI < HItemType > *, IteratorSTI < HItemType > *> 
      a = borderIteratorTT (fakeItem, link ); //ueber alle meine Slave-Knoten 
   
    IteratorSTI < HItemType > & iter = *(a.second);
    for (iter.first (); ! iter.done (); iter.next ()) 
    {
      HItemType & item = iter.item();

      // gather all data on slaves 
      if( dataHandle.containsItem(item) ) 
      {
        // write marker that show data is transmitted 
        sendBuff.writeObject( transmittedData );

        // reset read and write position
        osTmp.clear();
        // write data to fake buff to determine size of data package
        dataHandle.sendData(osTmp,item);

        int s = osTmp.size();
        // first write size 
        sendBuff.writeObject(s);
        // then write bytes 
        sendBuff.writeStream(osTmp);
      } 
      else 
      {
        // write noData marker 
        sendBuff.writeObject( noData );
      }
    }

    delete a.first;
    delete a.second;      

    return;
  }

  template <class HItemType, class CommMapType>
  GitterDunePll::DataBufferType& 
  GitterDunePll::
  getCommunicationBuffer( HItemType& item, CommMapType& commMap, const int nCommBuff )
  {
    DataBufferType& commBuff = commMap[ &item ];
    if( (int) commBuff.size() != nCommBuff ) 
      commBuff.resize( nCommBuff );
    return commBuff;
  }

  template <class ObjectStreamType, class HItemType, class CommBuffMapType> 
  void GitterDunePll::unpackOnMaster (
      ObjectStreamType & recvBuff, 
      CommBuffMapType& commBuffMap,
      HItemType * determType,
      GatherScatterType & dataHandle ,
      const int nl, const int link )
  {
    int hasdata;

    typedef SmallObjectStream BufferType;
    typedef std::vector< BufferType > DataBufferType;

    std::pair< IteratorSTI < HItemType > *, IteratorSTI < HItemType > *> 
      a = borderIteratorTT (determType, link);
   
    IteratorSTI < HItemType > & iter = *(a.first);

    // for all master items 
    for (iter.first (); ! iter.done (); iter.next ()) 
    {
      HItemType & item = iter.item();
     
      // read data marker 
      recvBuff.readObject(hasdata);
      
      // get comm buffers 
      DataBufferType & data = getCommunicationBuffer( item, commBuffMap, nl + 1 );

      // only gather master data once 
      if ( dataHandle.containsItem( item ) ) 
      {
        // pack master data 
        BufferType & mData = data[ nl ]; 
        // reset read and write position
        mData.clear();
          
        // write master data to fake buffer 
        dataHandle.sendData(mData,item);
      }

      // if data has been send, read data 
      if (hasdata != noData) 
      {
        // pack slave data to tmnp buffer 
        BufferType & slaveBuff = data[link]; 
        // reset read and write position
        slaveBuff.clear();

        int dataSize; 
        recvBuff.readObject(dataSize);
        // read dataSize bytes from recvBuff and write to slaveStream 
        recvBuff.readStream(slaveBuff, dataSize);
      }
    }

    delete a.first;
    delete a.second;

    return;
  }

  template <class ObjectStreamType, class HItemType, class CommBuffMapType > 
  void GitterDunePll::sendMaster (
      ObjectStreamType & sendBuff, 
      CommBuffMapType& commBuffMap,
      HItemType * determType,
      GatherScatterType & dataHandle ,
      const int nl , 
      const int myLink )
  {
    typedef SmallObjectStream BufferType;
    typedef std::vector< BufferType > DataBufferType;

    std::pair< IteratorSTI < HItemType > *, IteratorSTI < HItemType > *> 
      a = borderIteratorTT (determType , myLink ); //ueber alle meine Slave-Knoten
   
    IteratorSTI < HItemType > & iter = *(a.first);

    // create new link vector 
    std::vector< int > newLink( nl );
    for(int link=0; link<nl; ++link) 
    {
      newLink[ link ] = link;
    }

    // if myLink == link then write master data
    // instead of data of link 
    // we do not send link i its own data
    newLink[myLink] = nl;

    // for all master items 
    for (iter.first (); ! iter.done (); iter.next ()) 
    {
      HItemType & item = iter.item();

      // get comm buffer 
      //DataBufferType & dataBuff = item.commBuffer();
      DataBufferType & dataBuff = getCommunicationBuffer( item, commBuffMap, nl + 1);
      
      // scatter on master 
      if ( dataHandle.containsItem( item ) ) 
      {
        for(int link = 0; link<nl; ++link)
        {
          BufferType & localBuff = dataBuff[link];

          // check if stream has been read, if not scatter data 
          // this will unpack data on master only once 
          if( localBuff.validToRead() ) 
          {
            dataHandle.recvData(localBuff, item);
          }
        }
      } 
     
      // pack for slaves 
      {
        // write data marker 
        sendBuff.writeObject(transmittedData);

        for(int link = 0; link<nl; ++link)
        {
          // use new link to send master data to link we are sending for 
          BufferType & localBuff = dataBuff[ newLink[link] ];
          // get size 
          int s = localBuff.size();
          sendBuff.writeObject(s);
          // if buffer size > 0 write hole buffer to stream 
          if( s > 0 ) sendBuff.writeStream( localBuff );
        }
      } 
    }

    delete a.first;
    delete a.second;     

    return;
  }

  template <class ObjectStreamType, class HItemType> 
  void GitterDunePll::unpackOnSlaves (
      ObjectStreamType & recvBuff, 
      HItemType * determType,
      GatherScatterType & dataHandle ,
      const int nOtherlinks, const int myLink)
  {
    int hasdata;

    std::pair< IteratorSTI < HItemType > *, IteratorSTI < HItemType > *> 
      a = borderIteratorTT (determType, myLink );

    // get slave iterator 
    IteratorSTI < HItemType > & iter = *(a.second);
    
    for (iter.first (); ! iter.done (); iter.next ()) 
    {
      // read data marker 
      recvBuff.readObject(hasdata);

      if (hasdata != noData) 
      {
        HItemType & item = iter.item();
        if( dataHandle.containsItem( item ) )
        {
          // for number of recived data, do scatter 
          for(int link = 0; link<nOtherlinks; ++link)
          {
            int s;
            recvBuff.readObject(s);
            if(s > 0) dataHandle.recvData(recvBuff, item );
          }
        }
        else 
        {
          // for number of recived data, do remove  
          for(int link = 0; link<nOtherlinks; ++link)
          {
            int s;
            recvBuff.readObject(s); 
            // if no data for link exists, s == 0
            // otherwise remove s bytes from stream by increasing 
            // read byte counter 
            if(s > 0) recvBuff.removeObject( s );
          }
        }
      }
    }
    delete a.first;
    delete a.second;
  }

  void GitterDunePll::sendFaces (
      ObjectStream & sendBuff, 
      IteratorSTI < hface_STI > * iter , 
      GatherScatterType & faceData )
  {
    // temporary object buffer  
    SmallObjectStream osTmp; 
    
    for (iter->first (); ! iter->done (); iter->next ()) 
    {
      hface_STI & face = iter->item();
      if ( faceData.containsItem( face ) ) 
      {
        sendBuff.writeObject(transmittedData);
        osTmp.clear();
        faceData.sendData(osTmp, face );

        int size = osTmp.size();
        // determin size of data to be able to remove 
        sendBuff.writeObject(size);
        if( size > 0 ) sendBuff.writeStream( osTmp );
      }
      else 
      {
        sendBuff.writeObject(noData);
      }
    }
  }

  void GitterDunePll::unpackFaces (
      ObjectStream & recvBuff, 
      IteratorSTI < hface_STI > * iter , 
      GatherScatterType & faceData )
  {
    int hasdata;
    for (iter->first (); ! iter->done (); iter->next ()) 
    {
      recvBuff.readObject(hasdata);
      if (hasdata != noData) 
      {
        hface_STI & face = iter->item();
        int size; 
        recvBuff.readObject(size); 
        if( size > 0 )
        {
          // if entity is not contained just remove data from stream 
          if ( faceData.containsItem( face ) ) 
            faceData.recvData(recvBuff , face );
          else 
            recvBuff.removeObject( size );
        }
      }
    }
  }

  ////////////////////////////////////////////////////////
  //
  // communication of higher codim data (vertices,edges,faces)
  //
  ////////////////////////////////////////////////////////
  void GitterDunePll::doBorderBorderComm( 
    std::vector< ObjectStream > & osvec ,
    GatherScatterType & vertexData , 
    GatherScatterType & edgeData,  
    GatherScatterType & faceData )
  {
    const int nl = mpAccess ().nlinks ();
    
    const bool containsVertices = vertexData.contains(3,3);
    const bool containsEdges    = edgeData.contains(3,2);
    const bool containsFaces    = faceData.contains(3,1);

    const bool haveVerticesOrEdges = containsVertices || containsEdges;
     
    alugrid_assert ((debugOption (5) && containsVertices) ? (std::cout << "**INFO GitterDunePll::borderBorderComm (): (containsVertices)=true " << std::endl, 1) : 1);
    alugrid_assert ((debugOption (5) && containsEdges)    ? (std::cout << "**INFO GitterDunePll::borderBorderComm (): (containsEdges)=true " << std::endl, 1) : 1);
    alugrid_assert ((debugOption (5) && containsFaces)    ? (std::cout << "**INFO GitterDunePll::borderBorderComm (): (containsFaces)=true " << std::endl, 1) : 1);
     
    // buffers for vertex and edge master-slave communication
    std::map< vertex_STI*, DataBufferType > vertexCommMap;
    std::map< hedge_STI* , DataBufferType > edgeCommMap;

    {
      // gather all data from slaves 
      for( int link = 0; link < nl; ++link )
      {
        ObjectStream & sendBuff = osvec[link];
        sendBuff.clear();

        if (containsVertices)
        {
          vertex_STI * determType = 0;
          sendSlaves(sendBuff,determType,vertexData , link);
        }
        
        if (containsEdges) 
        {
          hedge_STI * determType = 0;
          sendSlaves(sendBuff,determType, edgeData , link);
        }
        
        if (containsFaces) 
        {
          hface_STI * determType = 0;
          std::pair< IteratorSTI < hface_STI > * , IteratorSTI < hface_STI > * >
            iterpair = borderIteratorTT(determType , link );
         
          // pack all faces that we are master on 
          sendFaces( sendBuff, iterpair.first  , faceData ); 
          // pack also all faces that we are not master on 
          sendFaces( sendBuff, iterpair.second , faceData ); 

          delete iterpair.first;
          delete iterpair.second;
        } 
      }
     
      /////////////////////////////////////////////////////
      // den anderen Partitionen die Slave-Daten senden
      /////////////////////////////////////////////////////
      osvec = mpAccess ().exchange (osvec);

      // now get all sended data and store on master item in local buffers
      for (int link = 0; link < nl; ++link) 
      { 
        ObjectStream & recvBuff = osvec[link];
        
        if (containsVertices) 
        {
          vertex_STI * determType = 0;
          unpackOnMaster(recvBuff,vertexCommMap,determType,vertexData,nl,link);
        }

        if (containsEdges) 
        {
          hedge_STI * determType = 0;
          unpackOnMaster(recvBuff,edgeCommMap,determType,edgeData,nl,link);
        }

        if (containsFaces) 
        {
          hface_STI * determType = 0;
          std::pair< IteratorSTI < hface_STI > * , IteratorSTI < hface_STI > * >
            iterpair = borderIteratorTT( determType , link );

          // first unpack slave data 
          unpackFaces(recvBuff,iterpair.second,faceData);
          // then unpack all master data 
          unpackFaces(recvBuff,iterpair.first ,faceData);

          delete iterpair.first;
          delete iterpair.second;
        }
      }
    }

    // now get all data from the local buffer of the master 
    // and send this data to the slaves (only for vertices and edges)
    if( haveVerticesOrEdges )
    {
      for (int link = 0; link < nl; ++link ) 
      {
        ObjectStream & sendBuff = osvec[link];
        sendBuff.clear();
        
        // write Number of my links 
        sendBuff.writeObject(nl); 

        if (containsVertices) 
        {
          vertex_STI * determType = 0;
          sendMaster(sendBuff,vertexCommMap,determType,vertexData,nl, link );
        }
        
        if (containsEdges) 
        {
          hedge_STI * determType = 0;
          sendMaster(sendBuff,edgeCommMap,determType,edgeData,nl, link );
        }
      }

      // clear buffers to save memory 
      vertexCommMap.clear();
      edgeCommMap.clear();
     
      ///////////////////////////////////////////////////
      // exchange all gathered data 
      ///////////////////////////////////////////////////
      osvec = mpAccess ().exchange (osvec);
     
      // now unpack all data on slave items 
      for (int link = 0; link < nl; ++link) 
      { 
        ObjectStream & recvBuff = osvec[link];
        
        int nOtherlinks;
        recvBuff.readObject(nOtherlinks); // read number of links 

        if (containsVertices) 
        {
          vertex_STI * determType = 0;
          unpackOnSlaves(recvBuff,determType,vertexData, nOtherlinks, link );
        }
        
        if (containsEdges) 
        {
          hedge_STI * determType = 0;
          unpackOnSlaves(recvBuff,determType, edgeData, nOtherlinks, link );
        }
      }

    } // end second loop over vertices and edges 

    return;
  }


  // pack element data to stream 
  void GitterDunePll::sendInteriorGhostElementData (
      ObjectStream & sendBuff, 
      IteratorSTI < hface_STI > * iter , 
      GatherScatterType & elementData)
  {
#ifdef ALUGRIDDEBUG
    const bool containsElements = elementData.contains(3,0);
    alugrid_assert ( containsElements );
#endif
    const int transmit = 1;
    for (iter->first (); ! iter->done (); iter->next ()) 
    {
      hface_STI & face = iter->item(); 
      
      // check ghost leaf 
      std::pair< ElementPllXIF_t *, int > inner = face.accessInnerPllX ();

      if ( elementData.containsInterior(face, *(inner.first) ) ) 
      { 
        sendBuff.writeObject(transmit);

        // first interior elements are packed 
        inner.first->writeDynamicState (sendBuff , elementData);
      }     
      else 
      {
        sendBuff.writeObject(noData);
      }
    }

    return;
  }

  // unpack all data from stream 
  void GitterDunePll::unpackInteriorGhostElementData (
      ObjectStream & recvBuff, 
      IteratorSTI < hface_STI > * iter , 
      GatherScatterType & elementData )
  {
#ifdef ALUGRIDDEBUG
    const bool containsElements = elementData.contains(3,0);
    alugrid_assert ( containsElements );
#endif
    
    for (iter->first (); ! iter->done (); iter->next ()) 
    {
      int hasdata = 0;        
      recvBuff.readObject(hasdata);

      if( hasdata ) 
      {
        std::pair< ElementPllXIF_t *, int > p = iter->item ().accessOuterPllX ();
        p.first->readDynamicState ( recvBuff , elementData);
      }
    }
    return;
  }

  // pack all data to stream 
  void GitterDunePll::sendInteriorGhostAllData (
      ObjectStream & sendBuff, 
      IteratorSTI < hface_STI > * iter , 
      GatherScatterType & vertexData , 
      GatherScatterType & edgeData,  
      GatherScatterType & faceData, 
      GatherScatterType & elementData ,
      const bool packInterior ,
      const bool packGhosts )
  {
    const bool containsVertices = vertexData.contains(3,3);
    const bool containsEdges    = edgeData.contains(3,2);
    const bool containsFaces    = faceData.contains(3,1);
    
    const bool haveHigherCodimData = containsVertices || 
      containsEdges ||  
      containsFaces;

    const bool containsElements = elementData.contains(3,0);
    
    std::pair< ElementPllXIF_t *, int > bnd( ( ElementPllXIF_t * ) 0 , -1);

    // temporary object buffer  
    for (iter->first (); ! iter->done (); iter->next ()) 
    {
      hface_STI & face = iter->item(); 
      
      // check ghost leaf 
      std::pair< ElementPllXIF_t *, int > inner = face.accessInnerPllX ();

      int interiorLeaf = 0;
      int ghostLeaf = 0;

      if(packInterior) 
      {
        interiorLeaf = (elementData.containsInterior(face, *(inner.first) )) ? 1 : 0;
      }

      if(packGhosts)
      {
        bnd = face.accessOuterPllX ();
        ghostLeaf = (elementData.containsGhost(face , *(bnd.first))) ? 2 : 0;
      }

      const int transmit = interiorLeaf + ghostLeaf;
      // transmit = 1 interior, transmit = 2 ghost, transmit = 3 both 
      // if at least one of this possibilities is true then send data
      if ( transmit > 0 ) 
      { 
        sendBuff.writeObject(transmit);

        // first interior elements are packed 
        if( interiorLeaf  > 0 )
        {
          if( haveHigherCodimData )
          {
            if (containsVertices) 
              inner.first->VertexData2os(sendBuff, vertexData, inner.second );
            if (containsEdges)    
              inner.first->EdgeData2os  (sendBuff, edgeData  , inner.second );
            if (containsFaces)    
              inner.first->FaceData2os  (sendBuff, faceData  , inner.second );
          }

          if (containsElements) 
            inner.first->writeDynamicState (sendBuff , elementData);
        }

        // then ghost elements 
        if( ghostLeaf > 0 ) 
        {
          if( haveHigherCodimData )
          {
            // get std::pair< ghost, local face num > 
            Gitter::ghostpair_STI gpair = bnd.first->getGhost();
            alugrid_assert ( gpair.first );

            if (containsVertices) 
              gpair.first->VertexData2os( sendBuff , vertexData, gpair.second );
            if (containsEdges)    
              gpair.first->EdgeData2os  ( sendBuff , edgeData, gpair.second );
            if (containsFaces)    
              gpair.first->FaceData2os  ( sendBuff , faceData, gpair.second );
          }
          
          if( containsElements ) 
          {
            alugrid_assert ( bnd.first );
            bnd.first->writeDynamicState(sendBuff, elementData );
          }

          // reset bnd pointer 
          bnd.first = 0;
        }
      }     
      else 
      {
        sendBuff.writeObject(noData);
      }
    }

    return;
  }

  // unpack all data from stream 
  void GitterDunePll::unpackInteriorGhostAllData (
      ObjectStream & recvBuff, 
      IteratorSTI < hface_STI > * iter , 
      GatherScatterType & vertexData , 
      GatherScatterType & edgeData,  
      GatherScatterType & faceData, 
      GatherScatterType & elementData )
  {
    const bool containsVertices = vertexData.contains(3,3);
    const bool containsEdges    = edgeData.contains(3,2);
    const bool containsFaces    = faceData.contains(3,1);
    const bool containsElements = elementData.contains(3,0);
    

    const bool haveHigherCodimData = containsVertices || 
                                     containsEdges ||  
                                     containsFaces;

    
    for (iter->first (); ! iter->done (); iter->next ()) 
    {
      int hasdata;        
      recvBuff.readObject(hasdata);

      if(hasdata != noData)
      {
        // interiorLeaf is true, if on other side ghostLeaf has been packed
        const bool interiorLeaf = (hasdata > 1);

        // ghostLeaf is true if on other side interior has been packed
        const bool ghostLeaf    = (hasdata == 1) || (hasdata == 3);

        // get face 
        hface_STI & face = iter->item ();
        // first unpack ghosts 
        if( ghostLeaf ) 
        {
          std::pair< ElementPllXIF_t *, int > p = face.accessOuterPllX ();

          // get std::pair< ghost, local face num > 
          Gitter::ghostpair_STI gpair = p.first->getGhost();
          alugrid_assert ( gpair.first );
        
          if( haveHigherCodimData )
          {
            if (containsVertices) 
              gpair.first->os2VertexData( recvBuff , vertexData, gpair.second );
            if (containsEdges)    
              gpair.first->os2EdgeData  ( recvBuff , edgeData, gpair.second );
            if (containsFaces)    
              gpair.first->os2FaceData  ( recvBuff , faceData, gpair.second );
          }

          if (containsElements) 
            p.first->readDynamicState ( recvBuff , elementData);
        }

        // then unpack interior 
        if( interiorLeaf )
        {
          std::pair< ElementPllXIF_t *, int > pll = face.accessInnerPllX ();
          std::pair< Gitter::helement_STI* , Gitter::hbndseg_STI * > 
            p ( (Gitter::helement_STI *) 0, (Gitter::hbndseg_STI *) 0);

          pll.first->getAttachedElement( p );
          alugrid_assert ( p.first );

          if( haveHigherCodimData )
          {
            if (containsVertices) 
              p.first->os2VertexData( recvBuff , vertexData, pll.second );
            if (containsEdges)    
              p.first->os2EdgeData( recvBuff , edgeData, pll.second );
            if (containsFaces)    
              p.first->os2FaceData( recvBuff , faceData, pll.second );
          }

          if (containsElements) 
            elementData.recvData( recvBuff , *(p.first) );
        }
      }
    }
    return;
  }


  /////////////////////////////////////////////////////
  //
  //  interior to ghost communication 
  //
  /////////////////////////////////////////////////////
  class PackUnpackInteriorGhostData : public MpAccessLocal::NonBlockingExchange::DataHandleIF
  {
    typedef Gitter::hface_STI hface_STI ;

    GitterDunePll& _gitter ;
    GatherScatterType& _vertexData ;
    GatherScatterType& _edgeData ;
    GatherScatterType& _faceData ;
    GatherScatterType& _elementData ;

    const bool _haveHigherCodimData;

    const GitterDunePll::CommunicationType _commType;

    PackUnpackInteriorGhostData( const PackUnpackInteriorGhostData& );
  public:
    PackUnpackInteriorGhostData( GitterDunePll& gitter,
                                 GatherScatterType& vertexData, 
                                 GatherScatterType& edgeData, 
                                 GatherScatterType& faceData, 
                                 GatherScatterType& elementData, 
                                 const GitterDunePll::CommunicationType commType )
      : _gitter( gitter ),
        _vertexData ( vertexData ),
        _edgeData   ( edgeData ),
        _faceData   ( faceData ),
        _elementData( elementData ),
        _haveHigherCodimData( vertexData.contains(3,3) || edgeData.contains(3,2) || faceData.contains(3,1) ),
        _commType( commType )
    {
    }

    // return true if at least one codimension is contained 
    bool containsSomeThing () const { return (_haveHigherCodimData || _elementData.contains(3,0) ); }

    void pack( const int link, ObjectStream& sendBuff ) 
    {
      const bool packInterior = (_commType == GitterDunePll::All_All_Comm) || 
                                (_commType == GitterDunePll::Interior_Ghost_Comm);
      
      const bool packGhosts   = (_commType == GitterDunePll::All_All_Comm) || 
                                (_commType == GitterDunePll::Ghost_Interior_Comm);

      // clear buffer 
      sendBuff.clear();
      
      const hface_STI * determType = 0; // only for type determination 
      std::pair< IteratorSTI < hface_STI > * , IteratorSTI < hface_STI > * > 
        iterpair = _gitter.borderIteratorTT( determType , link );

      if( _haveHigherCodimData || packGhosts )
      {
        // write all data belong to interior of master faces 
        _gitter.sendInteriorGhostAllData( sendBuff, iterpair.first , 
                                          _vertexData, _edgeData,
                                          _faceData, _elementData, 
                                          packInterior , packGhosts );
      
        // write all data belong to interior of slave faces 
        _gitter.sendInteriorGhostAllData( sendBuff, iterpair.second , 
                                          _vertexData, _edgeData,
                                          _faceData, _elementData ,
                                          packInterior , packGhosts );
      }
      else 
      {
        // write all data belong to interior of master faces 
        _gitter.sendInteriorGhostElementData( sendBuff, iterpair.first, _elementData);
      
        // write all data belong to interior of slave faces 
        _gitter.sendInteriorGhostElementData( sendBuff, iterpair.second, _elementData);
      }

      delete iterpair.first; 
      delete iterpair.second; 
    }

    void unpack( const int link, ObjectStream& recvBuff ) 
    {
      const hface_STI * determType = 0; // only for type determination 
      std::pair< IteratorSTI < hface_STI > * , IteratorSTI < hface_STI > * > 
        iterpair = _gitter.borderIteratorTT( determType , link );

      const bool packGhosts = (_commType == GitterDunePll::All_All_Comm) || 
                              (_commType == GitterDunePll::Ghost_Interior_Comm);

      if( _haveHigherCodimData || packGhosts )
      {
        // first unpack slave data, because this has been pack from master
        // first , see above 
        _gitter.unpackInteriorGhostAllData( recvBuff, iterpair.second , 
                                            _vertexData, _edgeData,
                                            _faceData, _elementData );
                       
        // now unpack data sended from slaves to master 
        _gitter.unpackInteriorGhostAllData( recvBuff, iterpair.first , 
                                            _vertexData, _edgeData,
                                            _faceData, _elementData );
      }
      else 
      {
        // first unpack slave data, because this has been pack from master
        // first , see above 
        _gitter.unpackInteriorGhostElementData( recvBuff, iterpair.second, _elementData );
       
        // now unpack data sended from slaves to master 
        _gitter.unpackInteriorGhostElementData( recvBuff, iterpair.first, _elementData );
      }

      delete iterpair.first;
      delete iterpair.second;
    }
  };

  void GitterDunePll::doInteriorGhostComm( 
    std::vector< ObjectStream > & osvec ,
    GatherScatterType & vertexData , 
    GatherScatterType & edgeData,  
    GatherScatterType & faceData, 
    GatherScatterType & elementData ,
    const CommunicationType commType )
  {
    // create data handle 
    PackUnpackInteriorGhostData data( *this, vertexData, edgeData, faceData, elementData, commType );

    if(! data.containsSomeThing() ) 
    {
      std::cerr << "WARNING: communication called with empty data set, all contains methods returned false! \n";
      return;
    }

    ///////////////////////////////////////////
    // exchange data 
    ///////////////////////////////////////////
    mpAccess ().exchange ( data );     
  } 

  ////////////////////////////////////////////////////////
  //
  // communicate data
  // 
  ////////////////////////////////////////////////////////
  void GitterDunePll::doCommunication ( 
               GatherScatterType & vertexData , 
               GatherScatterType & edgeData,  
               GatherScatterType & faceData ,
               GatherScatterType & elementData, 
               const CommunicationType commType ) 
  {
    const int nl = mpAccess ().nlinks ();

    const bool containsVertices = vertexData.contains(3,3);
    const bool containsEdges    = edgeData.contains(3,2);
    const bool containsFaces    = faceData.contains(3,1);
    const bool containsElements = elementData.contains(3,0);

    const bool haveHigherCodimData = containsVertices || 
      containsEdges ||  
      containsFaces;

    const bool containsSomeThing = containsElements || haveHigherCodimData;

    if(!containsSomeThing) 
    {
      std::cerr << "WARNING: communication called with empty data set, all contains methods returned false! \n";
      return;
    }
     
    alugrid_assert ((debugOption (5) && containsVertices) ? (std::cout << "**INFO GitterDunePll::doCommunication (): (containsVertices)=true " << std::endl, 1) : 1);
    alugrid_assert ((debugOption (5) && containsEdges)    ? (std::cout << "**INFO GitterDunePll::doCommunication (): (containsEdges)=true " << std::endl, 1) : 1);
    alugrid_assert ((debugOption (5) && containsFaces)    ? (std::cout << "**INFO GitterDunePll::doCommunication (): (containsFaces)=true " << std::endl, 1) : 1);
    alugrid_assert ((debugOption (5) && containsElements) ? (std::cout << "**INFO GitterDunePll::doCommunication (): (containsElements)=true " << std::endl, 1) : 1);
     
    // create vector of message buffers 
    // this vector is created here, that the buffer is allocated only once 
    std::vector< ObjectStream > vec (nl);
   
    // if data on entities of higer codim exists
    // then communication if more complicated 
    if ( haveHigherCodimData )
    {
      doBorderBorderComm( vec, vertexData, edgeData, faceData );
    }

    // this communication only makes sense if ghost cells are present 
    const bool ghostCellsAvailable = ghostCellsEnabled ();

    if( (commType != Border_Border_Comm) && ghostCellsAvailable ) // otherwise only border border 
    {
      doInteriorGhostComm( vec, vertexData, edgeData, faceData, elementData , commType ); 
    }

    return;
  }

  // border border comm 
  void GitterDunePll::borderBorderCommunication ( 
               GatherScatterType & vertexData , 
               GatherScatterType & edgeData,  
               GatherScatterType & faceData ,
               GatherScatterType & elementData )
  {
    doCommunication(vertexData,edgeData,faceData,elementData,Border_Border_Comm);
  }

  // interior ghost comm 
  void GitterDunePll::interiorGhostCommunication ( 
               GatherScatterType & vertexData , 
               GatherScatterType & edgeData,  
               GatherScatterType & faceData ,
               GatherScatterType & elementData )
  {
    doCommunication(vertexData,edgeData,faceData,elementData,Interior_Ghost_Comm);
  }

  // ghost to interior comm 
  void GitterDunePll::ghostInteriorCommunication ( 
               GatherScatterType & vertexData , 
               GatherScatterType & edgeData,  
               GatherScatterType & faceData ,
               GatherScatterType & elementData )
  {
    doCommunication(vertexData,edgeData,faceData,elementData,Ghost_Interior_Comm);
  }

  // all all comm 
  void GitterDunePll::allAllCommunication ( 
               GatherScatterType & vertexData , 
               GatherScatterType & edgeData,  
               GatherScatterType & faceData ,
               GatherScatterType & elementData )
  {
    doCommunication(vertexData,edgeData,faceData,elementData,All_All_Comm);
  }

  ////////////////////////////////////////////////////////////////////////////
  ////////////////////////////////////////////////////////////////////////////

  // rebuild ghost cells 
  void GitterDunePll::rebuildGhostCells() 
  {
    // only do this if ghost cells are enabled 
    if( ! ghostCellsEnabled() ) return;

    // compute the graph vertices new, even if already computed 
    GitterPll :: computeGraphVertexIndices  ();

    try 
    {
      // get number of links
      const int nl = mpAccess ().nlinks ();

      std::vector< ObjectStream > osv (nl);
      
      const hface_STI* determType = 0;

      // pack all elements neighbouring to internal boundary 
      // as ghost elements 
      {
        for (int link = 0; link < nl; ++link ) 
        {
          std::pair< IteratorSTI < hface_STI > *, IteratorSTI < hface_STI > *> 
                w = levelBorderIteratorTT (determType, link, 0 );
          
          ObjectStream & os = osv[link];
          
          {
            IteratorSTI < hface_STI > & inner = *w.first;
          
            for ( inner.first (); ! inner.done (); inner.next ()) 
            {
              std::pair< ElementPllXIF_t *, int > p = inner.item ().accessInnerPllX ();
              p.first->packAsGhost(os, p.second);
            }
          }

          {
            IteratorSTI < hface_STI > & outer = *w.second;
            for (outer.first (); ! outer.done (); outer.next ()) 
            {
              std::pair< ElementPllXIF_t *, int > p = outer.item ().accessInnerPllX ();
              p.first->packAsGhost(os, p.second);
            }
          }

          delete w.first;
          delete w.second;
        } 
      }
      
      // exchange gathered data 
      osv = mpAccess ().exchange (osv);
      
      // unpack all data on internal boundary and create 
      // ghost cells 
      {
        for (int link = 0; link < nl; ++link ) 
        {
          std::pair< IteratorSTI < hface_STI > *, IteratorSTI < hface_STI > *> 
                w = levelBorderIteratorTT (determType, link, 0 );
          
          ObjectStream & os = osv[link];

          {
            IteratorSTI < hface_STI > & outer = *w.second;
            for (outer.first (); ! outer.done (); outer.next ()) 
            {
              std::pair< ElementPllXIF_t *, int > p = outer.item ().accessOuterPllX ();
              p.first->insertGhostCell(os, p.second);
            }
          }
          {
            IteratorSTI < hface_STI > & inner = *w.first;
            for (inner.first (); ! inner.done (); inner.next ()) 
            {
              std::pair< ElementPllXIF_t *, int > p = inner.item ().accessOuterPllX ();
              p.first->insertGhostCell(os, p.second);
            }
          }

          delete w.first;
          delete w.second;
        } 
      }
    } 
    catch (Parallel:: AccessPllException) 
    {
      std::cerr << "  FEHLER Parallel::AccessPllException entstanden in: " << __FILE__ << " " << __LINE__ << std::endl;
    }

    return;
  }

  void GitterDunePll::checkGhostIndices() 
  {
    // only do this if ghost cells are enabled 
    if( ! ghostCellsEnabled() ) return;

    // get number of links 
    const int nl = mpAccess ().nlinks ();
    
    const hface_STI* determType = 0;
    {
      // for all links check all ghost elements 
      for (int link = 0; link < nl; ++link ) 
      {
        std::pair< IteratorSTI < hface_STI > *, IteratorSTI < hface_STI > *> 
              w = levelBorderIteratorTT (determType, link , 0);
        
        {
          IteratorSTI < hface_STI > & outer = *w.second;
          for (outer.first (); ! outer.done (); outer.next ()) 
          {
            std::pair< ElementPllXIF_t *, int > p = outer.item ().accessOuterPllX ();

            // get std::pair< ghost, local face num > 
            Gitter::ghostpair_STI gpair = p.first->getGhost();
            alugrid_assert ( gpair.first );

            gpair.first->resetGhostIndices();
          }
        }

        {
          IteratorSTI < hface_STI > & inner = *w.first;
          for (inner.first (); ! inner.done (); inner.next ()) 
          {
            std::pair< ElementPllXIF_t *, int > p = inner.item ().accessOuterPllX ();

            // get std::pair< ghost, local face num > 
            Gitter::ghostpair_STI gpair = p.first->getGhost();
            alugrid_assert ( gpair.first );

            gpair.first->resetGhostIndices();
          }
        }

        // free interators 
        delete w.first;
        delete w.second;
      } 
    } 
  }

  void GitterDunePll::duneBackup ( const char *filename ) 
  {
    // backup grid, same as in serial case 
    GitterDuneBasis::duneBackup( filename );
  }

  // wird von Dune verwendet 
  void GitterDunePll::restore ( std::istream &in ) 
  {
    typedef Gitter::Geometric::BuilderIF BuilderIF;
    alugrid_assert (debugOption (20) ? (std::cout << "**INFO GitterDunePll::restore (istream & = " << in << ") " << std::endl, 1) : 1);
    {
      AccessIterator < hedge_STI >::Handle ew (container ());
      for (ew.first (); !ew.done (); ew.next ()) ew.item ().restore (in);
    }
    {
      AccessIterator < hface_STI >:: Handle fw(container());
      for ( fw.first(); !fw.done (); fw.next()) fw.item().restore (in);
    }
    {
      AccessIterator < helement_STI >:: Handle ew(container());
      for ( ew.first(); !ew.done(); ew.next()) ew.item().restore (in);
    }

    // restore indices before ghosts are created 
    // otherwise indices of ghost will be wrong 
    this->restoreIndices (in);
   
#ifdef ALUGRIDDEBUG 
    const int maxIndexBefore = this->indexManager(BuilderIF::IM_Elements).getMaxIndex();
#endif

    // set ghost indices new for level 0 ghosts 
    checkGhostIndices ();

    // now restore faces and by this ghosts 
    // will be refined 
    {
      AccessIterator < hbndseg_STI >::Handle bw (container ());
      for (bw.first (); ! bw.done (); bw.next ()) bw.item ().restoreFollowFace ();
    }
    
    // max index should not be largen than before
    alugrid_assert ( (this->indexManager(BuilderIF::IM_Elements).getMaxIndex() != maxIndexBefore) ?
        (std::cout << maxIndexBefore << " vor | nach " << this->indexManager(BuilderIF::IM_Elements).getMaxIndex() << "\n",0) : 1);

  }

  void GitterDunePll::duneRestore ( const char *fileName )
  {
    alugrid_assert (debugOption (20) ? 
        (std::cout << "**INFO GitterDuneBasis::duneRestore (const char * = \""
              << fileName << "\") " << std::endl, 1) : 1);
                   
    std::ifstream in( fileName );
    if (!in) {
      std::cerr << "**WARNUNG (IGNORIERT) GitterDunePll::";
      std::cerr <<" duneRestore (const char *, double & ) Fehler beim \"Offnen von < "
           << (fileName ? fileName : "null") << " > " << std::endl;
    } 
    else 
    {
      duneRestore( in );
    }
  }

  void GitterDunePll::tovtk ( const std::string &fn ) 
  {
    const int myrank = mpAccess ().myrank ();
    const int nProc = mpAccess ().psize ();

    std::ostringstream ss;
    ss << "p" << myrank << "-" << fn;

    // openfile
    std::ofstream vtuFile;
    Gitter::tovtk( ss.str() );

    if( myrank == 0 )
    {
      std::ostringstream pllss;
      if( fn.substr(fn.find_last_of(".") + 1) == "vtu" )
        pllss << fn.substr(0,fn.find_last_of(".") + 1) << "pvtu";
      else
        pllss << fn << ".pvtu";
      std::ofstream pvtuFile;
      pvtuFile.open( pllss.str().c_str() );
    
      pvtuFile << "<?xml version=\"1.0\"?>" << std::endl;
      pvtuFile << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">" << std::endl;
      pvtuFile << "  <PUnstructuredGrid GhostLevel=\"0\">" << std::endl;
      pvtuFile << "    <PCellData Scalars=\"cell-id\">" << std::endl;
      pvtuFile << "      <PDataArray type=\"Float32\" Name=\"cell-id\" />" << std::endl;
      pvtuFile << "    </PCellData>" << std::endl;
      pvtuFile << "    <PPoints>" << std::endl;
      pvtuFile << "      <PDataArray type=\"Float32\" NumberOfComponents=\"3\" />" << std::endl;
      pvtuFile << "    </PPoints>" << std::endl;
      for( int p = 0; p < nProc; ++p )
        pvtuFile << "    <Piece Source=\"p" << p << "-" << fn << "\" />" << std::endl;
      pvtuFile << "  </PUnstructuredGrid>" << std::endl;
      pvtuFile << "</VTKFile>" << std::endl;

      pvtuFile.close();

      std::cout << "data written to " << pllss.str() << std::endl;
    }
  }
} // namespace ALUGrid
