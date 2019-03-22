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

#ifndef __HEADER__VTX_BTREE
#define __HEADER__VTX_BTREE

#include <algorithm>

namespace ALU2DGrid
{

  template < int N > class Vertex;
  template < int N, int NV > class Thinelement;
  template < int N, int NV > class Element;

  // ============================================================
  // Klasse Vtx_btree: Binaerer Baum von Vertex-Instanzen.
  // Ordnet die Knoten nach ihrem Abstand zu einem Referenzknoten,
  // der dem Konstruktor uebergeben wird.
  // ============================================================

  template< int N, int NV >
  class Vtx_btree
  {
  public:
    typedef Vertex < N > vertex_t;
    typedef Thinelement < N, NV > thinelement_t;
    typedef Element < N, NV > element_t;

    // Implementation eines Knotens des binaeren Baumes
    struct Node {
      vertex_t *vtx;
      thinelement_t *lnb;
      thinelement_t *rnb;
      Node* next;
      Node* prev;
      int _lidx,_ridx;
      Node(vertex_t *invtx, thinelement_t *plnb, thinelement_t *prnb)
      : vtx(invtx), lnb(plnb), rnb(prnb), next(0), prev(0)
      {alugrid_assert (invtx);}

      ~Node() {
        if(next) delete next;
        if(prev) delete prev;
      }

      Node* leftNode() {
        return prev;
      }
      Node* rightNode() {
        return next;
      }
      thinelement_t *leftElement() {
        return lnb;
      }
      thinelement_t *rightElement() {
        return rnb;
      }
      vertex_t *vertex() {
        return vtx;
      }

      // return depth of this tree
      int deepestLevel ( int prevLvl = 0 ) const
      {
        const int left  = (prev ? prev->deepestLevel() : 0);
        const int right = (next ? next->deepestLevel() : 0);
        return std::max(left,right) + prevLvl + 1;
      }

      // return number of hanging nodes in this tree
      int count() const {
        return 1 + (next ? next->count() : 0) + (prev ? prev->count() : 0);
      }

      int remove(vertex_t *pvtx);

      void nbconnect(int,thinelement_t *,int);
    }* head;

 public:
    vertex_t *rvtx;
    thinelement_t *lnb;
    thinelement_t *rnb;

    void insertNode(Node* node, Node* newNode);

    double dist(Vertex < N > *invtx);

    Vtx_btree* left() const;

    Vtx_btree* right() const;

  public:
    Vtx_btree(vertex_t *invtx, thinelement_t *plnb, thinelement_t *prnb)
      : head(0), rvtx(invtx), lnb(plnb), rnb(prnb) {
      alugrid_assert (rvtx);
      alugrid_assert (plnb);
      alugrid_assert (prnb);
    }
    
    ~Vtx_btree() {
      if( head ) delete head;
    }

    vertex_t *getHead() { return head->vtx; }

    thinelement_t *getlnb() { return head->lnb; }
    thinelement_t *getrnb() { return head->rnb; }

    void insert(vertex_t *invtx, thinelement_t *plnb, thinelement_t *prnb);
    
    void splitTree(Vtx_btree*& inleft, Vtx_btree*& inright);
    
    void merge(Vtx_btree* inleft, Vtx_btree* inright);

    void nbconnect ( int, thinelement_t *, int );

    int deepestLevel ()
    {
      return (head ? head->deepestLevel() : 0);
    }

    int count () const { return head->count(); }
    
    bool remove(Vertex < N > *vtx) {
      alugrid_assert (head->prev || head->next);
      return (head->remove(vtx)==1);
    }
  };

} // namespace ALU2DGrid

#endif // VTXBTREE_H
