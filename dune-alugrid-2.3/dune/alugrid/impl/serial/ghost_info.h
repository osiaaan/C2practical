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

// (c) Robert Kloefkorn 2007
#ifndef GHOSTINFO_H_INCLUDED
#define GHOSTINFO_H_INCLUDED

#include "gitter_sti.h"
#include "serialize.h"

namespace ALUGrid
{

  // interface class for macro ghost point
  class MacroGhostInfo
  : public MyAlloc 
  {
    public:
      virtual ~MacroGhostInfo () {}

      virtual const alucoord_t (& getPoint (int i) const )[3] = 0;
      virtual int nop () const = 0;
      virtual void inlineGhostElement(ObjectStream & os) const = 0; 
      virtual int orientation () const = 0;
  };

  typedef MacroGhostInfo MacroGhostInfo_STI;

  // little storage class for points and vertex numbers 
  // of transmitted macro elements to become ghosts 
  template <int points>
  class MacroGhostInfoStorage : public MacroGhostInfo
  {
  public:
    // number of points of element
    enum { noVx     = (points == 4) ? 8 : 4 };
    // number of all non-internal points 
    enum { noFaceVx = (points == 4) ? 4 : 1 };  

    static const int invalidFace = -11;
    
  protected:
    // coordiante of all non-internal points 
    alucoord_t _p[points][3];  // 24 or 96 bytes

    // vertex idents of all vertices of element 
    int _vx[noVx]; // 16 or 32 bytes 

    // vertex idents of all not internal vertices  
    int _vxface[noFaceVx]; // 4 or 16 bytes 

    // face number of internal face 
    int _fce; // 4 bytes 

    // do not allow copying
    MacroGhostInfoStorage(const MacroGhostInfoStorage & );

    MacroGhostInfoStorage() : _fce( invalidFace ) {}
  public:  
    // destructor 
    virtual ~MacroGhostInfoStorage () {}

    // return reference to _p
    const alucoord_t (& getPoints () const )[points][3] 
    { 
      alugrid_assert ( _fce != invalidFace );
      return _p; 
    }
    
    // return idents of ghost element  
    int (& vertices () )[noVx]
    {
      alugrid_assert ( _fce != invalidFace );
      return _vx; 
    }
    
    // return reference to vector with non-internal vertex idents 
    const int (& getOuterVertices () const )[noFaceVx]
    {
      alugrid_assert ( _fce != invalidFace );
      return _vxface;
    }

    // return local number of internal face 
    int internalFace () const 
    {
      alugrid_assert ( _fce != invalidFace );
      return _fce < 0 ? (-_fce) - 1 : _fce;
    }

    int orientation () const { return _fce < 0 ? 1 : 0; } 

    /////////////////////////////////////
    // interface of MacroGhostInfo_STI 
    ///////////////////////////////////// 
    virtual const alucoord_t (& getPoint (int i) const )[3] 
    {
      alugrid_assert ( _fce != invalidFace );
      alugrid_assert ( i>= 0 && i < points );
      return _p[i];
    }

    // return number of non-internal points 
    virtual int nop () const { return points; }

    // write internal data to stream 
    virtual void inlineGhostElement(ObjectStream&) const;

  protected:
    // read internal data from stream 
    void readData(ObjectStream&);
  };

  // macro ghost info for tetras 
  class MacroGhostInfoTetra
  : public MacroGhostInfoStorage< 1 >
  {
    enum { points = 1 };
    // do not copy 
    MacroGhostInfoTetra(const MacroGhostInfoTetra&);
  public:  
    // create storage by reading data from stream 
    explicit MacroGhostInfoTetra(ObjectStream& os) 
    { 
      this->readData(os); 
    }

    // contructor for tetras 
    MacroGhostInfoTetra(const Gitter:: Geometric :: tetra_GEO *, 
                        const int fce);
  };

  // macro ghost info for tetras 
  class MacroGhostInfoHexa
  : public MacroGhostInfoStorage< 4 >
  {
    enum { points = 4 };
    // no copying 
    MacroGhostInfoHexa(const MacroGhostInfoHexa&);
  public:  
    // create storage by reading data from stream 
    explicit MacroGhostInfoHexa(ObjectStream& os) 
    { 
      this->readData(os); 
    }

    // contructor for tetras 
    MacroGhostInfoHexa(const Gitter:: Geometric :: hexa_GEO * , 
                       const int fce);
  };

} // namespace ALUGrid

#endif // #ifndef GHOSTINFO_H_INCLUDED
