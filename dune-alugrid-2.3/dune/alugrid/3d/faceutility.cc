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
#include "faceutility.hh"

namespace Dune
{

  // Implementation of ALU3dGridSurfaceMappingFactory
  // ------------------------------------------------

  template< class Comm >
  typename ALU3dGridSurfaceMappingFactory< tetra, Comm >::SurfaceMappingType *
  ALU3dGridSurfaceMappingFactory< tetra, Comm >::buildSurfaceMapping ( const CoordinateType &coords ) const
  {
    return new SurfaceMappingType( fieldVector2alu3d_ctype(coords[0]) , 
                                   fieldVector2alu3d_ctype(coords[1]) , 
                                   fieldVector2alu3d_ctype(coords[2]) );
  }


  template< class Comm >
  typename ALU3dGridSurfaceMappingFactory< hexa, Comm >::SurfaceMappingType *
  ALU3dGridSurfaceMappingFactory< hexa, Comm >::buildSurfaceMapping ( const CoordinateType &coords ) const
  {
    return new SurfaceMappingType( coords[0], coords[1], coords[2], coords[3] );
  }


  template< class Comm >
  typename ALU3dGridSurfaceMappingFactory< tetra, Comm >::SurfaceMappingType *
  ALU3dGridSurfaceMappingFactory< tetra, Comm >::buildSurfaceMapping ( const GEOFaceType &face ) const
  {
    return new SurfaceMappingType( face.myvertex(0)->Point(), face.myvertex(1)->Point(), face.myvertex(2)->Point() );
  }


  template< class Comm >
  typename ALU3dGridSurfaceMappingFactory< hexa, Comm >::SurfaceMappingType *
  ALU3dGridSurfaceMappingFactory< hexa, Comm >::buildSurfaceMapping ( const GEOFaceType &face ) const
  {
    typedef FaceTopologyMapping< hexa > FaceTopo;
    // this is the new implementation using FieldVector 
    // see mappings.hh 
    // we have to swap the vertices, because 
    // local face numbering in Dune is different to ALUGrid (see topology.cc)
    return new SurfaceMappingType(
        face.myvertex( FaceTopo::dune2aluVertex(0) )->Point(),
        face.myvertex( FaceTopo::dune2aluVertex(1) )->Point(),
        face.myvertex( FaceTopo::dune2aluVertex(2) )->Point(),
        face.myvertex( FaceTopo::dune2aluVertex(3) )->Point() );
  }



  // Helper Functions
  // ----------------

  template< int m, int n >
  inline void
  alu3dMap2World ( const ALU3DSPACE LinearSurfaceMapping &mapping,
                   const FieldVector< alu3d_ctype, m > &x,
                   FieldVector< alu3d_ctype, n > &y )
  {
    mapping.map2world( fieldVector2alu3d_ctype( x ), fieldVector2alu3d_ctype( y ) );
  }

  template< int m, int n >
  inline void
  alu3dMap2World ( const BilinearSurfaceMapping &mapping,
                   const FieldVector< alu3d_ctype, m > &x,
                   FieldVector< alu3d_ctype, n > &y )
  {
    mapping.map2world( x, y );
  }



  //- class ALU3dGridGeometricFaceInfoBase
  template< ALU3dGridElementType type, class Comm >
  void ALU3dGridGeometricFaceInfoBase< type, Comm >
    ::referenceElementCoordinatesUnrefined ( SideIdentifier side, CoordinateType &result ) const
  {
    // get the parent's face coordinates on the reference element (Dune reference element)
    CoordinateType cornerCoords;
    referenceElementCoordinatesRefined ( side, cornerCoords );

    typename Base::SurfaceMappingType *referenceElementMapping = Base::buildSurfaceMapping( cornerCoords );

    NonConformingMappingType faceMapper( connector_.face().parentRule(), connector_.face().nChild() );

    const ReferenceFaceType& refFace = getReferenceFace();
    // do the mappings
    const int numCorners = refFace.size( 2 );
    for( int i = 0; i < numCorners; ++i )
    {
      const FieldVector< alu3d_ctype, 2 > &childLocal = refFace.position( i, 2 );
      alu3dMap2World( *referenceElementMapping, faceMapper.child2parent( childLocal ), result[ i ] );
    }

    delete referenceElementMapping;
  }



  // Explicit Template Instatiation
  // ------------------------------

  template struct ALU3dGridSurfaceMappingFactory< tetra, ALUGridNoComm >;
  template struct ALU3dGridSurfaceMappingFactory< hexa, ALUGridNoComm >;

  template class ALU3dGridGeometricFaceInfoBase< tetra, ALUGridNoComm >;
  template class ALU3dGridGeometricFaceInfoBase< hexa, ALUGridNoComm >;

  template struct ALU3dGridSurfaceMappingFactory< tetra, ALUGridMPIComm >;
  template struct ALU3dGridSurfaceMappingFactory< hexa, ALUGridMPIComm >;

  template class ALU3dGridGeometricFaceInfoBase< tetra, ALUGridMPIComm >;
  template class ALU3dGridGeometricFaceInfoBase< hexa, ALUGridMPIComm >;

} // end namespace Dune
