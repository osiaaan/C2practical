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

#ifndef DUNE_ALU_BNDPROJECTION_HH
#define DUNE_ALU_BNDPROJECTION_HH

namespace Dune {

  //! \brief ALUGrid boundary projection implementation 
  //!  DuneBndProjection has to fulfil the DuneBoundaryProjection interface 
  template <class GridImp, class ctype = double > 
  class ALUGridBoundaryProjection 
    : public GridImp :: ALUGridVertexProjectionType
  {
    typedef GridImp GridType;
    // type of double coordinate vector 
    typedef ctype coord_t[ GridType :: dimensionworld ];
  protected:

    //! reference to boundary projection implementation 
    const GridType& grid_;
  public: 
    //! type of boundary projection 
    typedef typename GridType :: DuneBoundaryProjectionType DuneBoundaryProjectionType;

    //! type of coordinate vector 
    typedef typename DuneBoundaryProjectionType :: CoordinateType CoordinateType;

    //! constructor storing reference to boundary projection implementation 
    ALUGridBoundaryProjection(const GridType& grid) 
      : grid_( grid ) 
    {
    }

    //! (old) method projection vertices defaults to segment 0
    int operator () (const coord_t &orig, 
                     coord_t &prj) const 
    {
      return this->operator()( orig, 0, prj);
    }

    //! projection operator  
    int operator () (const coord_t &orig, 
                     const int segmentIndex,
                     coord_t &prj) const 
    {
      // get boundary projection 
      const DuneBoundaryProjectionType* bndPrj = 
        grid_.boundaryProjection( segmentIndex );

      // if pointer is zero we do nothing, i.e. identity mapping
      if( bndPrj ) 
      {
        // call projection operator 
        reinterpret_cast<CoordinateType &> (* (&prj[0])) = 
          (*bndPrj)( reinterpret_cast<const CoordinateType &> (* (&orig[0])) );
      }

      // return 1 for success 
      return 1;
    }
  }; 

} // end namespace Dune 
#endif
