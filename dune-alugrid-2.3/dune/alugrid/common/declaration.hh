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

#ifndef DUNE_ALUGRID_FORWARDDECLARATION
#define DUNE_ALUGRID_FORWARDDECLARATION

#define ALU3DGRID_PARALLEL HAVE_MPI

#include <dune/common/parallel/collectivecommunication.hh>
#if ALU3DGRID_PARALLEL
#include <dune/common/parallel/mpicollectivecommunication.hh>
#endif // #if ALU3DGRID_PARALLEL


namespace Dune
{

  //! \brief basic element types for ALUGrid
  enum ALUGridElementType { simplex, cube };
  //! \brief available refinement types for ALUGrid
  enum ALUGridRefinementType { conforming, nonconforming };

  //! \brief type of class for specialization of serial ALUGrid (No_Comm as communicator) 
  struct ALUGridNoComm 
  {
    No_Comm noComm_;
    ALUGridNoComm() : noComm_() {}
    ALUGridNoComm( const No_Comm& comm ) : noComm_( comm ) {}
    operator No_Comm () const { return noComm_; }
  };

  //! \brief type of class for specialization of parallel ALUGrid (MPI_Comm as communicator) 
  struct ALUGridMPIComm {
#if ALU3DGRID_PARALLEL
    MPI_Comm mpiComm_;
    ALUGridMPIComm() : mpiComm_( MPI_COMM_WORLD ) {}
    ALUGridMPIComm( MPI_Comm comm ) : mpiComm_( comm ) {}
    operator MPI_Comm () const { return mpiComm_; }
#endif
  } ;

/**
   \brief [<em> provides \ref Dune::Grid </em>]
   \brief grid with support for quadrilateral and hexahedral grid (template parameter cube) 
   and simplicial meshes (template parameter simplex) in 2d and 3d.
   @ingroup GridImplementations
   @ingroup ALUGrid

   The ALUGrid implements the Dune GridInterface for 2d quadrilateral and 3d hexahedral 
   as well as 2d triangular and  3d tetrahedral meshes. 
   This grid can be locally adapted (non-conforming and conforming bisection) 
   and used in parallel computations using dynamic load balancing. 

   @note
   (see ALUGrid homepage: http://www.mathematik.uni-freiburg.de/IAM/Research/alugrid/)

   \li Available Implementations 
        - quadrilateral and hexahedral elements only nonconforming refinement
          - Dune::ALUGrid< 2, 2, cube, nonconforming >  
          - Dune::ALUGrid< 2, 3, cube, nonconforming >
          - Dune::ALUGrid< 3, 3, cube, nonconforming >
        - simplicial elements and nonconforming refinement  
          - Dune::ALUGrid< 2, 2, simplex, nonconforming >  
          - Dune::ALUGrid< 2, 3, simplex, nonconforming >
          - Dune::ALUGrid< 3, 3, simplex, nonconforming >
        - simplicial elements and bisection refinement  
          - Dune::ALUGrid< 2, 2, simplex, conforming >  
          - Dune::ALUGrid< 2, 3, simplex, conforming >
          - Dune::ALUGrid< 3, 3, simplex, conforming > (work in progress)

   \note template parameter Comm defaults to MPI_Comm, if MPI is available, No_Comm  otherwise.

   For installation instructions see http://www.dune-project.org/external_libraries/install_alugrid.html .
*/
  template <int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, 
            class Comm = 
#if ALU3DGRID_PARALLEL
              ALUGridMPIComm
#else 
              ALUGridNoComm 
#endif
           >
  class ALUGrid;

  //- traits class for declaring base class for ALUGrid
  template <int dim, int dimw, ALUGridElementType elType, class Comm >
  struct ALUGridBaseGrid ;
}
#endif // #ifndef DUNE_ALUGRID_FORWARDDECLARATION
