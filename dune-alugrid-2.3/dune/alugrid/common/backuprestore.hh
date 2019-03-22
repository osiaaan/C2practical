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

#ifndef DUNE_GRID_ALUGRID_BACKUPRESTORE_HH
#define DUNE_GRID_ALUGRID_BACKUPRESTORE_HH

//- system headers 
#include <fstream>

//- Dune headers 
#include <dune/common/exceptions.hh>
#include <dune/grid/common/backuprestore.hh>
#include <dune/alugrid/common/declaration.hh>

namespace Dune
{

  /** \copydoc Dune::BackupRestoreFacility */
  template< int dim, int dimworld, ALUGridElementType elType, ALUGridRefinementType refineType, class Comm >
  struct BackupRestoreFacility< ALUGrid< dim, dimworld, elType, refineType, Comm > >
  {
    // type of grid 
    typedef ALUGrid< dim, dimworld, elType, refineType, Comm > Grid;

    static std::string createFilename( const std::string &path, const std::string &fileprefix )
    {
      std::string filename( path );
      if( fileprefix.size() > 0 ) 
      {
        filename += "/" + fileprefix ;
      }
      else if( filename[ filename.size() - 1 ] == char('/') )
      {
        filename += "/alugrid"; 
      }
      return filename;
    }

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,filename)  */
    static void backup ( const Grid &grid, const std::string &filename )
    {
      std::ofstream file( filename.c_str() );
      if( file ) 
      {
        // call backup on grid 
        backup( grid, file );
      }
      else 
      {
        std::cerr << "ERROR: BackupRestoreFacility::backup: couldn't open file `" << filename << "'" << std::endl;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::backup(grid,stream)  */
    static void backup ( const Grid &grid, std::ostream &stream )
    {
      // call backup on grid 
      grid.backup( stream );
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(filename) */
    static Grid *restore ( const std::string &filename )
    {
      // Problem here: how to pass boundary projections 
      std::ifstream file( filename.c_str() );
      if( file ) 
      {
        return restore( file );
      }
      else 
      {
        std::cerr << "ERROR: BackupRestoreFacility::restore: couldn't open file `" << filename << "'" << std::endl;
        return 0;
      }
    }

    /** \copydoc Dune::BackupRestoreFacility::restore(stream) */
    static Grid *restore ( std::istream &stream )
    {
      // Problem here: how to pass boundary projections 
      Grid* grid = new Grid();
      grid->restore( stream );
      return grid;
    }
  };

} // namespace Dune

#endif // #ifndef DUNE_GRID_ALUGRID_BACKUPRESTORE_HH
