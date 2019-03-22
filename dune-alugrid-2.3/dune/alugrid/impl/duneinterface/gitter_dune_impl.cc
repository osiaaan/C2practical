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

// Robert Kloefkorn (c) 2004 - 2005 
#include <config.h>

#include <fstream>

#include "gitter_dune_impl.h"

namespace ALUGrid
{

  IteratorSTI < Gitter::helement_STI > * GitterDuneImpl::
  leafIterator (const helement_STI *) 
  {
    return new Insert < PureElementAccessIterator < Gitter::helement_STI >::Handle,
    TreeIterator < Gitter::helement_STI, is_leaf < Gitter::helement_STI> > > (container ());
  }

  IteratorSTI < Gitter::helement_STI > * GitterDuneImpl::
  leafIterator (const IteratorSTI < helement_STI > * p) 
  {
    return new Insert < PureElementAccessIterator < Gitter::helement_STI >::Handle,
    TreeIterator < Gitter::helement_STI, is_leaf < Gitter::helement_STI> > > 
    (*(const Insert < PureElementAccessIterator < Gitter::helement_STI >::Handle,
    TreeIterator < Gitter::helement_STI, is_leaf < Gitter::helement_STI> > > *) p);
  }

  // wird von Dune verwendet 
  void GitterDuneBasis::duneBackup (const char * fileName) 
  {
    // diese Methode wird von der Dune Schnittstelle aufgerufen und ruft
    // intern lediglich backup (siehe oben) und backupCMode des Makrogitters
    // auf, allerdings wird hier der path und filename in einer variablen
    // uebergeben 

    alugrid_assert (debugOption (20) ? (std::cout << "**INFO GitterDuneBasis::duneBackup (const char * = \""
                         << fileName << "\") " << std::endl, 1) : 1);

    std::ofstream out( fileName );
    if( !out )
    {
      std::cerr << "WARNING (ignored): Error creating '" << (fileName ? fileName : "") << "' in GitterDuneBasis::duneBackup( const char *, double )." << std::endl;
    }
    else
    {
      Gitter::backup (out);

      backupIndices ( out );

      {
        std::string fullName = std::string( fileName ) + std::string( ".macro" );
        std::ofstream macro( fullName.c_str() );
        if( !macro )
          std::cerr << "WARNING (ignored): Cannot create '" << fullName << "' in GitterDuneBasis::duneBackup( const char *, const char * )." << std::endl;
        else
          container ().backupCMode( macro );
      }
    }
  }

  // wird von Dune verwendet 
  void GitterDuneBasis::duneRestore (const char * fileName)
  {
    // diese Methode wird von der Dune Schnittstelle aufgerufen 
    // diese Methode ruft intern restore auf, hier wird lediglich 
    // der path und filename in einer variablen uebergeben

    alugrid_assert (debugOption (20) ? (std::cout << "**INFO GitterDuneBasis::duneRestore (const char * = \""
                   << fileName << "\") " << std::endl, 1) : 1);

    std::ifstream in( fileName );
    if( !in )
      std::cerr << "WARNING (ignored): Cannot open '" << (fileName ? fileName : "") << "' in GitterDuneBasis::duneRestore( const char *, double & )." << std::endl;
    else
    {
      Gitter::restore( in );
      this->restoreIndices( in );
    }
  }

  int GitterDuneBasis::preCoarsening  (Gitter::helement_STI & elem)
  {
    // if _arp is set then the extrenal preCoarsening is called 
    return (_arp) ? (*_arp).preCoarsening(elem) : 0;
  }

  int GitterDuneBasis::postRefinement (Gitter::helement_STI & elem)
  {
    // if _arp is set then the extrenal postRefinement is called 
    return (_arp) ? (*_arp).postRefinement(elem) : 0;
  }

  int GitterDuneBasis::preCoarsening  (Gitter::hbndseg_STI & bnd)
  {
    // if _arp is set then the extrenal preCoarsening is called 
    return (_arp) ? (*_arp).preCoarsening( bnd ) : 0;
  }

  int GitterDuneBasis::postRefinement (Gitter::hbndseg_STI & bnd)
  {
    // if _arp is set then the extrenal postRefinement is called 
    return (_arp) ? (*_arp).postRefinement( bnd ) : 0;
  }

  void GitterDuneBasis ::
  setAdaptRestrictProlongOp( AdaptRestrictProlongType & arp )
  {
    if( _arp ) 
      std::cerr << "WARNING (ignored): _arp not null in GitterDuneBasis::setAdaptRestrictProlongOp." << std::endl;
    _arp = &arp;
  }

  void GitterDuneBasis::removeAdaptRestrictProlongOp()
  {
    _arp = 0;
  }

  bool GitterDuneBasis::duneAdapt (AdaptRestrictProlongType & arp) 
  {
    alugrid_assert (debugOption (20) ? (std::cout << "**INFO GitterDuneBasis::duneAdapt ()" << std::endl, 1) : 1);
    // set restriction-prolongation callback obj
    setAdaptRestrictProlongOp(arp); 
    // call adapt method 
    const bool refined = this->adapt();
    // sets pointer to zero 
    removeAdaptRestrictProlongOp ();
    return refined;
  }

} // namespace ALUGrid
