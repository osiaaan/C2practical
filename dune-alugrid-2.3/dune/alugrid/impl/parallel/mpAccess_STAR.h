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

// (c) bernhard schupp 1997 - 1998
#ifndef MPACCESS_STARMPI_H_INCLUDED
#define MPACCESS_STARMPI_H_INCLUDED

#include "mpAccess_MPI.h"

class MpAccessSTAR_MPI : public MpAccessMPI 
{
  typedef MpAccessMPI BaseType;
protected:  
  using BaseType :: _mpiCommPtr;
  using BaseType :: _psize;
  using BaseType :: _myrank;

  int star_allgather (int *, int , int *, int) const ;
  int star_allgather (char *, int, char *, int) const ;
  int star_allgather (double *, int, double *, int ) const ;

  void initStarMPI();
public :
  // constructor taking MPI_Comm 
  // to avoid MPI types here this is a template constructor 
  template <class MPICommunicator>  
  inline MpAccessSTAR_MPI (MPICommunicator mpicomm ) 
    : BaseType( mpicomm )
  {
    initStarMPI();
  }

  // destructor 
  ~MpAccessSTAR_MPI();

  // copy constructor 
  inline MpAccessSTAR_MPI (const MpAccessSTAR_MPI &a ) 
    : BaseType( a ) 
  {}

public:  
  int gmax (int) const ;
  int gmin (int) const ;
  int gsum (int) const ;
  long gmax (long) const ;
  long gmin (long) const ;
  long gsum (long) const ;
  double gmax (double) const ;
  double gmin (double) const ;
  double gsum (double) const ;
  void gmax (double*,int,double*) const ;
  void gmin (double*,int,double*) const ;
  void gsum (double*,int,double*) const ;
  void gmax (int*,int,int*) const ;
  void gmin (int*,int,int*) const ;
  void gsum (int*,int,int*) const ;
  pair<double,double> gmax (pair<double,double>) const ;
  pair<double,double> gmin (pair<double,double>) const ;
  pair<double,double> gsum (pair<double,double>) const ;
  vector < int > gcollect (int) const ;
  vector < double > gcollect (double) const ;
  vector < vector < int > > gcollect (const vector < int > &) const ;
  vector < vector < double > > gcollect (const vector < double > &) const ;
  vector < ObjectStream > gcollect (const ObjectStream &, const vector<int>& ) const ;
  using MpAccessMPI :: gcollect ;
  using MpAccessMPI :: gmax ;
} ;
#endif
