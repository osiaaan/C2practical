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

#ifndef __MAPP_TETRA_3D_EXT_HH__
#define __MAPP_TETRA_3D_EXT_HH__

// BSGridALUVecType is defined in BSGrid interface to dune or in gitter_sti.hh 
class BSGridLinearSurfaceMapping
{
protected:
  double _n [3] ;
  const double (&_p0)[3], (&_p1)[3], (&_p2)[3] ;
public: 
  inline BSGridLinearSurfaceMapping 
    (const double (&)[3], const double (&)[3], const double (&)[3]);
    
  // same as method normal of LinearSurfaceMapping, just for Dune Vecs 
  inline void normal(double * normal) const;
  
  // same as method normal of LinearSurfaceMapping, just for Dune Vecs 
  template <class aluvec_t> 
  inline void normal(aluvec_t & normal) const;
};

class BSGridBilinearSurfaceMapping {
protected:
  double _n[3][3];
  double _b[4][3];
  const double (&_p0)[3];
  const double (&_p1)[3];
  const double (&_p2)[3];
  const double (&_p3)[3];
public:
  inline BSGridBilinearSurfaceMapping(const double(&)[3], const double(&)[3],
                                      const double(&)[3], const double(&)[3]);

  inline void normal(const double* local, double* normal) const;

  template <class aluvec_t>
  inline void normal(const aluvec_t& local, aluvec_t& normal) const;
};

inline BSGridLinearSurfaceMapping :: 
BSGridLinearSurfaceMapping (const double (&x0)[3], 
                            const double (&x1)[3], const double (&x2)[3])
  //: LinearSurfaceMapping (x0,x1,x2)
  : _p0 (x0), _p1 (x1), _p2 (x2)
{
  // copied from LinearSurfaceMapping, ist a bit faster 
  _n[0] = -0.5 * ((_p1[1]-_p0[1]) *(_p2[2]-_p1[2]) - (_p2[1]-_p1[1]) *(_p1[2]-_p0[2])) ;
  _n[1] = -0.5 * ((_p1[2]-_p0[2]) *(_p2[0]-_p1[0]) - (_p2[2]-_p1[2]) *(_p1[0]-_p0[0])) ;
  _n[2] = -0.5 * ((_p1[0]-_p0[0]) *(_p2[1]-_p1[1]) - (_p2[0]-_p1[0]) *(_p1[1]-_p0[1])) ;
}

inline void BSGridLinearSurfaceMapping :: normal (double * normal) const
{
  normal[0] = this->_n[0];
  normal[1] = this->_n[1];
  normal[2] = this->_n[2];
  return ;
}

template <class aluvec_t>
inline void BSGridLinearSurfaceMapping :: normal (aluvec_t & normal) const
{
  normal[0] = this->_n[0];
  normal[1] = this->_n[1];
  normal[2] = this->_n[2];
  return ;
}

inline BSGridBilinearSurfaceMapping::
BSGridBilinearSurfaceMapping(const double(&x0)[3], const double(&x1)[3],
                             const double(&x2)[3], const double(&x3)[3]) :
  _p0 (x0), _p1 (x1), _p2 (x2), _p3 (x3) {
  _b [0][0] = _p0 [0] ;
  _b [0][1] = _p0 [1] ;
  _b [0][2] = _p0 [2] ;
  _b [1][0] = _p3 [0] - _p0 [0] ;
  _b [1][1] = _p3 [1] - _p0 [1] ;
  _b [1][2] = _p3 [2] - _p0 [2] ;
  _b [2][0] = _p1 [0] - _p0 [0] ;
  _b [2][1] = _p1 [1] - _p0 [1] ;
  _b [2][2] = _p1 [2] - _p0 [2] ;
  _b [3][0] = _p2 [0] - _p1 [0] - _b [1][0] ;
  _b [3][1] = _p2 [1] - _p1 [1] - _b [1][1] ;
  _b [3][2] = _p2 [2] - _p1 [2] - _b [1][2] ;
  _n [0][0] = _b [1][1] * _b [2][2] - _b [1][2] * _b [2][1] ;
  _n [0][1] = _b [1][2] * _b [2][0] - _b [1][0] * _b [2][2] ;
  _n [0][2] = _b [1][0] * _b [2][1] - _b [1][1] * _b [2][0] ;
  _n [1][0] = _b [1][1] * _b [3][2] - _b [1][2] * _b [3][1] ;
  _n [1][1] = _b [1][2] * _b [3][0] - _b [1][0] * _b [3][2] ;
  _n [1][2] = _b [1][0] * _b [3][1] - _b [1][1] * _b [3][0] ;
  _n [2][0] = _b [3][1] * _b [2][2] - _b [3][2] * _b [2][1] ;
  _n [2][1] = _b [3][2] * _b [2][0] - _b [3][0] * _b [2][2] ; 
  _n [2][2] = _b [3][0] * _b [2][1] - _b [3][1] * _b [2][0] ;
  return ;
}

inline void 
BSGridBilinearSurfaceMapping::normal(const double* local, double* normal) const {
  double x = local [0];
  double y = local [1];
  normal [0] = -( _n [0][0] + _n [1][0] * x + _n [2][0] * y) ;
  normal [1] = -( _n [0][1] + _n [1][1] * x + _n [2][1] * y) ;
  normal [2] = -( _n [0][2] + _n [1][2] * x + _n [2][2] * y) ;
  return ;
}

template <class aluvec_t>
inline void BSGridBilinearSurfaceMapping::normal(const aluvec_t& local, aluvec_t& normal) const {
  double x = local [0];
  double y = local [1];
  normal [0] = -( _n [0][0] + _n [1][0] * x + _n [2][0] * y) ;
  normal [1] = -( _n [0][1] + _n [1][1] * x + _n [2][1] * y) ;
  normal [2] = -( _n [0][2] + _n [1][2] * x + _n [2][2] * y) ;
  return ;
}



#endif
