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

#ifndef Problem_BALL_HH
#define Problem_BALL_HH

#include <cmath>
#include <cassert>
#include <iostream>
#include <sstream>

#include <dune/common/fvector.hh>

#include "problem-transport.hh"
#include "problem.hh"

template< int dimD >
struct BallData
: public ProblemData< dimD,1 >
{
  typedef ProblemData< dimD,1 > Base;

  typedef typename Base::DomainType DomainType;
  typedef typename Base::RangeType RangeType;

  const static int dimDomain = DomainType::dimension;

  explicit BallData (const int problem) : c_(0.5), r0_(0.3), problem_( problem )
  {}

  //! \copydoc ProblemData::gridFile
  std::string gridFile ( const std::string &path ) const
  { 
    std::ostringstream dgfFileName;
    if( problem_ == 1 ) 
      dgfFileName << path << "/dgf/cube_hc_512.dgf";
    else if ( problem_ == 2 ) 
      dgfFileName << path << "/dgf/cube_hc_4096.dgf";
    else if ( problem_ == 3 ) 
      dgfFileName << path << "/dgf/cube_hc_32768.dgf";
    else if ( problem_ == 4 ) 
      dgfFileName << path << "/dgf/cube_hc_262144.dgf";
    else 
      dgfFileName << path << "/dgf/unitcube3d.dgf";
    return dgfFileName.str();
  }

  //! \copydoc ProblemData::endTime
  double endTime () const
  {
    return 1.0;
  }

  int bndType( const DomainType &normal, const DomainType &x, const double time) const 
  {
    return 0;
  }

  //! \copydoc ProblemData::adaptationIndicator
  double adaptationIndicator ( const DomainType& x, double time,
                               const RangeType& uLeft, const RangeType &uRight ) const 
  { 
    DomainType xx(x);
    xx -= c_;
    DomainType y(0);
    y[0] = std::cos(time*2.*M_PI)*r0_;
    y[1] = std::sin(time*2.*M_PI)*r0_;
    xx -= y;
    double r = xx.two_norm();
    return ( (r>0.15 && r<0.25)? 1 : 0 );
  } 

  //! \copydoc ProblemData::refineTol
  double refineTol () const
  {
    return 0.1;
  }

  //! \copydoc ProblemData::saveInterval
  double saveInterval() const 
  {
    return 0.1;
  }

private:
  DomainType c_;
  double r0_;
  const int problem_;
};

// BallModel
// ------------

/** \brief Problem describing the Euler equation of gas dynamics 
 */
template <int dimD>
struct BallModel : public TransportModel<dimD>
{
  typedef ProblemData< dimD,1 > Problem;

  typedef typename Problem::DomainType DomainType;
  typedef typename Problem::RangeType RangeType;

  static const int dimDomain = DomainType::dimension;
  static const int dimRange = RangeType::dimension;
  static const bool hasFlux = false;

  /** \brief constructor 
   *  \param problem switch between different data settings 
   */
  BallModel( unsigned int problem )
  : problem_( 0 )
  {
    switch( problem )
    {
    case 0:
    case 1:
    case 2:
    case 3:
      problem_ = new BallData< dimDomain >( problem );
      break;

    default:
      std::cerr << "ProblemData not defined - using problem 1!" << std::endl;
      problem_ = new BallData< dimDomain >( problem );
    }
  }

  /** \brief destructor */
  ~BallModel() 
  {
    delete problem_;
  }

  /** \copydoc TransportProblem::data */
  const Problem &problem () const
  {
    return *problem_;
  }

  double fixedDt () const
  {
    return 0.01;
  }

  /** \copydoc TransportProblem::indicator */
  double indicator ( const DomainType &normal,
                     const double time,
                     const DomainType &xGlobal,
                     const RangeType &uLeft, const RangeType &uRight) const 
  {
    return problem().adaptationIndicator( xGlobal, time, uLeft, uRight );
  }

  /** \copydoc TransportProblem::boundaryIndicator */
  double boundaryIndicator ( const int bndId, 
                             const DomainType &normal, 
                             const double time,
                             const DomainType &xGlobal,
                             const RangeType& uLeft) const
  {
    return indicator( normal,time,xGlobal, uLeft, problem().boundaryValue(xGlobal,time) );
  }

  Problem *problem_;
};

#endif 
