//============================================================================
//
//  CG method after Hestenes und Stiefel
//
//  see also:
//  Ashby, Manteuffel, Saylor
//     A taxononmy for conjugate gradient methods
//     SIAM J Numer Anal 27, 1542-1568 (1990)
//
//  ----------------------------
//  revised to work without system BLAS routines 
//      Tobias Malkmus, Juli 2008
//
//============================================================================


#include <iostream>
#include <cassert>
#include <cmath>

using namespace std;

#include "vector.hh"

/** 
   @addtogroup Solver 
   @{
*/


/** \brief CG scheme after Hestenes und Stiefel 
    \param N number of unknowns 
    \param A matrix (implement mult method for this matrix)
    \param b right hand side 
    \param x unknown 
    \param eps error bound 
    \param detailed flag for toggeling output 
 */
template< class MATRIX >
inline int
cghs( const unsigned int N, 
      const MATRIX &A, 
      const double *b, 
      double *x, 
      const double eps,
      const bool detailed = false );

/** \brief CG scheme after Hestenes und Stiefel 
    \param A matrix (implement mult method for this matrix)
    \param b right hand side 
    \param x unknown 
    \param eps error bound 
    \param detailed flag for toggeling output (default is \b false)
 */
template< class T, class MATRIX >
inline int
cghs( const MATRIX &A, 
      const Vector<T> &b, 
      Vector<T>& x, 
      const double eps,
      const bool detailed = false )
{
  assert( b.size() == x.size() );
  return cghs(b.size(), A, b.raw(), x.raw(), eps, detailed );
}

namespace { // anonymous name space only so that doxygen does not add to documentation
//============================================================================
// some BLAS routines 
//============================================================================

inline void addscalar(const int N, const double alpha, 
                      const double *x,double *y)
{
  for(int k=0; k<N; ++k)
  {
    y[k] += alpha * x[k];
  }
}
inline double scalarproduct(const int N,const double *x,const double *y)
{
  double ret =0.;
  for(int k=0; k<N; ++k)
  {
    ret += x[k] * y[k];
  }
  return ret;
}
inline void rescale(const int N, const double alpha, double *x)
{
  for(int k=0; k<N; ++k)
  {
    x[k] *= alpha;
  }
}

inline void copy(const int N, const double *x, double *y)
{
  for(int k=0; k<N; ++k)
  {
    y[k] = x[k];
  }
}
}

template <class MATRIX>
inline int cghs( const unsigned int N, 
                 const MATRIX &A, 
                 const double *b, 
                 double *x, 
                 const double eps,
                 const bool detailed ) 
{
  // for N zero do nothing 
  if ( N==0 ) return -1;
  
  double *g = new double[N];
  double *r = new double[N];
  double *p = new double[N];
  double t, tau, sig, rho, gam;
  double err=eps*eps; // *scalarproduct(N,b,b);    //scalar produkt
  
  int its=0;

  // initial multiplication 
  A.mult(x,g);
  addscalar(N,-1.,b,g);   //g=-1b+g
  rescale(N,-1.,g);       //g = -g
  copy(N,g,r);            //r=g
  
  // loop until scalar product of residual is small 
  while ( scalarproduct(N,g,g) > err ) 
  {
    // matrix multiplication 
    A.mult(r,p);
    
    // evaluate scalar products 
    rho=scalarproduct(N,p,p);
    sig=scalarproduct(N,r,p);
    tau=scalarproduct(N,g,r);
    assert(sig == sig);
    assert(tau == tau);
    t=tau/sig;
    
    addscalar(N,t,r,x);    //x= t*r+x
    addscalar(N,-t,p,g);   //g= -t*p +g
    
    gam=(t*t*rho-tau)/tau;
    
    // apply correction 
    rescale(N,gam,r);      //r=r*gam
    addscalar(N,1.,g,r);   //r=g+r
    
    if ( detailed )
    {
      cout<<"cghs "<<its<<"\t"<<sqrt(scalarproduct(N,g,g))<<endl;
    }
    
    ++its;
  }
  
  if ( detailed )
  {
    cout<<"cghs "<<its<<"\t"<<sqrt(scalarproduct(N,g,g))<<endl;
  }
  
  delete[] g;
  delete[] r;
  delete[] p;
  
  return its;
}

/** 
   @}
*/

