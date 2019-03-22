#include <config.h>
#include <cassert>

#include "../common/grid.hh"
#include "../common/sparsematrix.hh"
#include "../common/quadrature.hh"
#include "../common/quadrature5.hh"
#include "../common/discretefunction.hh"

#include "../problem/problem6.hh"
#include "quadraticdiscretefunction.hh"
#include "CRdiscretefunction.hh"

// use for example 'make CXXFLAGS="-O3 -DLINEAR=1" laplace'
#ifdef LINEAR
typedef DiscreteFunction<LinearBaseFunction> DiscreteFunctionType;
#elif defined QUADRATIC
typedef DiscreteFunction<QuadraticBaseFunction> DiscreteFunctionType;
#elif defined CR
typedef DiscreteFunction<CRBaseFunction> DiscreteFunctionType;
#else
#error NO FE DEFINED
#endif

using namespace std;

void solve(const GridType& grid,
           const Problem& problem,
           DiscreteFunctionType &uh)
{
  // number of local dofs
  const int locNrDof = DiscreteFunctionType::BaseFunctionSetType::locNrDof;
  const DiscreteFunctionType::BaseFunctionSetType&
        uhbase = uh.baseFunctionSet();

  // set all dofs to zero
  uh.clear();

  // right hand side
  DiscreteFunctionType F(grid);
  F.clear();

  // system matrix M
  SparseMatrix<double> M( uh.size(), 35 );

  // create quadrature
  ThirdQuadrature quad;

  // start grid traversal
  const ElementIteratorType endit = grid.end();
  for(ElementIteratorType it = grid.begin();
      it != endit; ++it)
  {
    // get reference to element
    const ElementType& element = *it;

    /////////////////////////////////////////
    //  assemble right hand side
    /////////////////////////////////////////
    {
      double I[ locNrDof ];
      for (int k=0; k<locNrDof; ++k)
      {
        I[k] = 0.0;
      }

      if( problem.Q4() )
      {
        GlobalCoordType gradPhi[ locNrDof ];
        GlobalCoordType du;
        GlobalCoordType gradBase;

          // right hand side
          for (int qp=0; qp<quad.nop(); ++qp)
          {
            // get integration weight
            const double weight = quad.weight( qp )
                                * element.integrationElement();

            // evaluate right hand side function
            const double fx = problem.f( element.global( quad.point(qp) ) );
            problem.du( element.global(quad.point(qp) ), du );

            for (int k=0; k<locNrDof; ++k)
            {
              gradBase = uh.gradient( element, k , quad.point(qp ) );
              // for all quadrature points evaluate integrand
              I[k] += fx * weight *
                      uh.evaluate(element, k, quad.point(qp) );
              I[k] += weight * ( du[0]*gradBase[0] + du[1]*gradBase[1] );
              assert( I[k] == I[k] );
            }
          }
      }
      else
      {

        // right hand side
        for (int qp=0; qp<quad.nop(); ++qp)
        {
          // get integration weight
          const double weight = quad.weight( qp )
                              * element.integrationElement();

          // evaluate right hand side function
          const double fx = problem.f( element.global( quad.point(qp) ) );

          for (int k=0; k<locNrDof; ++k)
          {
            // for all quadrature points evaluate integrand
            I[k] += fx * weight *
                    uh.evaluate(element, k, quad.point(qp) );
            assert( I[k] == I[k] );
          }
        }
      }

      // add to right hand side
      for (int k=0; k<locNrDof; ++k)
      {
        F.dof(element,k) += I[k];
      }
    }


    /////////////////////////////////////
    // setup system matrix
    ////////////////////////////////////

    {
      // local matrix storage
      double matrix[ locNrDof ][ locNrDof ];
      // global index of dofs on element
      int idx[ locNrDof ] ;
      // value of base functions
      double phi[ locNrDof ];
      //added this in for the assignment
      GlobalCoordType gradPhi[ locNrDof ];


      for (int k=0; k< locNrDof; ++k)
      {
        idx[k] = uhbase.map(element,k);
        for(int j=0; j< locNrDof; ++j)
        {
          matrix[k][j] = 0.0;
        }
      }

      for (int qp=0; qp<quad.nop(); ++qp)
      {
        // global coordinate of quadrature point
        const double weight = quad.weight( qp )
                            * element.integrationElement();

        // get gradients of basis function
        for (int k=0; k< locNrDof; ++k)
        {
          phi[k] = uh.evaluate(element, k, quad.point(qp) );
          gradPhi[k] = uh.gradient( element, k , quad.point(qp) );
        }
        // define the gradients of the phi[k]

        for (int k=0; k< locNrDof; ++k)
        {
          // eval diagonal element
          matrix[k][k] += weight * (phi[k] * phi[k]);
          // add grad of the phi


          matrix[k][k] += weight * ( gradPhi[k][0]*gradPhi[k][0] + gradPhi[k][1]*gradPhi[k][1] );

          // eval others
          for (int j=0; j< k; ++j )
          {
            double I = weight * (phi[k] * phi[j]);
            I += weight * ( gradPhi[k][0]*gradPhi[j][0] + gradPhi[k][1]*gradPhi[j][1] );
            matrix[k][j] += I;
            matrix[j][k] += I;
          }
        }
      }
      // add to global matrix
      for (int k=0; k< locNrDof; ++k)
      {
        for (int j=0; j< locNrDof; ++j )
        {
          M.add( idx[k] , idx[j] , matrix[k][j] );
        }
      }
    }
  }
  // boundary treatment
  {
    // loop over all elements
    for(ElementIteratorType it = grid.begin(); it != endit; ++it)
    {
      const ElementType& element = * it;
      for ( int edge = 0; edge < element.nCorners; ++edge )
      {
        if(problem.useDirichlet() )
        {
          for(int i = 0; i < locNrDof; i++)
              {
                if( element.onBoundary(edge) )
                    {
                      if( uhbase.onBnd(edge,i) )
                      {
                           F.dof(element,i) = problem.g(element.global( uhbase.point(i) ) );
                           uh.dof(element,i) = problem.g( element.global( uhbase.point(i) ));
                           M.setDiag( uhbase.map(element,i) , 1);
                      }
                    }
              }
        }
        // use method 'Element::onBoundary' to determine if the edge is part of the
        // boundary of Omega...
        // If it is: determine which degrees of freedom are on this edge -
        //           use the method 'onBnd' from the shape function set
        //           For those degrees of freedom set the components of the
        //           right hand side F and the solution vector uh to the boundary values
        //           (use problem.g, F.dof, and uh.dof). Finally with M.setDiag you
        //           can change the row in the sparse matrix to delta_{ij}...
      }
    } // end element iterator
  }
  //call cghssolver
  cghs( M, F, uh , 1e-12, false );
}

//////////////////////////////////////////////////////////////
//
// calculate the error of between the exact solution and a
// linear discrete function
//
//////////////////////////////////////////////////////////////
void error(const GridType& grid,
           const Problem& problem,
           const DiscreteFunctionType &uh,
           double &l2error,
           double &h1error)

{
  l2error = 0.;
  h1error = 0.;

  // create quadrature
  QuinticQuadrature quad;

  // loop over all elements
  const ElementIteratorType endit = grid.end();
  for(ElementIteratorType it = grid.begin();
      it != endit; ++it)
  {
    // get element
    const ElementType& element = *it;

    double locl2error=0;
    double loch1error=0.;

    // for all quadrature points evaluate integrand
    for (int qp=0; qp<quad.nop(); ++qp)
    {
      // global coordinate of quadrature point
      GlobalCoordType x = element.global( quad.point(qp) );

      // u( x )
      double locError = problem.u( x );

      // u_h( x )
      locError -= uh.evaluate( element, quad.point( qp ) );

      // add to local L2 error
      locl2error += quad.weight( qp ) * element.integrationElement()
                    * (locError * locError);

      // grad u ( x )
      GlobalCoordType du;
      problem.du( x, du );

      // minus grad u_h ( x )
      du -= uh.gradient( element, quad.point(qp) ) ;

      // H1-Fehler
      loch1error += quad.weight(qp) * element.integrationElement()
                    * ( du * du ) ;
    }

    // add local error to global error
    l2error += locl2error;
    h1error += loch1error;
  }
  l2error = sqrt(l2error);
  h1error = sqrt(h1error);
}

void compute(GridType& grid,
             const Problem& problem,
             const int refines)
{
  const int dimension = grid.dimension;
  double l2error[refines];
  double h1error[refines];
  double eocl2[refines];
  double eoch1[refines];

  double hold=1./sqrt(grid.size(dimension));

  // start computation
  for (int n=0; n < refines; ++n)
  {
    double hnew=1./sqrt(grid.size(dimension));

    // approximate solution
    DiscreteFunctionType uh(grid);

    // compute approximate solution
    solve(grid,problem,uh);

    // compute the approximation error
    error( grid, problem, uh, l2error[n], h1error[n] );

    if (n>0)
    {
      const double hfraction=hold/hnew;
      eocl2[n]=(log(l2error[n-1])-log(l2error[n]))/
                log(hfraction);

      eoch1[n]=(log(h1error[n-1])-log(h1error[n]))/
                log(hfraction);
    }
    else
    {
      eocl2[0]=-1;
      eoch1[0]=-1;
    }

    cout << uh.size() << "\t  "
         << hnew << " "
         << " \t\t "
         << l2error[n] << " \t " << eocl2[n]
         << " \t\t "
         << h1error[n] << " \t " << eoch1[n]
         << " \t\t "
         << endl;

    // new h is old h now
    hold = hnew;

    // optput data
    Output output( grid, n );
    output.add( uh );
    output.write();

    // refine all elements of grid twice to half grid width
    grid.globalRefine( GridType :: dimension );

  } // end of loop for (n=..)
}

/**************************************************************/
int main(int argc, char ** argv)
{
  // read in command line parameters
  string gridname ( "../problem/cube.dgf" );
  int loops = 2;
  int prob = 1;

  if (argc<4)
  {
    cout << "Usage: "<< argv[0] << " <grid file> <loops> <problem> " << endl;
    cout << "                         loops > 0" << endl;
    cout << "                         problem = 0" << endl;
  }
  else
  {
    gridname = argv[1];
    loops    = atoi( argv[2] );
    prob     = atoi( argv[3] );
  }

  Problem *problem=0;
  switch (prob)
  {
    case 0: problem = new Problem2019_Ass2_1();
            break;
    case 1: problem = new Problem2019_Ass2_2();
            break;
    case 2: problem = new Problem2019_Ass2_3();
            break;
    default: std::cerr << "Wrong problem number " << prob
                       << ". Should be less than 3." << std::endl;
             return 1;
  }

  // initialize grid
  GridType grid ( gridname.c_str() ); // read the macro grid file
  // make initial grid
  grid.globalRefine( 4 );

  compute( grid, *problem, loops );

  delete problem;

  return 0;
}
