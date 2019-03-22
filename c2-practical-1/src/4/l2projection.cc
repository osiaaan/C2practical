#include <config.h>
#include <cassert>

#include "../common/grid.hh"
#include "../common/sparsematrix.hh"
#include "../common/quadrature.hh"
#include "../common/discretefunction.hh"

#include "../problem/problem4.hh"

typedef DiscreteFunction<LinearBaseFunction> DiscreteFunctionType;

using namespace std;

//////////////////////////////////////////////////////////////
//
// calculate the linear interpolation of a given function
// (this is taken from the previous assignment and just here to
// allow an easy comparison of interpolation and L2-projection
//
//////////////////////////////////////////////////////////////
void interpolate(const GridType& grid,
                 const Problem& problem,
                 DiscreteFunctionType &uh)
{
  // number of local dofs
  const int locNrDof = DiscreteFunctionType::BaseFunctionSetType::locNrDof;
  const DiscreteFunctionType::BaseFunctionSetType&
        uhbase = uh.baseFunctionSet();

  // set all dofs to zero
  uh.clear();

  // start grid traversal
  const ElementIteratorType endit = grid.end();
  for(ElementIteratorType it = grid.begin();
      it != endit; ++it)
  {
    // get reference to element
    const ElementType& element = *it;

    // for all vertices do interpolation
    for(int i=0; i < locNrDof; ++i)
    {
      uh.dof(element, i) = problem.u( element.global( uhbase.point(i) ) );
    }
  }
}
//////////////////////////////////////////////////////////////
//
// calculate the l2 projection of a given function into the
// linear lagrange space
//
//////////////////////////////////////////////////////////////
void l2projection(const GridType& grid,
                  const Problem& problem,
                  DiscreteFunctionType &uh)
{
  // number of local dofs (hatN)
  const int locNrDof = DiscreteFunctionType::BaseFunctionSetType::locNrDof;

  // set all dofs to zero
  uh.clear();

  // right hand side (\int_\Omega f\varphi)
  DiscreteFunctionType F(grid);
  F.clear();

  // system matrix M (\int_\Omega \varphi_j\varphi_i)
  SparseMatrix<double> M( uh.size(), 35 );

  // create quadrature (of second order, one for phi and one for f)
  LineQuadrature quad;

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
      // storage for local right hand side (\int_T u\phi^T_k)
      double I[ locNrDof ];
      for (int k=0; k<locNrDof; ++k)
        I[k] = 0.0;

      // right hand side assembly using a quadrature rule
      const int nop = quad.nop();
      for (int qp=0; qp<nop; ++qp)
      {
        // get integration weight
        const double weight = quad.weight( qp )
                            * element.integrationElement();

        // evaluate function to interpolate
        const double ux = problem.u( element.global( quad.point(qp) ) );

        for (int k=0; k<locNrDof; ++k)
        {
          // for all quadrature points evaluate integrand (u(x)phi^k_T(x))
          I[k] += ux * weight *
                  uh.evaluate(element, k, quad.point(qp) );
          assert( I[k] == I[k] );
        }
      }

      // add to right hand side
      for (int k=0; k<locNrDof; ++k)
      {
        // add to F_{\vu(K,k)}
        F.dof(element,k) += I[k];
      }
    }

    /////////////////////////////////////
    // setup system matrix
    ////////////////////////////////////

    {
      // local matrix storage
      double matrix[ locNrDof ][ locNrDof ];

      double sum = 0.0;
      double const nop = quad.nop();

      for(int i = 0 ; i < locNrDof ; ++i)
      {
          for(int j = 0; j < locNrDof : j++)
          {
            for(int qp = 0; qp < nop; ++qp)
            {
                const double weight = quad.weight(qp) * element.integrationElement();
                sum += weight * uh.evaluate(element, i, quad.point(gp)) * uh.evaluate(element, j, quad.point(gp));
            }
            matrix[i][j] = sum;
          }
      }

      //global indexing
      int idx[ locNrDof ];

      for (int k=0; k< locNrDof; ++k)
      {
        idx[k] = uhbase.map(element,k);
      }

      for (int k=0; k< locNrDof; ++k)
      {
        for (int j=0; j< locNrDof; ++j )
        {
          M.add( idx[k] , idx[j] , matrix[k][j] );
        }
      }


      const DiscreteFunctionType::BaseFunctionSetType&
    uhbase = uh.baseFunctionSet();


        add(const int row , const int column, const T &value)

      //////////////////////////////////////////////////////////////////////////
      //
      //  Exercise 4:
      //
      //  Fill in the local mass matrix (matrix) and distribute to global
      //  mass matrix M: m_ij = int_Omega phi_i phi_j
      //
      //  Step 1: set of local mass matrix m^T_kl = int_T phi^T_k phi^T_l (compare with
      //          assembly of right hand side b_i = int_Omega u phi_i
      //  Step 2: add entries m^T_kl into correct position of global mass matrix M
      //
      //  Look at documentation of SparseMatrix class:
      //  Methods needed here: add(...)
      //  Look at documentation of DiscreteFunction and LinearBaseFunction
      //  Further new class to look are the Quadrature classes
      //
      //////////////////////////////////////////////////////////////////////////
      // EXERCISE("Fill in contribution to mass matrix from element");
    }
  }

  // call cghssolver (parameters are
  // the matrix, rhs, solution vector (containing the initial guess),
  // the tolerance for the residual, and a bool flag for switching between
  // verbose and non verbose output
  // (change to true to see how the iterations)
  cghs( M, F, uh , 1e-8, false );
}

//////////////////////////////////////////////////////////////
//
// calculate the error between the exact solution and a
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
  ThirdQuadrature quad;

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
  double l2error[refines];
  double h1error[refines];
  double eocl2[refines];
  double eoch1[refines];

  double hold=1./sqrt(grid.size(2));

  // start computation
  for (int n=0; n < refines; ++n)
  {
    double hnew=1./sqrt(grid.size(2));

    // approximate solution
    DiscreteFunctionType uh(grid);

    // compute approximate solution
    interpolate(grid,problem,uh);
    // l2projection(grid,problem,uh);

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
    Output out( grid, n );
    out.addVertexData( uh, "I_h(u)" );
    out.write();

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
    cout << "                         problem = 0..2" << endl;
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
    case 0: problem = new ProblemSin();
            break;
    case 1: problem = new Problem2019_Ass1_1();
            break;
    case 2: problem = new Problem2019_Ass1_2_3(1./M_PI);
            break;
    case 3: problem = new Problem2019_Ass1_2_3(0.5);
            break;
    case 4: problem = new Problem2019_Ass1_4();
            break;
    case 5: problem = new Problem2019_Ass1_5(sqrt(2.)/2.);
            break;
    case 6: problem = new Problem2019_Ass1_5(0.5);
            break;
    default: std::cerr << "Wrong problem number " << prob
                       << ". Should be less than 7." << std::endl;
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
