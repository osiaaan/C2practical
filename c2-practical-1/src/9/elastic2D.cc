#include <config.h>
#include <cassert>

#include "../common/grid.hh"
#include "cgsolver.hh"
#include "../common/sparsematrix.hh"
#include "../common/quadrature.hh"
#include "../common/quadrature5.hh"
#include "../common/discretefunction.hh"

#include "problem_elastic2D.hh"  // TODO: change to problem_elastic2D.hh
#include "quadraticdiscretefunction.hh"
#include "crouzeixraviartdiscretefunction.hh"
#include "discretevectorfield.hh"  // TODO: needs implementing

// use for example 'make CXXFLAGS="-O3 -DLINEAR=1" laplace'
#ifdef LINEAR
typedef LinearBaseFunction BaseFunctionType;
#elif defined QUADRATIC
typedef QuadraticBaseFunction BaseFunctionType;
#elif defined CR
typedef CRBaseFunction BaseFunctionType;
#else
#error NO FE DEFINED
#endif
typedef DiscreteVectorfield<BaseFunctionType> DiscreteFunctionType;
// TODO use vector field: typedef DiscreteVectorField<BaseFunctionType> DiscreteFunctionType;



using namespace std;

void solve(const GridType& grid,
           const Problem& problem,
           DiscreteFunctionType &uh)
{ ////////////////////////////////////////////////////////////
  // TODO: switch to vector valued functions and then
  //       to elasticity
  //       Tip: could first test a vector valued
  //       laplace, i.e.,
  //         -laplace u_1 = f_1
  //         -laplace u_2 = f_2
  //       where even f_1=f_2=f could be used at first
  //       If this works then switch to ellasticity
  ////////////////////////////////////////////////////////////
  // number of local dofs
  const int locNrDof = DiscreteFunctionType::BaseFunctionSetType::locNrDof;
  const DiscreteFunctionType::BaseFunctionSetType& uhbase = uh.baseFunctionSet();
  const int dow = 2;

  // set all dofs to zero
  uh.clear();

  // right hand side
  DiscreteFunctionType F(grid);
  F.clear();

  // system matrix M
  SparseMatrix<double> M( uh.size(), 100 );

  // create quadrature
  ThirdQuadrature quad;

  // local rhs storage
  double I[locNrDof][dow];
  // local matrix storage
  double matrix[locNrDof][locNrDof][dow][dow];
  // global index of dofs on element
  // int idx[locNrDof];
  int vidx[locNrDof][dow];

  // start grid traversal
  const ElementIteratorType endit = grid.end();
  for(ElementIteratorType it = grid.begin();
      it != endit; ++it)
  {
    // get reference to element
    const ElementType& element = *it;

    for (int k=0; k<locNrDof; ++k)
    {
      I[k][0] = 0.0;
      I[k][1] = 0.0;
      vidx[k][0] = uhbase.map(element,k);
      vidx[k][1] = uhbase.map(element,k) + (uh.size())/2.0;
      for (int j=0; j<locNrDof; ++j)
      {
        matrix[k][j][0][0] = 0.0;
        matrix[k][j][0][1] = 0.0;
        matrix[k][j][1][0] = 0.0;
        matrix[k][j][1][1] = 0.0;
      }
    }

    // value of base functions
    GlobalCoordType gradphi[locNrDof];

    /////////////////////////////////////////
    //  assemble volume terms
    /////////////////////////////////////////
    {
      for (int qp=0; qp<quad.nop(); ++qp)
      {
        // global coordinate of quadrature point
        const double weight = quad.weight( qp ) * element.integrationElement();

        // evaluate right hand side function
        const GlobalCoordType fqp = problem.f( element.global( quad.point(qp) ) );

        for (int k=0; k<locNrDof; ++k)
        {
	        double phik = uh.evaluate(element, k, quad.point(qp) );
          I[k][0] += fqp[0] * weight * phik;
          I[k][1] += fqp[1] * weight * phik;
          assert( I[k] == I[k] );
	      }

        // get gradients of basis function_
        for (int k=0; k<locNrDof; ++k)
        {
          gradphi[k] = uh.gradient(element, k, quad.point(qp));
        }

        for (int k=0; k<locNrDof; ++k)
        {
         for (int j=0; j<locNrDof; ++j)
         {
           const double v_00 = (2.0*problem.get_mu()*gradphi[k][0] + problem.get_lambda()*gradphi[k][0])*gradphi[j][0] +
                                problem.get_mu()*gradphi[k][1]*gradphi[j][1];
           const double v_01 = problem.get_mu()*gradphi[k][1]*gradphi[j][0] + problem.get_lambda()*gradphi[k][0]*gradphi[j][1];
           const double v_10 = problem.get_mu()*gradphi[k][0]*gradphi[j][1] + problem.get_lambda()*gradphi[k][1]*gradphi[j][0];
           const double v_11 = (2.0*problem.get_mu()*gradphi[k][1] + problem.get_lambda()*gradphi[k][1])*gradphi[j][1] +
                                problem.get_mu()*gradphi[k][0]*gradphi[j][0];
           matrix[k][j][0][0] += v_00*weight;
           matrix[k][j][0][1] += v_01*weight;
           matrix[k][j][1][0] += v_10*weight;
           matrix[k][j][1][1] += v_11*weight;
           assert( matrix[k][j] == matrix[k][j] );
         }
       }
      } // end of element quadrature loop
    } // end of volume assembly

    /////////////////////////////////////////
    //  assemble (skeleton) and boundary terms
    /////////////////////////////////////////

    {
      // iterate over all 'intersection'
      const IntersectionIteratorType enditi= element.iend();
      for(IntersectionIteratorType iti = element.ibegin(); iti != enditi; ++iti)
      {
        const IntersectionType& intersection = *iti;
        GlobalCoordType normal; // outer unit normal scaled with integration element
        normal = intersection.outerNormal();
        double intersection_integrationelement = normal.two_norm();
        // 2 point Gauss quadrature
        const double sqrt3inv = 1./sqrt(3.);
        double x[2] = {0.5*(1.+sqrt3inv),0.5*(1.-sqrt3inv)};
        double w[2] = {0.5,0.5};
        if(intersection.boundary()) // boundary term
        {
          // test in midpoint of edge if is dirichlet boundary
          const GlobalCoordType xp = element.global( intersection.geometryInSelf( {0.5} ) );
          if (! problem.isDirichlet(xp))
          { // add non trivial Neuman contribution
            for (int qp=0; qp<2; ++qp)
            {
              const GlobalCoordType xqpIn = intersection.geometryInSelf( x[qp] );
              double weight = intersection_integrationelement*w[qp];
              GlobalCoordType hqp = problem.h( element.global( xqpIn ) );
              for (int k=0; k<locNrDof; ++k)
              {
                const double phikIn = uh.evaluate( element, k, xqpIn );
                I[k][0] += hqp[0] * phikIn * weight;
                I[k][1] += hqp[1] * phikIn * weight;
              }
            }
          }
        } // end was boundary intersection
        else
        { // skeleton terms
          const ElementType &neighbor = intersection.outside();
          double h2 = std::max( neighbor.area(), element.area() );
          double h  = sqrt(h2);
          for (int qp=0; qp<2; ++qp)
          {
            const GlobalCoordType xqpIn  = intersection.geometryInSelf( x[qp] );
            const GlobalCoordType xqpOut = intersection.geometryInNeighbor( x[qp] );
            double weight = w[qp] * intersection_integrationelement;
            // TODO: add terms for CR penalty


            double IN[locNrDof];
            double OUT[locNrDof];
            for (int k=0; k<locNrDof; ++k)
            {
              IN[k] = uh.evaluate( element, k, xqpIn );
              OUT[k] = uh.evaluate( element, k, xqpOut );
            }

            for (int k=0; k<locNrDof; ++k)
            {
             for (int l=0; l<locNrDof; ++l)
             {
               double phi_k = OUT[k] - IN[k];
               double phi_l = OUT[l] - IN[l];
               double constant = 2*problem.get_mu()*problem.get_beta();
               constant = constant/h;
               const double v_00 = constant * phi_k * phi_l * normal[0] * normal[0];
               const double v_01 = constant * phi_k * phi_l * normal[0] * normal[1];
               const double v_10 = constant * phi_k * phi_l * normal[1] * normal[0];
               const double v_11 = constant * phi_k * phi_l * normal[1] * normal[1];

               matrix[k][l][0][0] += v_00*weight;
               matrix[k][l][0][1] += v_01*weight;
               matrix[k][l][1][0] += v_10*weight;
               matrix[k][l][1][1] += v_11*weight;

               assert( matrix[k][l] == matrix[k][l] );
             }
           }

          }
        }
      }
    } // end of skeleton assembly


    // add to right hand side
    for (int k=0; k<locNrDof; ++k)
    {
      //vidx[k][0] = uhbase.map(element,k);
      //vidx[k][1] = uhbase.map(element,k) + (uh.size()/2.0);
      F.dof(element,k,0) += I[k][0];
      F.dof(element,k,1) += I[k][1];

      //F[vidx[k][0]] += I[k][0];
      //F[vidx[k][1]] += I[k][1];
    }

    // add to global matrix
    for (int k=0; k<locNrDof; ++k)
     for (int l=0; l<locNrDof; ++l)
     {
      M.add( vidx[k][0] , vidx[l][0] , matrix[k][l][0][0] );
      M.add( vidx[k][0] , vidx[l][1] , matrix[k][l][0][1] );
      M.add( vidx[k][1] , vidx[l][0] , matrix[k][l][1][0] );
      M.add( vidx[k][1] , vidx[l][1] , matrix[k][l][1][1] );
    }
  } // end element iterator

  // diagnostic:
  // M.print();
  // assert(M.symmetric());

  // Dirichlet boundary treatment
  if (problem.useDirichlet())
  {
    // loop over all elements
    for(ElementIteratorType it = grid.begin(); it != endit; ++it)
    {
      const ElementType& element = * it;
      // iterate over intersections
      const IntersectionIteratorType enditi= element.iend();
      for(IntersectionIteratorType iti = element.ibegin(); iti != enditi; ++iti)
      {
        const IntersectionType& intersection = *iti;
        if(intersection.boundary())
        {
          // test in midpoint of edge if is dirichlet boundary
          const GlobalCoordType xp = element.global( intersection.geometryInSelf( {0.5} ) );
          if (problem.isDirichlet(xp))
          {
            const int iidx = intersection.numberInSelf();
            for (int k=0; k< locNrDof; ++k)
            {
              // if dof is on boundary
              if(uhbase.onBnd(iidx,k))
              {
                const GlobalCoordType xp = element.global( uhbase.point(k) );
                // TODO: switch to vector valued 'g'
                const GlobalCoordType gxp = problem.g(xp);
                assert(gxp == gxp);
                // set right hand side and solution and adjust matrix
                F.dof(element,k,0)  = gxp[0];
                F.dof(element,k,1)  = gxp[1];
                uh.dof(element,k,0) = gxp[0];
                uh.dof(element,k,1) = gxp[1];
                M.setDiag( uhbase.map(element,k) );
                M.setDiag( uhbase.map(element,k) + (uh.size())/2.0 );
              }
            }
          }
        }
      }
    } // end element iterator
  }

  //call cghssolver
  int iterations = cgSolver( M, F, uh , 1e-8, 10000, false );
  std::cout << "cg iterations: " << iterations << std::endl;
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
           double &h1error,
           double &hEerror)

{
  l2error = 0.;
  h1error = 0.;
  hEerror = 0.;

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
    double locEerror=0.;

    // for all quadrature points evaluate integrand
    for (int qp=0; qp<quad.nop(); ++qp)
    {
      // global coordinate of quadrature point
      GlobalCoordType x = element.global( quad.point(qp) );

      // u( x )
      GlobalCoordType locError = problem.u( x );

      // u_h( x )
      locError[0] -= uh.evaluate( element, quad.point( qp ) )[0];
      locError[1] -= uh.evaluate( element, quad.point( qp ) )[1];

      // add to local L2 error
      locl2error += quad.weight( qp ) * element.integrationElement()
                    * (locError[0]*locError[0] + locError[1]*locError[1]);

      // grad u ( x )
      JacobianType du;
      du = problem.du( x );

      // TODO: switch to elasticity energy norm
      // minus grad u_h ( x )
      du -= uh.gradient( element, quad.point(qp) ) ;

      // H1-Fehler
      loch1error += quad.weight(qp) * element.integrationElement()
                    * ( du[0][0]*du[0][0] + du[0][1]*du[0][1]
                      + du[1][0]*du[1][0] + du[1][1]*du[1][1] ) ;

      locEerror += 2*problem.get_mu()*loch1error +
                   problem.get_lambda()*quad.weight(qp)*element.integrationElement() *
                   ( du[0][0] + du[1][1] )* ( du[0][0] + du[1][1] );
    }

    // add local error to global error
    l2error += locl2error;
    h1error += loch1error;
    hEerror += locEerror;
  }
  l2error = sqrt(l2error);
  h1error = sqrt(h1error);
  hEerror = sqrt(hEerror);
}

void compute(GridType& grid,
             const Problem& problem,
             const int refines)
{
  const int dimension = grid.dimension;
  double l2error[refines];
  double h1error[refines];
  double hEerror[refines];
  double eocl2[refines];
  double eoch1[refines];
  double eochE[refines];

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
    error( grid, problem, uh, l2error[n], h1error[n], hEerror[n] );

    if (n>0)
    {
      const double hfraction=hold/hnew;
      eocl2[n]=(log(l2error[n-1])-log(l2error[n]))/
                log(hfraction);

      eoch1[n]=(log(h1error[n-1])-log(h1error[n]))/
                log(hfraction);

      eochE[n]=(log(hEerror[n-1])-log(hEerror[n]))/
                log(hfraction);
    }
    else
    {
      eocl2[0]=-1;
      eoch1[0]=-1;
      eochE[0]=-1;
    }

    cout << uh.size() << "\t  "
         << hnew << " "
         << " \t\t "
         << l2error[n] << " \t " << eocl2[n]
         << " \t\t "
         << h1error[n] << " \t " << eoch1[n]
         << " \t\t "
         << hEerror[n] << " \t " << eochE[n]
         << " \t\t "
         << endl;

    // new h is old h now
    hold = hnew;

    // optput data
    OutputVF output( grid, n );
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
    cout << "Usage: "<< argv[0] << " <grid file> <loops> <problem>" << endl;
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
    case 1: problem = new ProblemMixed();
            break;
   case 2: problem = new Problem2a();
           break;
   case 3: problem = new Problem2b();
           break;
   case 4: problem = new Problem31();
           break;
   case 5: problem = new Problem32();
           break;
   case 6: problem = new Problem41();
           break;
   case 7: problem = new Problem42();
           break;
   case 8: problem = new Problem5();
           break;
   default: std::cerr << "Wrong problem number " << prob << std::endl;
            return 1;
  }

  if (gridname == "default" || gridname == "Default")
    gridname = problem->gridName();

  // initialize grid
  GridType grid ( gridname.c_str() ); // read the macro grid file
  // make initial grid
  grid.globalRefine( 4 );

  compute( grid, *problem, loops );

  delete problem;

  return 0;
}
