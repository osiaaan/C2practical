//Problem Interface and different Problems

/**
    @addtogroup Problem
    @{
*/

/** \brief Problem Interface to setup different problems
  */
typedef Dune::FieldMatrix<double,2,2> JacobianType;
struct Problem
{
    Problem (double mu, double lam)
    {
      Lame_mu = mu;
      Lame_lambda = lam;
    }

    /** \brief destructor */
    virtual ~Problem() {}

    /** \brief return value of exact solution
        \param[in] x  global coordinate to evaluate u

        \return u(x)
     */
    virtual GlobalCoordType u(const GlobalCoordType &x) const
    {
      return GlobalCoordType(0);
    }

    virtual JacobianType du(const GlobalCoordType &x) const
    {
      return JacobianType(0);
    }

    /** \brief return value of right hand side
        \param[in] x  global coordinate to evaluate f

        \return f(x)
    */
    virtual GlobalCoordType f(const GlobalCoordType &x) const
    {
      return GlobalCoordType(0);
    }

     /** \brief return value of Neumann boundary conditions;
                evaluated everywhere,
                make sure it vanishes on the Dirichlet boundary!
         \param[in] x  global coordinate to evaluate g

         \return g(x)
    */
    virtual GlobalCoordType h(const GlobalCoordType &x) const
    {
      return GlobalCoordType(0);
    }

    /** \brief return true if dirichlet data is to be used on full boundary. In case
     *         that false is returned Neumann zero boundary connditions are used - then
     *         \f$ \lambda \f$ should be greater than zero.
    */
    virtual bool useDirichlet() const
    {
      return true;
    }

    /** \brief return true if a boundary point is on the Dirichlet boundary
    */
    virtual bool isDirichlet(const GlobalCoordType &x) const
    {
      return useDirichlet();
    }

     /** \brief return value of Dirichlet boundary conditions
         \param[in] x  global coordinate to evaluate g

         \return ub(x)
    */
    virtual GlobalCoordType g(const GlobalCoordType &x) const
    {
      return u(x);
    }

    /** \brief return Lame constant mu
        \return \f$ \mu \f$
    */
    double get_mu() const
    {
      return Lame_mu;
    }

    /** \brief return Lame constant lambda
        \return \f$ \lambda \f$
    */
    double get_lambda() const
    {
      return Lame_lambda;
    }

    virtual const char* gridName() const
    {
      return nullptr;
    }

  protected:
    double Lame_mu;
    double Lame_lambda;
};

class ProblemMixed : public Problem {
  public:
    ProblemMixed()
    : Problem(1., 1.) {};
    virtual ~ProblemMixed() {}

    GlobalCoordType u(const GlobalCoordType &x) const
    {
      double factor_ = 4*M_PI;
      GlobalCoordType yo(0);
      yo[0] = sin(factor_*x[0]*x[1]);
      yo[1] = yo[0];
      return yo;
    }

    JacobianType du(const GlobalCoordType &x) const
    {
      double factor_ = 4*M_PI;
      JacobianType yo(0);

      double c = cos(factor_*x[0]*x[1]);
      yo[0][0] = factor_*c*x[1];
      yo[0][1] = factor_*c*x[0];
      yo[1][0] = yo[0][0];
      yo[1][1] = yo[0][1];
      return yo;
    }

    // set to -laplace u
    GlobalCoordType f(const GlobalCoordType &x) const
    {
      double factor_ = 4*M_PI;
      GlobalCoordType yo(0);
      double s = -sin(factor_*x[0]*x[1]);
      yo[0] = -factor_*factor_*s*(x[1]*x[1]+x[0]*x[0]);
      yo[1] = -factor_*factor_*s*(x[1]*x[1]+x[0]*x[0]);
      return yo;
    }
    double lambda() const
    {
      return 0;
    }
    bool useDirichlet() const override
    {
      return true;
    }
    virtual bool isDirichlet(const GlobalCoordType &x) const override
    {
      // neither the left nor the lower edge is Dirichlet
      // but both of the others are
      return (x[0] > 1e-12 && x[1] > 1e-12);
    }

    GlobalCoordType h(const GlobalCoordType &x) const override
        {
          JacobianType grad = du(x);
          GlobalCoordType ret(0);
          if (x[0] < 1e-12) {     // left boundary n=(-1,0)
            ret[0] = -grad[0][0];
            ret[1] = -grad[1][0];
            return ret;
          }
          else if (x[1] < 1e-12) { // lower boundary n=(0,-1)
            ret[0] = -grad[0][1];
            ret[1] = -grad[1][1];
            return ret;
          }
          else                   // should never get here!
            assert(0);
        }

    virtual const char* gridName() const
    {
      return "../problem/cube01.dgf";
    }
    //private:
    //double factor_;
};
/**
  @}
 */
