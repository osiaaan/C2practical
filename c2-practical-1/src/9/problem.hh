//Problem Interface and diffrent Problems

/**
    @addtogroup Problem
    @{
*/

/** \brief Problem Interface to setup different problems
  */
struct Problem 
{
    /** \brief destructor */
    virtual ~Problem() {}
    
    /** \brief return value of exact solution
        \param[in] x  global coordinate to evaluate u 

        \return u(x) 
     */
    virtual double u(const GlobalCoordType &x) const = 0;

    /** \brief return value of right hand side 
        \param[in] x  global coordinate to evaluate f 

        \return f(x)
    */
    virtual double f(const GlobalCoordType &x) const = 0;

     /** \brief return value of Dirichlet boundary conditions 
                (default implementation return u(x) ) 
         \param[in] x  global coordinate to evaluate g

         \return g(x)
    */
    virtual double g(const GlobalCoordType &x) const 
    {
      return u(x);
    }

    //! NEW
    /** \brief return value for Neuman boundary conditions
                (default implementation return 0 ) 
         \param[in] x  global coordinate to evaluate g

         \return h(x)
    */
    virtual double h(const GlobalCoordType &x) const 
    {
      return 0;
    }

    /** \brief return derivative of u  
        \param[in] x  global coordinate to evaluate du 
        \param[out] du the gradient of the function u
    */
    virtual void du(const GlobalCoordType &x, GlobalCoordType &du) const = 0;

    /** \brief return constant in front of zero order term

        \return \f$ \lambda \f$
    */
    virtual double lambda() const 
    {
      return 0.0;
    }

    /** \brief return true if dirichlet data is to be used on full boundary. In case
     *         that false is returned Neumann zero boundary connditions are used - then
     *         \f$ \lambda \f$ should be greater than zero.
    */
    virtual bool useDirichlet() const
    {
      return true;
    }

    //! NEW
    /** \brief return if point x is on Dirichlet boundary
        Default implementation is true if useDirichlet() returns true,
        false otherwise
        \param[in] x  global coordinate
        \param[out] true if x is on Dirichlet boundary
    */
    virtual bool isDirichlet(const GlobalCoordType &x) const
    {
      return useDirichlet();
    }

    //! NEW
    /** \brief return the name of the grid file to use for this
        problem. If a grid file name is passed in on command line then that overrides this
        setting - pass 'default' on the command line to use the default file for this problem.
    */
    virtual const char* gridName() const
    {
      return nullptr;
    }
};

/** \brief simple problem using mixed BCs
 */
class ProblemMixed : public Problem {
  public:
    ProblemMixed()
    : factor_(4*M_PI) {}
    virtual ~ProblemMixed() {}

    double u(const GlobalCoordType &x) const
    {
      return sin(factor_*x[0]*x[1]);
    }

    void du(const GlobalCoordType &x,GlobalCoordType &du) const
    {
      double c = cos(factor_*x[0]*x[1]);
      du[0] = factor_*c*x[1];
      du[1] = factor_*c*x[0];
    }

    // set to -laplace u
    double f(const GlobalCoordType &x) const
    {
      double s = -sin(factor_*x[0]*x[1]);
      return -factor_*factor_*s*(x[1]*x[1]+x[0]*x[0]);
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
      // neither the right nor the lower edge is Dirichlet
      // but both of the others are
      return (x[0] > 1e-12 && x[1] > 1e-12);
    }
    virtual double h(const GlobalCoordType &x) const override
    {
      GlobalCoordType grad;
      du(x,grad);
      if (x[0] < 1e-12)      // right boundary n=(-1,0)
        return -grad[0];
      else if (x[1] < 1e-12) // lower boundary n=(0,-1)
        return -grad[1];
      else                   // should never get here!
        assert(0);
    }
    virtual const char* gridName() const
    {
      return "../problem/cube01.dgf";
    }
    private:
    double factor_;
};

/**
  @}
 */

