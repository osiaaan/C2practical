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

    /** \brief return constant in front of advection term

        \return \f$ b \f$
    */
    virtual GlobalCoordType b() const
    {
      return GlobalCoordType(0);
    }

    /** \brief return true if dirichlet data is to be used on full boundary. In case
     *         that false is returned Neumann zero boundary connditions are used - then
     *         \f$ \lambda \f$ should be greater than zero.
    */
    virtual bool useDirichlet() const
    {
      return false;
    }

    //Q4
    virtual bool Q4() const
    {
      return false;
    }
};

/**
  @}
 */
