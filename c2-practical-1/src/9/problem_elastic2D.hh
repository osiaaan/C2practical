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
    Problem (double mu, double lam, double beta)
    {
      Lame_mu = mu;
      Lame_lambda = lam;
      beta_ = beta;
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

    double get_beta() const
    {
      return beta_;
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
    double beta_;
    double Lame_lambda;
};

class ProblemMixed : public Problem {
  public:
    ProblemMixed()
    : Problem(1., 1., 1.) {};
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

class Problem2a : public Problem {
  public:
    Problem2a()
    : Problem(10., 100., 0.) {};
    virtual ~Problem2a() {}

    GlobalCoordType u(const GlobalCoordType &x) const
    {
      GlobalCoordType yiw(0);
      yiw[0] = 0.0;
      yiw[1] = x[0]*x[0]*x[1]*x[1]*(1.-x[0])*(1.-x[0])*(1.-x[1])*(1.-x[1]);
      return yiw;
    }

    JacobianType du(const GlobalCoordType &x) const
    {
      JacobianType yo(0);

      yo[0][0] = 0;
      yo[0][1] = 0;
      yo[1][0] = x[1]*x[1]*(1-x[1])*(1-x[1])*2*x[0]*(2*x[0]*x[0]-3*x[0]+1);
      yo[1][1] = x[0]*x[0]*(1-x[0])*(1-x[0])*2*x[1]*(2*x[1]*x[1]-3*x[1]+1);;
      return yo;
    }

    // set to -laplace u
    GlobalCoordType f(const GlobalCoordType &x) const
    {
      GlobalCoordType yo(0);

      yo[0] = -400*x[1]*(2*x[1]*x[1] -3*x[1] +1)*x[0]*(2*x[0]*x[0] -3*x[0] +1) -
                40*x[1]*(2*x[1]*x[1] -3*x[1] +1)*x[0]*(2*x[0]*x[0] -3*x[0] +1);
      yo[1] = -10*x[1]*x[1]*(1.-x[1])*(1.-x[1]) * 2*(6*x[0]*x[0] -6*x[0] + 1) -
                (2*10 + 100) * x[0]*x[0]*(1.-x[0])*(1.-x[0]) * 2*(6*x[1]*x[1] -6*x[1] + 1);
      return yo;
    }
    bool useDirichlet() const override
    {
      return true;
    }
    virtual bool isDirichlet(const GlobalCoordType &x) const override
    {
      // neither the left nor the lower edge is Dirichlet
      // but both of the others are
      return true;//(x[0] > 1e-12 && x[1] > 1e-12);
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

class Problem2b : public Problem {
  public:
    Problem2b()
    : Problem(200., 10., 0.) {};
    virtual ~Problem2b() {}

    double const mu = 200.0;
    double const lambda = 10.0;
    GlobalCoordType u(const GlobalCoordType &x) const
    {
      GlobalCoordType yiw(0);
      yiw[0] = 0.0;
      yiw[1] = (x[0]*x[1])/10.0;
      return yiw;
    }

    JacobianType du(const GlobalCoordType &x) const
    {
      JacobianType yo(0);

      yo[0][0] = 0;
      yo[0][1] = 0;
      yo[1][0] = x[1]/10.0;
      yo[1][1] = x[0]/10.0;
      return yo;
    }

    // set to -laplace u
    GlobalCoordType f(const GlobalCoordType &x) const
    {
      GlobalCoordType yo(0);

      yo[0] =  -(mu + lambda)/10.0;
      yo[1] = 0.0;
      return yo;
    }
    bool useDirichlet() const override
    {
      return true;
    }
    virtual bool isDirichlet(const GlobalCoordType &x) const override
    {
      // neither the left nor the lower edge is Dirichlet
      // but both of the others are
      return (x[0] < 1e-12);
    }

    GlobalCoordType h(const GlobalCoordType &x) const override
        {
          JacobianType grad = du(x);
          GlobalCoordType ret(0);
          if (1 - x[1] < 1e-12) {     // upper boundary n=(0,1)
            ret[0] = mu*(grad[0][1] + grad[1][0]);
            ret[1] = (2*mu + lambda)*grad[1][1] + lambda*grad[0][0];
            return ret;
          }
          else if (1 - x[0] < 1e-12) { // right boundary n=(1,0)
            ret[0] = (2*mu + lambda)*grad[0][0] + lambda*grad[1][1];
            ret[1] = mu*(grad[0][1] + grad[1][0]);
            return ret;
          }
          else if (x[1] < 1e-12) { // lower boundary n=(0,-1)
            ret[0] = -mu*(grad[0][1] + grad[1][0]);
            ret[1] = -(2*mu + lambda)*grad[1][1] - lambda*grad[0][0];
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

class Problem31 : public Problem {
  public:
    Problem31()
    : Problem(1., 1., 0.) {}

    virtual ~Problem31() {}

    double h(const double &x) const
    {
      return x*x*(1-x)*(1-x);
    }

    double dh(const double &x) const
    {
      return 2*x*(1-x)*(1-x) - 2*x*x*(1-x);
    }

    double ddh(const double &x) const
    {
      return 2*(1-x)*(1-x) - 4*x*(1-x) - 4*x*(1-x) + 2*x*x;
    }

    double dddh(const double &x) const
    {
      return 12*(2*x-1);
    }

    GlobalCoordType u(const GlobalCoordType &x) const
    {
      GlobalCoordType u(0);
      u[0] = -100*h(x[0])*dh(x[1]);
      u[1] = 100*h(x[1])*dh(x[0]);
      return u;
    }

    JacobianType du(const GlobalCoordType &x) const
    {
      JacobianType du(0);
      du[0][0] = -100*dh(x[0])*dh(x[1]);
      du[0][1] = -100*h(x[0])*ddh(x[1]);
      du[1][0] = 100*h(x[1])*ddh(x[0]);
      du[1][1] = 100*dh(x[1])*dh(x[0]);
      return du;
    }

    // set to -laplace u
    GlobalCoordType f(const GlobalCoordType &x) const
    {
      GlobalCoordType ret;
      ret[0] = -100*Lame_mu*(-ddh(x[0])*dh(x[1])-h(x[0])*dddh(x[1]));
      ret[1] = -100*Lame_mu*(ddh(x[1])*dh(x[0])+h(x[1])*dddh(x[0]));
      return ret;
    }

    GlobalCoordType g(const GlobalCoordType &x) const
    {
      return u(x);
    }

    // double lambda() const
    // {
    //   return 0;
    // }

    bool useDirichlet() const override
    {
      return true;
    }

    virtual bool isDirichlet(const GlobalCoordType &x) const override
    {
      // neither the left nor the lower edge is Dirichlet
      // but both of the others are
      return true;
    }


    virtual const char* gridName() const
    {
      return "../problem/cube01.dgf";
    }

};

class Problem32 : public Problem {
  public:
    Problem32()
    : Problem(1., 1e8, 0.) {}

    virtual ~Problem32() {}

    double h(const double &x) const
    {
      return x*x*(1-x)*(1-x);
    }

    double dh(const double &x) const
    {
      return 2*x*(1-x)*(1-x) - 2*x*x*(1-x);
    }

    double ddh(const double &x) const
    {
      return 2*(1-x)*(1-x) - 4*x*(1-x) - 4*x*(1-x) + 2*x*x;
    }

    double dddh(const double &x) const
    {
      return 12*(2*x-1);
    }

    GlobalCoordType u(const GlobalCoordType &x) const
    {
      GlobalCoordType u(0);
      u[0] = -100*h(x[0])*dh(x[1]);
      u[1] = 100*h(x[1])*dh(x[0]);
      return u;
    }

    JacobianType du(const GlobalCoordType &x) const
    {
      JacobianType du(0);
      du[0][0] = -100*dh(x[0])*dh(x[1]);
      du[0][1] = -100*h(x[0])*ddh(x[1]);
      du[1][0] = 100*h(x[1])*ddh(x[0]);
      du[1][1] = 100*dh(x[1])*dh(x[0]);
      return du;
    }

    // set to -laplace u
    GlobalCoordType f(const GlobalCoordType &x) const
    {
      GlobalCoordType ret;
      ret[0] = -100*Lame_mu*(-ddh(x[0])*dh(x[1])-h(x[0])*dddh(x[1]));
      ret[1] = -100*Lame_mu*(ddh(x[1])*dh(x[0])+h(x[1])*dddh(x[0]));
      return ret;
    }

    GlobalCoordType g(const GlobalCoordType &x) const
    {
      return u(x);
    }

    // double lambda() const
    // {
    //   return 0;
    // }

    bool useDirichlet() const override
    {
      return true;
    }

    virtual bool isDirichlet(const GlobalCoordType &x) const override
    {
      // neither the left nor the lower edge is Dirichlet
      // but both of the others are
      return true;
    }


    virtual const char* gridName() const
    {
      return "../problem/cube01.dgf";
    }

};

class Problem41 : public Problem {
public:
    Problem41()
    : Problem(1., 1., 1.) {}

    virtual ~Problem41() {}




    double h(const double &x) const
    {
      return x*x*(1-x)*(1-x);
    }

    double dh(const double &x) const
    {
      return 2*x*(1-x)*(1-x) - 2*x*x*(1-x);
    }

    double ddh(const double &x) const
    {
      return 2*(1-x)*(1-x) - 4*x*(1-x) - 4*x*(1-x) + 2*x*x;
    }

    double dddh(const double &x) const
    {
      return 12*(2*x-1);
    }

    GlobalCoordType u(const GlobalCoordType &x) const
    {
      GlobalCoordType u(0);
      u[0] = -100*h(x[0])*dh(x[1]);
      u[1] = 100*h(x[1])*dh(x[0]);
      return u;
    }

    JacobianType du(const GlobalCoordType &x) const
    {
      JacobianType du(0);
      du[0][0] = -100*dh(x[0])*dh(x[1]);
      du[0][1] = -100*h(x[0])*ddh(x[1]);
      du[1][0] = 100*h(x[1])*ddh(x[0]);
      du[1][1] = 100*dh(x[1])*dh(x[0]);
      return du;
    }

    // set to -laplace u
    GlobalCoordType f(const GlobalCoordType &x) const
    {
      GlobalCoordType ret;
      ret[0] = -100*Lame_mu*(-ddh(x[0])*dh(x[1])-h(x[0])*dddh(x[1]));
      ret[1] = -100*Lame_mu*(ddh(x[1])*dh(x[0])+h(x[1])*dddh(x[0]));
      return ret;
    }

    GlobalCoordType g(const GlobalCoordType &x) const
    {
      return u(x);
    }

    // double lambda() const
    // {
    //   return 0;
    // }

    bool useDirichlet() const override
    {
      return true;
    }

    virtual bool isDirichlet(const GlobalCoordType &x) const override
    {
      // neither the left nor the lower edge is Dirichlet
      // but both of the others are
      return true;
    }


    virtual const char* gridName() const
    {
      return "../problem/cube01.dgf";
    }

};

class Problem42 : public Problem {
  public:
    Problem42()
    : Problem(1., 1., 0.0) {
    }

    virtual ~Problem42() {}

    double get_beta() const
    {
      return 0.0;
    }

    double h(const double &x) const
    {
      return x*x*(1-x)*(1-x);
    }

    double dh(const double &x) const
    {
      return 2*x*(1-x)*(1-x) - 2*x*x*(1-x);
    }

    double ddh(const double &x) const
    {
      return 2*(1-x)*(1-x) - 4*x*(1-x) - 4*x*(1-x) + 2*x*x;
    }

    double dddh(const double &x) const
    {
      return 12*(2*x-1);
    }

    GlobalCoordType u(const GlobalCoordType &x) const
    {
      GlobalCoordType u(0);
      u[0] = -100*h(x[0])*dh(x[1]);
      u[1] = 100*h(x[1])*dh(x[0]);
      return u;
    }

    JacobianType du(const GlobalCoordType &x) const
    {
      JacobianType du(0);
      du[0][0] = -100*dh(x[0])*dh(x[1]);
      du[0][1] = -100*h(x[0])*ddh(x[1]);
      du[1][0] = 100*h(x[1])*ddh(x[0]);
      du[1][1] = 100*dh(x[1])*dh(x[0]);
      return du;
    }

    // set to -laplace u
    GlobalCoordType f(const GlobalCoordType &x) const
    {
      GlobalCoordType ret;
      ret[0] = -100*Lame_mu*(-ddh(x[0])*dh(x[1])-h(x[0])*dddh(x[1]));
      ret[1] = -100*Lame_mu*(ddh(x[1])*dh(x[0])+h(x[1])*dddh(x[0]));
      return ret;
    }


    GlobalCoordType g(const GlobalCoordType &x) const
    {
      return u(x);
    }

    // double lambda() const
    // {
    //   return 0;
    // }

    bool useDirichlet() const override
    {
      return true;
    }

    virtual bool isDirichlet(const GlobalCoordType &x) const override
    {
      // neither the left nor the lower edge is Dirichlet
      // but both of the others are
      return true;
    }


    virtual const char* gridName() const
    {
      return "../problem/cube01.dgf";
    }

};

class Problem5 : public Problem {
  public:
    Problem5()
    : Problem(10., 10., 0.) {};
    virtual ~Problem5() {}

    GlobalCoordType u(const GlobalCoordType &x) const
    {
      GlobalCoordType yiw(0);
      yiw[0] = 0.0;
      yiw[1] = (x[0]*x[1])/10.0;
      return yiw;
    }

    JacobianType du(const GlobalCoordType &x) const
    {
      JacobianType yo(0);

      yo[0][0] = 0;
      yo[0][1] = 0;
      yo[1][0] = x[1]/10.0;
      yo[1][1] = x[0]/10.0;
      return yo;
    }

    // set to -laplace u
    GlobalCoordType f(const GlobalCoordType &x) const
    {
      GlobalCoordType yo(0);

      yo[0] = 0.0;
      yo[1] = 0.009;
      return yo;
    }
    bool useDirichlet() const override
    {
      return true;
    }
    virtual bool isDirichlet(const GlobalCoordType &x) const override
    {
      // neither the left nor the lower edge is Dirichlet
      // but both of the others are
      return (x[0] < 1e-12);
    }

    GlobalCoordType h(const GlobalCoordType &x) const override
        {
          JacobianType grad = du(x);
          GlobalCoordType ret(0);
          return ret;
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
