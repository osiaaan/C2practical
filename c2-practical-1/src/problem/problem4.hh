#include "problem3.hh"

/**
    @addtogroup Problem
    @{
*/

/** \brief C0 problem
 */

/** \brief C1 problem
 */
class ProblemC1 : public Problem {

  public:
    ProblemC1() {}
    virtual ~ProblemC1() {}

    double u(const GlobalCoordType &x) const
    {
      double r2 = x.two_norm2();
      return (r2<0.5)?cos(r2*2.*M_PI)+1.:0;
    }

    // not needed here
    void du(const GlobalCoordType &x,GlobalCoordType &du) const
    {
      double r2 = x.two_norm2();
      double factor = (r2<0.5)?-4.*M_PI*sin(r2*2.*M_PI):0;
      du[0] = factor*x[0];
      du[1] = factor*x[1];
    }

    // not needed here
    double f(const GlobalCoordType &x) const
    {
      return 0.;
    }
};

/** \brief Sqrt Problem for Interpolation
 */
class ProblemSqrt : public Problem {

  public:
    ProblemSqrt() {}
    virtual ~ProblemSqrt() {}

    double u(const GlobalCoordType &x) const
    {
      return sqrt(x[0] * x[1]);
    }

    // not needed here
    void du(const GlobalCoordType &x,GlobalCoordType &du) const
    {
      double factor = 1./(2.*sqrt(x[0]*x[1]));
      du[0] = factor*x[0];
      du[1] = factor*x[1];
    }

    // not needed here
    double f(const GlobalCoordType &x) const
    {
      return (x[0]*x[0]+x[1]*x[1])/(4.*pow(x[0]*x[1],3./2.));
    }
};


class Problem2019_Ass1_1 : public Problem {
  public:
    Problem2019_Ass1_1() {}
    virtual ~Problem2019_Ass1_1() {}

    double u(const GlobalCoordType &x) const
    {
      double q = M_PI*x[0]/(0.25+x[0]*x[1]);
      return sin(q);
    }

    // not needed here
    void du(const GlobalCoordType &x,GlobalCoordType &du) const
    {
      double q    = M_PI*x[0]/(0.25+x[0]*x[1]);
      double dqdx = M_PI*( (0.25+x[0]*x[1]) - x[0]*x[1] ) / pow(0.25+x[0]*x[1],2);
      double dqdy = M_PI*( -x[0]*x[0] ) / pow(0.25+x[0]*x[1],2);
      du[0] = cos(q)*dqdx;
      du[1] = cos(q)*dqdy;
    }

    // not needed here
    double f(const GlobalCoordType &x) const
    {
      assert(0);
    }
};

class Problem2019_Ass1_2_3 : public Problem {
private:
  const double a_;
public:
  Problem2019_Ass1_2_3(double a) : a_(a) {}
  virtual ~Problem2019_Ass1_2_3() {}

  double u(const GlobalCoordType &x) const {
    if (x[0]==a_ && x[1]==a_) {
      return 0;
    }
    else
    {
      return (pow(x[0]-a_,3)-(x[0]-a_)*pow(x[1]-a_,2)) / (pow(x[0]-a_,2)+pow(x[1]-a_,2));
    }
  }

  void du(const GlobalCoordType &x, GlobalCoordType &du) const {
    if(x[0]==a_ && x[1]==a_) {
      du[0] = du[1] = 0.0;
    }
    else
    {
      double den = pow(pow(x[0]-a_,2)+pow(x[1]-a_,2),2);
      du[0] = (pow(x[0]-a_,4)+4.0*pow(x[0]-a_,2)*pow(x[1]-a_,2)-pow(x[1]-a_,4)) / den;
      du[1] = (-4.0*pow(x[0]-a_,3)*(x[1]-a_)) / den;
    }
  }

  double f(const GlobalCoordType &x) const
  {
    return 0;
  }

};

/** \brief Linear Problem for Interpolation
 */
class Problem2019_Ass1_4 : public Problem {

  public:
    Problem2019_Ass1_4() {}
    virtual ~Problem2019_Ass1_4() {}

    double u(const GlobalCoordType &x) const
    {
      return sqrt(2.)*x[0] + sqrt(3.)*x[1];
    }

    // not needed here
    void du(const GlobalCoordType &x,GlobalCoordType &du) const
    {
      du[0] = sqrt(2.);
      du[1] = sqrt(3.);
    }

    // not needed here
    double f(const GlobalCoordType &x) const
    {
      return 0.0 ;
    }
};

/** \brief H^2\H^1 Problem for Interpolation
 */
class Problem2019_Ass1_5 : public Problem {

  double a_;
  public:
    Problem2019_Ass1_5(double a) : a_(a) {}
    virtual ~Problem2019_Ass1_5() {}

    double u(const GlobalCoordType &x) const
    {
      return (x[0]>a_)? 0:sin(10.*(x[0]-a_));
    }

    // not needed here
    void du(const GlobalCoordType &x,GlobalCoordType &du) const
    {
      du[0] = (x[0]>a_)? 0:10.*cos(-10.*(x[0]-a_));
      du[1] = 0.;
    }

    // not needed here
    double f(const GlobalCoordType &x) const
    {
      return 0.0 ;
    }
};

class Problem2_1 : public Problem {

  double a_;
  public:
    Problem2_1(double a) : a_(a) {}
    virtual ~Problem2_1() {}


};


/**
double u(const GlobalCoordType &x) const
{
  double r2 = x.two_norm2();
  return (r2<0.5)?cos(r2*M_PI):0;
}

// not needed here
void du(const GlobalCoordType &x,GlobalCoordType &du) const
{
  double r2 = x.two_norm2();
  double factor = (r2<0.5)?-2.*M_PI*sin(r2*M_PI):0;
  du[0] = factor*x[0];
  du[1] = factor*x[1];
}

// not needed here
double f(const GlobalCoordType &x) const
{
  return 0.;
}
  @}
 */
