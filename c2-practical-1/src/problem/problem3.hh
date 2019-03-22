#include "problem.hh"

/**
    @addtogroup Problem
    @{
*/
////////////////////////////////////////////////////
//
//  Assignment 1:
//    - implement different test cases
//    - only methods u and du are required but the
//      method f has to be implemented (can return 0)
//
////////////////////////////////////////////////////
/** \brief Sinus Problem for Interpolation
 */
class ProblemSin : public Problem {

  public:
    ProblemSin() {}
    virtual ~ProblemSin() {}

    double u(const GlobalCoordType &x) const
    {
      return sin(2.*M_PI*x[0]*x[1]);
    }

    // not needed here
    void du(const GlobalCoordType &x,GlobalCoordType &du) const
    {
      double q    = 2.*M_PI*x[0]*x[1];
      double dqdx = 2.*M_PI*x[1];
      double dqdy = 2.*M_PI*x[0];
      du[0] = cos(q)*dqdx;
      du[1] = cos(q)*dqdy;
    }

    // not needed here
    double f(const GlobalCoordType &x) const
    {
      assert(0);
    }
};

class Problem1 : public Problem {

public:
  Problem1() {}
   virtual ~Problem1() {}

    double u(const GlobalCoordType &x) const override
   {
       double q = (M_PI*x[0])/(0.25+(x[0]*x[1]));
        return sin(q);
    }

    void du(const GlobalCoordType &x, GlobalCoordType &du) const
    {
      double q = (M_PI*x[0])/(0.25+x[0]*x[1]);
      double ddx = ((M_PI/4)*cos(q))/((0.25+x[0]*x[1])*(0.25+x[0]*x[1]));
      double ddy = cos(q)/((0.25+x[0]*x[1])*(0.25+x[0]*x[1]));
      du[0] = ddx;
      du[1] = -M_PI*x[0]*x[0]*ddy;
    }

    double f(const GlobalCoordType &x) const
    {
       assert(0);
    }
};

class Problem2 : public Problem {

public:
    Problem2() {}
    virtual ~Problem2() {}

    const double a = 1.0/M_PI;
    double u(const GlobalCoordType &x) const override
    {
      double num = (x[0]-a)*(x[0]-a)*(x[0]-a)-(x[0]-a)*(x[1]-a)*(x[1]-a);
      double den = (x[0]-a)*(x[0]-a)+(x[1]-a)*(x[1]-a);
      return num/den;
    }

    void du(const GlobalCoordType &x, GlobalCoordType &du) const
    {
        double num1 = -2 * (a-x[0]) * ((x[0]-a)*(x[1]-a)*(x[1]-a) - (x[0]-a)*(x[0]-a)*(x[0]-a)) +
        ((x[0]-a)*(x[0]-a) + (x[1]-a)*(x[1]-a))*(3*(x[0]-a)*(x[0]-a) - (x[1]-a)*(x[1]-a));
        double den1 = ( (a-x[0])*(a-x[0]) + (a-x[1])*(a-x[1]) );
        du[0] = num1/(den1*den1);

        double num2 = -4.*(x[0]-a)*(x[0]-a)*(x[0]-a)*(x[1]-a);
        du[1] = num2/(den1*den1);
    }

    double f(const GlobalCoordType &x) const
    {
       assert(0);
    }
};

class Problem3 : public Problem {

public:
    Problem3() {}
    virtual ~Problem3() {}

    const double a = 1.0/2.0;
    double u(const GlobalCoordType &x) const override
    {
      // Excluding the singularity
      if (x[0]==0.5 && x[1]==0.5)
      {
          return 0;
      }
        double num = (x[0]-a)*(x[0]-a)*(x[0]-a)-(x[0]-a)*(x[1]-a)*(x[1]-a);
        double den = (x[0]-a)*(x[0]-a)+(x[1]-a)*(x[1]-a);
        return num/den;
    }

    void du(const GlobalCoordType &x, GlobalCoordType &du) const
    {
        double num1 = -2 * (a-x[0]) * ((x[0]-a)*(x[1]-a)*(x[1]-a) - (x[0]-a)*(x[0]-a)*(x[0]-a)) +
        ((x[0]-a)*(x[0]-a) + (x[1]-a)*(x[1]-a))*(3*(x[0]-a)*(x[0]-a) - (x[1]-a)*(x[1]-a));
        double den1 = ( (a-x[0])*(a-x[0]) + (a-x[1])*(a-x[1]) );
        du[0] = num1/(den1*den1);

        double num2 = -4.*(x[0]-a)*(x[0]-a)*(x[0]-a)*(x[1]-a);
        du[1] = num2/(den1*den1);
    }

    double f(const GlobalCoordType &x) const
    {
       assert(0);
    }
};


class Problem4 : public Problem {

  public:
    Problem4() {}
    virtual ~Problem4() {}

    double u(const GlobalCoordType &x) const
    {
      return sqrt(2)*x[0] + sqrt(3)*x[1];
    }

    void du(const GlobalCoordType &x,GlobalCoordType &du) const
    {
      du[0] = sqrt(2);
      du[1] = sqrt(3);
    }

    // not needed here
    double f(const GlobalCoordType &x) const
    {
      assert(0);
    }
};



class Problem5 : public Problem {

  public:
    Problem5() {}
    virtual ~Problem5() {}

    double u(const GlobalCoordType &x) const
    {
      return pow(fabs(x[0]),1.0/2.0);
      //return pow(fabs(x[0]),1.0/2.0);
    }

    void du(const GlobalCoordType &x,GlobalCoordType &du) const
    {
      du[0] = (1.0/2.0)*pow(fabs(x[0]),-0.5);
      du[1] = 0;
    }

    // not needed here
    double f(const GlobalCoordType &x) const
    {
      assert(0);
    }
};


/** \brief Problem for exercise 3 */



/**
  @}
 */
