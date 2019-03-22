#include "problem4.hh"

/**
    @addtogroup Problem
    @{
*/

/** \brief C0 problem
 */


 class Problem2019_Ass2_1 : public Problem {

   public:
     Problem2019_Ass2_1() {}
     virtual ~Problem2019_Ass2_1() {}


     double u(const GlobalCoordType &x) const
     {
       return cos( 16.0 * M_PI * x[0] * x[1] * ( 1 - x[0] ) * ( 1 - x[1] ) );
     }

     // not needed here
     void du(const GlobalCoordType &x,GlobalCoordType &du) const
     {
       du[0] = -16 * M_PI * (2*x[0] -1) * (x[1] - 1) * x[1] * sin(16 * M_PI * (x[0] - 1) * x[0] * ( x[1] - 1) * x[1] );
       du[1] = -16 * M_PI * (2*x[1] -1) * (x[0] - 1) * x[0] * sin(16 * M_PI * (x[0] - 1) * x[0] * ( x[1] - 1) * x[1] );
     }

     double f(const GlobalCoordType &x) const
     {
       double t0_y = -32 * M_PI * ( x[0] - 1 ) * x[0] * sin( 16 * M_PI * ( x[0] - 1 ) * x[0] * ( x[1] - 1 ) * x[1] );
       double t1_y = 16 * M_PI * ( x[0] - 1 ) * x[0] * ( 2 * x[1] - 1 );
       double t2_y = 16 * M_PI * ( x[0] - 1 ) * x[0] * ( x[1] -1 ) + 16 *M_PI * ( x[0] - 1 ) * x[0] * x[1];
       double t3_y = cos( 16 * M_PI * ( x[0] - 1 ) * x[0] * ( x[1] - 1 ) * x[1] );
       double derivative2_y = t0_y - (t1_y * t2_y * t3_y);

       double t0_x = -32 * M_PI * (x[1] - 1) * x[1] * sin(16 * M_PI * (x[0] - 1) * x[0] * ( x[1] - 1) * x[1] );
       double t1_x = 16 * M_PI * (x[1] -1) * x[1] * (2 * x[0] - 1);
       double t2_x = 16 * M_PI * (x[0] - 1) * x[1] * ( x[1] -1 ) + 16 *M_PI * ( x[1] - 1 ) * x[0] * x[1];
       double t3_x = cos( 16 * M_PI * (x[0] - 1) * x[0] * (x[1] - 1) * x[1] );
       double derivative2_x = t0_x - (t1_x * t2_x * t3_x);

       return u(x) - derivative2_x - derivative2_y ;
     }

       double lambda() const
     {
         return 1;
     }

     bool useDirichlet() const
     {
         return false;
     }

 };

 class Problem2019_Ass2_2 : public Problem {

   public:
     Problem2019_Ass2_2() {}
     virtual ~Problem2019_Ass2_2() {}


     double u(const GlobalCoordType &x) const
     {
       return cos( 16.0 * M_PI * x[0] * x[1] * ( 1 - x[0] ) * ( 1 - x[1] ) );
     }

     // not needed here
     void du(const GlobalCoordType &x,GlobalCoordType &du) const
     {
       du[0] = -16 * M_PI * (2*x[0] -1) * (x[1] - 1) * x[1] * sin(16 * M_PI * (x[0] - 1) * x[0] * ( x[1] - 1) * x[1] );
       du[1] = -16 * M_PI * (2*x[1] -1) * (x[0] - 1) * x[0] * sin(16 * M_PI * (x[0] - 1) * x[0] * ( x[1] - 1) * x[1] );
     }

     double f(const GlobalCoordType &x) const
     {
       double t0_y = -32 * M_PI * ( x[0] - 1 ) * x[0] * sin( 16 * M_PI * ( x[0] - 1 ) * x[0] * ( x[1] - 1 ) * x[1] );
       double t1_y = 16 * M_PI * ( x[0] - 1 ) * x[0] * ( 2 * x[1] - 1 );
       double t2_y = 16 * M_PI * ( x[0] - 1 ) * x[0] * ( x[1] -1 ) + 16 *M_PI * ( x[0] - 1 ) * x[0] * x[1];
       double t3_y = cos( 16 * M_PI * ( x[0] - 1 ) * x[0] * ( x[1] - 1 ) * x[1] );
       double derivative2_y = t0_y - (t1_y * t2_y * t3_y);

       double t0_x = -32 * M_PI * (x[1] - 1) * x[1] * sin(16 * M_PI * (x[0] - 1) * x[0] * ( x[1] - 1) * x[1] );
       double t1_x = 16 * M_PI * (x[1] -1) * x[1] * (2 * x[0] - 1);
       double t2_x = 16 * M_PI * (x[0] - 1) * x[1] * ( x[1] -1 ) + 16 *M_PI * ( x[1] - 1 ) * x[0] * x[1];
       double t3_x = cos( 16 * M_PI * (x[0] - 1) * x[0] * (x[1] - 1) * x[1] );
       double derivative2_x = t0_x - (t1_x * t2_x * t3_x);

       return u(x) - derivative2_x - derivative2_y ;
     }

       double lambda() const
     {
         return 1;
     }

     bool useDirichlet() const
     {
         return true;
     }

 };

 class Problem2019_Ass2_3 : public Problem {

   public:
     Problem2019_Ass2_3() {}
     virtual ~Problem2019_Ass2_3() {}

     const double b = 0.6;
     const double a = 0.9;

     double u(const GlobalCoordType &x) const
     {
       if( x[1]  <= b*(x[0]+1) )
       {
         return x[0]*x[1];
       }
       else
       {
         return x[0]*x[1] + pow( x[1] - b*( x[0] + 1 ) , a );
       }
     }

     // not needed here
     void du(const GlobalCoordType &x,GlobalCoordType &du) const
     {
       if( x[1]  <= b*(x[0]+1) )
       {
        du[0] = x[1];
        du[1] = x[0];
       }
       else
       {
         du[0] = x[1] - a*b*pow( x[1] - b*( x[0] + 1 ) , a-1 );
         du[1] = x[0] + a*pow( x[1] - b*( x[0] + 1 ) , a-1 );
       }

     }

     double f(const GlobalCoordType &x) const
     {
       return u(x);
     }

       double lambda() const
     {
         return 1;
     }

     bool useDirichlet() const
     {
         return true;
     }

     bool Q4() const
     {
       return true;
     }
 };
