#include <fstream>
#include <iostream>
#include <stdexcept>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <vector>
#include "headers/tm.h"

using namespace std;


const double SAFETY = 0.9;
const double PGROW = -0.2;
const double PSHRNK = -0.25;
const double ERRCON = 1.89e-4;
const int MAXSTP = 30000;
const double TINY = 1.0e-30;


void derivs0(t1D, t1D, t1D, t1D &);
void derivs1(t1D, t1D, t1D, t1D &);
void derivs2(t1D, t1D, t1D, t1D &);
void derivs3(t1D, t1D, t1D, t1D &);
void initMemory(int, int, t1D &, t1D &, t1D &, t1D &, t1D &, t1D &, t1D &, t1D &, t2D &);
void odeint(t1D, t1D, t2D &, t1D &, t1D &, int, int &, double, double,
            void (*)(t1D, t1D, t1D, t1D &),
            double = 0.0, double = 1.0e-7, double = 1.0e-3);  //default value of eps = 1.0e-7
void setOutputData(string, t2D, t1D, t1D, t1D, t1D, t1D, int, int);
void readInputData(int, int, ifstream &, t1D &, t1D &);

//---------------------------------------------------------
string outputFileC  = "coverage.dat";
string inputFile    = "inputGrowthData.dat";

//---------------------------------------------------------
int main()
{
   void (*derivs[])(t1D, t1D, t1D, t1D &) = {derivs0, derivs1, derivs2, derivs3};
   t1D thetaStart, intensity,  sqrDelta, growthTime, theta, dThetaDt, C, gR;
   t2D coverage;
   int numIntervals, numLayers, numReturn, model;
   double tMax, dXsav;
   string temp;

   ifstream fileTh(Str(inputFile), ios::in);
   if(!fileTh) {
      cout << "Input file: " <<inputFile<<" not found";
      cin.get();
      exit(EXIT_FAILURE);
   }
   fileTh >> model;           //growth model
   if((model < 0 ) || (model > 3)) {
      cout << "Growth model not exists";
      cin.get();
      exit(EXIT_FAILURE);
   }
   fileTh >> temp;            //ignores a text line from inputGrowthData file
   fileTh >> numLayers;       //number of growing layers
   fileTh >> tMax;            //upper limit of growth time interval
   fileTh >> numIntervals;    //number of integration intervals
   fileTh >> temp;

   initMemory(numLayers, numIntervals, thetaStart, intensity,  sqrDelta, growthTime,
              theta, dThetaDt, C, gR, coverage);
   readInputData(numLayers, model, fileTh, C, gR);

   fileTh.close();

   dXsav = tMax/numIntervals;
   for (int n=1; n<=numLayers; n++)
      thetaStart[n] = 0.0;

   odeint(C, gR, coverage, thetaStart, growthTime, numIntervals, numReturn, dXsav, tMax, derivs[model]);

   cout << "numReturn = " << numReturn << endl;

   setOutputData(outputFileC, coverage, growthTime, thetaStart, intensity, sqrDelta, gR, numLayers, numReturn);

   cout << "End of calculations...";
   //cin.get();
   return 0;
}
//---------------------------------------------------------
void  initMemory(int numLayers, int numIntervals, t1D &thetaStart, t1D &intensity, t1D &sqrDelta,
                 t1D &growthTime, t1D &theta, t1D &dThetaDt, t1D &C, t1D &gR, t2D &coverage)
{
  //The initMemory() function allocates memory
   int nMax = numLayers+1;
   int iMax = numIntervals+1;
   try {
      thetaStart = t1D(nMax);
      intensity = t1D(iMax);
      sqrDelta = t1D(iMax);
      growthTime = t1D(iMax);
      theta = t1D(nMax);
      dThetaDt = t1D(nMax);
      C = t1D(nMax);
      gR = t1D(nMax);
      coverage = t2D(nMax, growthTime);
   }
   catch (out_of_range o) {
      cout << o.what() << endl;
   }
}
//---------------------------------------------------------
double th(int n, t1D &theta){
   static const int nMax = theta.size()-1;

   if(n<=0) return 1.0;
      else if (n > nMax) return 0.0;
         else return theta[n];
}
//----------------------------------------------------------
double dn1(int n, t1D &theta) {
   double tt;
   tt = th(n, theta);
   if (tt < 0.0)
      tt = 0.0;
   if (tt > 1.0)
      tt = 1.0;
   return tt*sqrt(1.0-tt);
};
//----------------------------------------------------------
double dn2(int n, t1D &theta) {
   double tt;
   tt = th(n, theta);
   if (tt < 0.0)
      tt = 0.0;
   if (tt > 1.0)
      tt = 1.0;
   return tt*(1.0-tt);
};
//----------------------------------------------------------
double dn3(int n, t1D &theta) {
   double tt;
   tt = th(n, theta);
   if (tt < 0.0)
      tt = 0.0;
   if (tt > 1.0)
      tt = 1.0;
   if (tt < 0.5)
      return sqrt(tt);
      else
        return sqrt(1.0-tt);
};
//----------------------------------------------------------
//  The derivsX() are the user-supplied functions that compute the derivatives dTheta/dt
//----------------------------------------------------------
void derivs0(t1D C, t1D gR, t1D theta, t1D &dThetaDt){
//Diffusive growth
   double dTheta;
   static const int nMax = theta.size();

   for(int n=1; n<nMax; n++){
      dTheta = th(n-1, theta)-th(n, theta);
      dThetaDt[n] = dTheta*gR[n] + C[n]*(th(n+1, theta)-th(n+2, theta))*dTheta -
                    C[n]*(th(n, theta)-th(n+1, theta))*(th(n-2, theta)-th(n-1, theta));
   }
}
//---------------------------------------------------------
void derivs1(t1D C, t1D gR, t1D theta, t1D &dThetaDt){
//Distributed growth - variant 1
   double dTheta;
   static const int nMax = theta.size();

   for (int n=1; n<nMax; n++){
      dTheta = th(n-1, theta)-th(n, theta);
      if (dn1(n-1, theta) != 0.0)
        dTheta -= C[n-1]*dTheta*dn1(n-1, theta)/(dn1(n-1, theta)+dn1(n, theta));
      if (dn1(n, theta) > 0.0)
        dTheta += C[n]*(th(n, theta)-th(n+1, theta))*dn1(n, theta)/(dn1(n, theta)+dn1(n+1, theta));
      dThetaDt[n] = dTheta*gR[n];
   }
}
//---------------------------------------------------------
void derivs2(t1D C, t1D gR, t1D theta, t1D &dThetaDt){
//Distributed growth - variant 2
   double dTheta;
   static const int nMax = theta.size();

   for (int n=1; n<nMax; n++){
      dTheta = th(n-1, theta)-th(n, theta);
      if (dn2(n-1, theta) != 0.0)
        dTheta -= C[n-1]*dTheta*dn2(n-1, theta)/(dn2(n-1, theta)+dn2(n, theta));
      if (dn2(n, theta) > 0.0)
        dTheta += C[n]*(th(n, theta)-th(n+1, theta))*dn2(n, theta)/(dn1(n, theta)+dn2(n+1, theta));
      dThetaDt[n] = dTheta*gR[n];
   }
}
//---------------------------------------------------------
void derivs3(t1D C, t1D gR, t1D theta, t1D &dThetaDt){
//Distributed growth - variant 3
   double dTheta;
   static const int nMax = theta.size();

   for (int n=1; n<nMax; n++){
      dTheta = th(n-1, theta)-th(n, theta);
      if (dn3(n-1,theta) != 0.0)
        dTheta -= C[n-1]*dTheta*dn3(n-1, theta)/(dn3(n-1, theta)+dn3(n, theta));
      if (dn3(n, theta) > 0.0)
        dTheta += C[n]*(th(n, theta)-th(n+1, theta))*dn3(n, theta)/(dn3(n, theta)+dn3(n+1, theta));
      dThetaDt[n] = dTheta*gR[n] ;
   }
}
//---------------------------------------------------------
void rkdp(t1D C, t1D gR, t1D y, t1D dydx, double h, t1D &yout, t1D &yerr,
          void (*derivs)(t1D C, t1D gR, t1D theta, t1D &dThetaDt))
{
// Dormand-Prince fifth-order Runge-Kutta method. For more details See Refs. [14,15]

   static const double //c2=0.2,c3=0.3,c4=0.8,c5=8.0/9.0,
                       a21=0.2, a31=3.0/40.0, a32=9.0/40.0, a41=44.0/45.0,
                       a42=-56.0/15.0, a43=32.0/9.0, a51=19372.0/6561.0,
                       a52=-25360.0/2187.0, a53=64448.0/6561.0, a54=-212.0/729.0,
                       a61=9017.0/3168.0, a62=-355.0/33.0, a63=46732.0/5247.0,
                       a64=49.0/176.0, a65=-5103.0/18656.0,
                       a71=35.0/384.0, a73=500.0/1113.0, a74=125.0/192.0,
                       a75=-2187.0/6784.0, a76=11.0/84.0,e1=71.0/57600.0,
                       e3=-71.0/16695.0, e4=71.0/1920.0, e5=-17253.0/339200.0,
                       e6=22.0/525.0, e7=-1.0/40.0;
   int i;
   static const int nMax = y.size();
   t1D ytemp, dydxnew, k2, k3, k4, k5, k6;
   k2=t1D(nMax);
   k3=t1D(nMax);
   k4=t1D(nMax);
   k5=t1D(nMax);
   k6=t1D(nMax);
   ytemp=t1D(nMax);
   dydxnew=t1D(nMax);

   for (i=0;i<nMax;i++)     //First step
      ytemp[i]=y[i]+h*a21*dydx[i];
   derivs(C, gR, ytemp,k2); //Second step
   for (i=0;i<nMax;i++)
      ytemp[i]=y[i]+h*(a31*dydx[i]+a32*k2[i]);
   derivs(C, gR, ytemp,k3); //Third step
   for (i=0;i<nMax;i++)
      ytemp[i]=y[i]+h*(a41*dydx[i]+a42*k2[i]+a43*k3[i]);
   derivs(C, gR, ytemp,k4); //Fourth step
   for (i=0;i<nMax;i++)
      ytemp[i]=y[i]+h*(a51*dydx[i]+a52*k2[i]+a53*k3[i]+a54*k4[i]);
   derivs(C, gR, ytemp,k5); //Fifth step
   for (i=0;i<nMax;i++)
      ytemp[i]=y[i]+h*(a61*dydx[i]+a62*k2[i]+a63*k3[i]+a64*k4[i]+a65*k5[i]);
   derivs(C, gR, ytemp,k6); //Sixth step
   for (i=0;i<nMax;i++)     //Accumulate increments with proper weights
      yout[i]=y[i]+h*(a71*dydx[i]+a73*k3[i]+a74*k4[i]+a75*k5[i]+a76*k6[i]);
   derivs(C, gR, yout, dydxnew);
   for (i=0;i<nMax;i++)     //Estimate error
      yerr[i]=h*(e1*dydx[i]+e3*k3[i]+e4*k4[i]+e5*k5[i]+e6*k6[i]+e7*dydxnew[i]);
}
//---------------------------------------------------------
void rkqs(t1D C, t1D gR, t1D &y, t1D &dydx, double &x, double htry, double eps,
          t1D yscal, double &hnext,
          void (*derivs)(t1D C, t1D gR, t1D theta, t1D &dThetaDt))

//Fifth-order Runge Kutta integration step with monitoring of local truncation error to ensure accuracy and adjust stepsize.
//For more details See W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling,
//Numerical Recipes in C: The Art of Scientific Computing, second ed., Cambridge University Press, 1992
{
   static const int nMax = y.size();
   int i;
   double errmax,h,htemp;
   t1D yerr,ytemp;

   yerr=t1D(nMax);
   ytemp=t1D(nMax);

   h=htry; //initial stepsize
   for (;;) {
      rkdp(C, gR, y, dydx, h, ytemp, yerr, derivs);
      errmax=0.0; //compute scaled maximum error
      for (i=1;i<nMax;i++)
         errmax=max(errmax,fabs(yerr[i]/yscal[i]));
      errmax /= eps;
      if (errmax <= 1.0)
         break;
      htemp=SAFETY*h*pow(errmax,PSHRNK);
      h =(h >= 0.0 ? max(htemp,0.1*h) : min(htemp,0.1*h));
   }
   if (errmax > ERRCON) hnext=SAFETY*h*pow(errmax,PGROW);
      else
         hnext=5.0*h;
   x += h;
   for (i=1;i<nMax;i++)
      y[i]=ytemp[i];
}
//---------------------------------------------------------
void odeint(t1D C, t1D gR, t2D &coverage, t1D &thetaStart, t1D &growthTime,
            int numIntervals, int &numReturn, double dXsav, double tMax,
            void (*derivs)(t1D C, t1D gR, t1D theta, t1D &dThetaDt),
            double t0, double eps, double h1)

//Runge-Kutta driver with adaptive stepsize control. For more details See
//W.H. Press, B.P. Flannery, S.A. Teukolsky, W.T. Vetterling,
//Numerical Recipes in C: The Art of Scientific Computing, second ed., Cambridge University Press, 1992
{
   static const int nMax = thetaStart.size();

   double xsav, x, hnext, h;
   t1D yscal, y, dydx;

   yscal = t1D(nMax);
   y = t1D(nMax);
   dydx = t1D(nMax);

   x = t0;
   h = h1;
   numReturn = 0;
   for (int i=1; i<nMax; i++)
      y[i] = thetaStart[i];
   xsav = x-dXsav*2.0;
   for (int nstp=1; nstp<=MAXSTP; nstp++) {
      derivs(C, gR, y, dydx);
      for (int i=1; i<nMax; i++)
         yscal[i] = fabs(y[i]) + fabs(dydx[i]*h)+TINY;
      if ((fabs(x-xsav) > fabs(dXsav)) && (numReturn < numIntervals-1)){
         growthTime[++numReturn] = x;
         for (int i=1; i<nMax; i++)
            coverage[i][numReturn] = y[i];
         xsav = x;
      };
      if (((x+h-tMax)*(x+h-t0)) > 0.0)
         h = tMax-x;
      rkqs(C, gR, y, dydx, x, h, eps, yscal, hnext, derivs);
      if (((x-tMax)*(tMax-t0)) >= 0.0) {
         for (int i=1; i<nMax; i++)
            thetaStart[i] = y[i];
         growthTime[++numReturn] = x;
         for (int i=1; i<nMax; i++)
            coverage[i][numReturn] = y[i];
         break;
     };
     h = hnext;
   };
};
//---------------------------------------------------------
void setOutputData(string outputFileC, t2D coverage, t1D growthTime, t1D thetaStart,
                   t1D intensity, t1D sqrDelta, t1D gR, int numLayers, int numReturn)
//The setOutputData() function saves the results in the current directory
{
   ofstream fileC(Str(outputFileC), ios::out);
   fileC.precision(5);
   fileC.setf(ios::fixed);

   fileC << numLayers << endl;
   fileC << numReturn << endl;
   for (int n=1 ; n<=numLayers; n++) {
      for (int t=1 ; t<=numReturn; t++) {
         //Writes the total coverages versus growth time
         fileC << growthTime[t] <<"\t"<< coverage[n][t] << endl;
      }
   }
   fileC.close();
}
//---------------------------------------------------------
void readInputData(int numLayers, int model, ifstream &fileTh, t1D &C, t1D &gR) {
//The vector C stores values of the parameters An or kn, respectively
//depending on the context of their use (distributed or diffusive growth model).
//The vector gR stores values of the growth rate [1/tau monolayer per unit time]
  for(int n=1; n<=numLayers; n++) {
     fileTh >> C[n] >> gR[n];
     if(C[n] < 0.0) C[n] = 0.0;
     if(gR[n] <= 0.0) gR[n] = 1.0;
     if ((model != 0) && (C[n] >= 1.0)) C[n] = 0.9999;
  }
}
//---------------------------------------------------------



