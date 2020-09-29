#include <fstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <iomanip>
#include "headers/dtp.h"
#include "headers/tm.h"
#include "headers/uc.h"

using namespace std;

//-- Function prototypes----------------------------------
void loadParam(string, double &, double &, double &, double &, double &);
void loadCoverages(string, int &, int &, t1D &, t1D &, t1D &, t1D &, t1D &, t1D &, t2D &, t2D &,
                   t2D &, t3D &, t3D &, t3D &, t3D &);
void thicknessZi_Substrate(t1D &, double &, double &);
void thicknessZn_Substrate(t1D &, double &);
void thicknessZi_GrowingLayers(t1D &, double &, double);
void thicknessZn_GrowingLayers(t1D &, double);
double sqrU(double, double, double);
void crystPotUgSubstrate(double, double, t1D, t1D, t2D &);
void crystPotUgGrowingLayers(double, double, t1D, t1D, t2D &, t2D &);
void matrix_Mn(int, double, double, double, const double, int, t1D, t1D &, t1D, t2D,
               t2D, t2D, t3D &, t3D &, t3D &, t3D &);
void matrix_M(int, t1D, t3D, t3D, t3D, t3D, t2D &, t2D &);
void amplitudeR0(const double, int, ofstream &, t1D, t2D, t2D);
void saveCrystallinePotentials(string, t1D, t1D);
void intensNorm(int, string, string, t1D, t1D);
//---------------------------------------------------------
t1D  a(4), aa(4), b(4), bb(4);
string coverages = "coverage.dat";         //input file with growth model
string RHEEDparam = "inputRheedData.dat";  //input file with RHEED parameters

int main()
{
   int numReturn;
   int numberOfGrowingLayers;
   double ti, s, h1, msda, alpha, beta, T, E, iAngle;
   t1D Z, Zn, Ugz, UgzAdd, gt, R0Norm;
   t2D Ug, UgAdd, coverage, ReR(2, t1D(2)), ImR(2, t1D(2));
   t3D ReM, ImM, ReImM, ImReM;

   loadParam(RHEEDparam, alpha, beta, T, E, iAngle);
   loadCoverages(coverages, numReturn, numberOfGrowingLayers, gt, Z, Zn, Ugz,
                 UgzAdd, R0Norm, coverage, Ug, UgAdd, ReM, ImM, ReImM, ImReM);

   const double Kappa = sqrt(2*m0c2*E*(1+E/(2*m0c2)))/(c*hk);
   const double KappaZ = Kappa*sin(M_PI*iAngle/180.0);  // electron wave vector (z-direction) in the vacuum
   const double sKz = sqr(KappaZ);

   thicknessZi_Substrate(Z, ti, s);
   thicknessZi_GrowingLayers(Z, ti, s);
   thicknessZn_Substrate(Zn, h1);
   msda = sqrU(T, mAtSub, TDSub);
   crystPotUgSubstrate(E, msda, Z, Zn, Ug);

   thicknessZn_GrowingLayers(Zn, h1);
   msda = sqrU(T, mAtGl, TDGl);
   crystPotUgGrowingLayers(E, msda, Z, Zn, Ug, UgAdd);

   string fName = coverages.erase(coverages.find(".dat")) + "-L" + Str(numberOfGrowingLayers) + \
                  " iA=" + Str(iAngle) +" E=" +Str(E) +" T="+ Str(T) +".dat";

   ofstream fileR0(Str(fName), ios::out);
   fileR0.precision(5);
   fileR0.setf(ios::fixed);

   for(int t = 0; t < numReturn; t++){
      matrix_Mn(numberOfGrowingLayers, alpha, beta, ti, sKz, t, Z, Ugz, UgzAdd, coverage, Ug,
                UgAdd, ReM, ImM, ReImM, ImReM);
      for(int jj = 0; jj < 2; jj++) {
         matrix_M(jj, Z, ReM, ImM, ReImM, ImReM, ReR, ImR);
      }
      amplitudeR0(KappaZ, t, fileR0, gt, ReR, ImR);
   }

   fileR0.close();

   saveCrystallinePotentials("potential.dat", Z, Ugz);

   string nfName = "norm-"+fName;
   intensNorm(numReturn, Str(fName), Str(nfName), R0Norm, gt);
   cout << "End of calculations...";
   //cin.get();
   return EXIT_SUCCESS;
}
//-- Bodies of functions-----------------------------------
double sqrU(double Tx, double mAtx, double TDx)
{
  //The sqrU() function computes the mean square displacement of atoms (msda) along the z-axis
  //according to the Debye's model

  double fi, cD, sum = 0.0;
  cD = 1e+20*3*sqr(h)/(mAtx*kb*TDx); //[square Angstrom]
  const double Th = Tx/TDx;
  for(int j = 1; j <= 3; j++) {
    sum += (j*Th+sqr(Th))*exp(-j*1.0/Th)/sqr(j);
  };
  if(Tx >= TDx)
    fi = Th+1/(36*Th)-pow(Th,-3.0)/3600.0;
    else
      fi = 0.25+sqr(M_PI*Th)/6.0-sum;
  return cD*(0.25+Th*fi);
}
//---------------------------------------------------------
void parametersTv(int k, double msda)
{
  //To interpret actual experimental data, it is important to take account the effect
  //of lattice thermal vibrations by replacing the Doyle and Turner coefficients bk
  //by bk+8*M_MP^2*u^2, where u^2 denotes the mean square displacement of atoms (msda) along the z-axis.
  //For more details See P.Mazurek, A.Daniluk, K.Paprocki, Vacuum 72 (2004) 363

  b[k] += 8*sqr(M_PI)*msda;
}
//---------------------------------------------------------
void  thicknessZi_Substrate(t1D &Z, double &ti, double &s)
{
  //The thicknessZi_Substrate() function computes the thin slice
  //position along the axis perpendicular to the substrate surface

  s = -d0Sub;
  ti = d0Sub/nsl;
  for(int i = 0; i < nSubstrate; i++) {
    Z[i] = s;
    s += ti;
  }
}
//---------------------------------------------------------
void  thicknessZi_GrowingLayers(t1D &Z, double &ti, double s)
{
  //The thicknessZi_GrowingLayers() function computes the thin slice
  //position along the axis perpendicular to the growing layers surface

  ti = d0Gl/nsl;
  for(int i = nSubstrate; i < Z.size(); i++) {
    Z[i] = s;
    s += ti;
  }
}
//---------------------------------------------------------
void  thicknessZn_Substrate(t1D &Zn, double &h1)
{
  //The thicknessZn_Substrate() function computes the layer position
  //along the axis perpendicular to the substrate surface
  h1=0.0;
  for (int n = 0; n < numberOfSubstrateLayers; n++){
    Zn[n] = h1;
    h1 += d0Sub;
  }
}
//---------------------------------------------------------
void  thicknessZn_GrowingLayers(t1D &Zn, double h1)
{
  //The thicknessZn_GrowingLayers() function computes the layer position
  //along the axis perpendicular to the growing layers surface

  h1 -= d0Sub;
  h1 += d0Gl;
  for (int n = numberOfSubstrateLayers; n < Zn.size(); n++){
    Zn[n] = h1;
    h1 += d0Gl;
  }
}
//---------------------------------------------------------
void  crystPotUgSubstrate(double E, double msda, t1D Z, t1D Zn, t2D &Ug)
{
  //The  crystPotUgSubstrate() function computes the zeroth Fourier
  //component of the scattering potential of the substrate

  double UU,U1,U2,U3,result;

  a = {sa0,sa1,sa2,sa3};

  for(int n = 0; n < numberOfSubstrateLayers; n++) {
     for(int i = 0; i < Z.size(); i++) {
        //crystalline Si has similar fcc lattice with two-atoms basis
        Ug[n][i] = -(1+E/m0c2)*8*M_PI/S0Sub;
        UU = 0.0f;
        for(int k = 0; k < 4; k++) {
           b = {sb0,sb1,sb2,sb3};
           parametersTv(k, msda);
           U1 = sqrt(M_PI/b[k])*a[k];
           result = (-4)*sqr(M_PI)*sqr(Z[i]-Zn[n])/b[k];
           U2 = exp(result);
           result = (-4)*sqr(M_PI)*sqr(Z[i]-Zn[n]-(d0Sub/4.0f))/b[k];
           U3 = exp(result);
           UU += U1*(U2+U3);
        }
        Ug[n][i] *= UU;
     }
  }
}
//---------------------------------------------------------
void  crystPotUgGrowingLayers(double E, double msda, t1D Z, t1D Zn, t2D &Ug, t2D &UgAdd)
{
  //The crystPotUgLayers() function computes the zeroth Fourier
  //component of the scattering potential of the growing layers

  double UU,U1,U2,U3,Uadd,result,nFactor;

  a = {ga0,ga1,ga2,ga3};

  nFactor =  -(1+E/m0c2)*8*M_PI/S0Gl;  //normalization factor for the growing layers

  for(int n = numberOfSubstrateLayers; n < Zn.size(); n++) {
     for(int i = 0; i < Z.size(); i++) {
       //crystalline Ag has fcc lattice with one atom basis
        Ug[n][i] = nFactor;
        U3 = (-4)*sqr(M_PI)*sqr(Z[i]-Zn[n])/100; //For Empirical Model of the scattering potential.
                                                 //For more details See: Z. Mitura, S.L. Dudarev, L.M. Peng, G. Gladyszewski, M.J. Whelan, J. Cryst. Growth 235 (2002) 79.
        UU = 0.0;
        Uadd = 0.0;
        for(int k = 0; k < 4; k++) {
           b = {gb0,gb1,gb2,gb3};
           parametersTv(k, msda);
           U1 = sqrt(M_PI/b[k])*a[k];
           result = (-4)*sqr(M_PI)*sqr(Z[i]-Zn[n])/b[k];
           U2 = exp(result);
           Uadd += (a[k]/sqrt(100.0))*exp(U3); //For Empirical Model of the scattering potential.
           UU += U1*U2;
        }
        Ug[n][i] *= UU;
        UgAdd[n][i] = nFactor*Uadd;
     }
  }
}
//---------------------------------------------------------
void  matrix_Mn(int numberOfGrowingLayers, double alpha, double beta, double ti,
                const double sKz, int t, t1D Z, t1D &Ugz, t1D UgzAdd, t2D coverage, t2D Ug, t2D UgAdd,
                t3D &ReM, t3D &ImM, t3D &ReImM, t3D &ImReM)
{
  //The  matrix_Mn() function computes the scattering matrix for relating electron wave function and its
  //surface normal derivative on the upper and lower surfaces of an assembly of thin slices,
  //each having thickness ti

  double UUU, UUUadd, P1, P2, P3, X, Y, A, B ,Xs, Ys, De, lcr;

  for(int i = 0; i < Z.size(); i++) {
    UUU = 0.0;
    UUUadd = 0.0;
    for(int n = 0; n < numberOfSubstrateLayers+numberOfGrowingLayers; n++)  {
       if(n < numberOfSubstrateLayers)
         lcr = 1.0;
         else
           lcr = coverage[n-numberOfSubstrateLayers][t];

       UUU += lcr*Ug[n][i];    //Proportional model of the scattering potential
       UUUadd += lcr*UgAdd[n][i];
    }
    Ugz[i] = UUU;
    UgzAdd[i] = UUUadd;

    //A - real part of the potential
    A = sKz -Ugz[i];
    B = -alpha*Ugz[i] - beta*UgzAdd[i];
    //B - imaginary part of the potential
    //For beta = 0.0 the imaginary part of the potential reduces to the proportional model

    X = sqrt((A+sqrt(sqr(A)+sqr(B)))/2.0);
    Y = B/(2*X);
    Ys = Y*ti;
    Xs = X*ti;
    De = 2*(sqr(X)+sqr(Y));
    P1 = cos(Xs);
    P2 = sin(Xs);
    P3 = -P2;
    ReM[0][0][i] = (exp(-Ys)*P1+exp(Ys)*P1 )/2.0;
    ImM[0][0][i] = (exp(-Ys)*P2+exp(Ys)*P3)/2.0;
    ReM[0][1][i] = (exp(-Ys)*(-Y*P1+X*P2)+exp(Ys)*(Y*P1-X*P3))/De;
    ImM[0][1][i] = (exp(-Ys)*(-X*P1-Y*P2)+exp(Ys)*(X*P1+Y*P3))/De;
    ReM[1][0][i] = (exp(-Ys)*(-X*P2-Y*P1)+exp(Ys)*(X*P3+Y*P1))/2.0;
    ImM[1][0][i] = (exp(-Ys)*(X*P1-Y*P2)+exp(Ys)*(-X*P1+Y*P3))/2.0;
    ReM[1][1][i] = ReM[0][0][i];
    ImM[1][1][i] = ImM[0][0][i];
  }
}
//---------------------------------------------------------
void  matrix_M(int jj, t1D Z, t3D ReM, t3D ImM, t3D ReImM, t3D ImReM, t2D &ReR, t2D &ImR)
{
  //The matrix_M() function computes the complex elements of the matrix Mn

  int j;
  double SumR, SumI, SumRI, SumIR;

  for(int i = 0; i < Z.size()-1; i++) {
     j = jj;
     for(int r = 0; r < 2; r++) {
        SumR = 0.0;
        SumI = 0.0;
        SumRI = 0.0;
        SumIR = 0.0;
        for(int l = 0; l < 2; l++) {
           SumR += ReM[j][l][i]*ReM[l][r][i+1];
           SumI += ImM[j][l][i]*ImM[l][r][i+1];
           SumRI += ReM[j][l][i]*ImM[l][r][i+1];
           SumIR += ImM[j][l][i]*ReM[l][r][i+1];
        }
        ReM[j][r][i+1] = SumR;
        ImM[j][r][i+1] = SumI;
        ReM[j][r][i+1] = ReM[j][r][i+1]-ImM[j][r][i+1];
        ReImM[j][r][i+1] = SumRI;
        ImReM[j][r][i+1] = SumIR;
        ImM[j][r][i+1] = ReImM[j][r][i+1]+ImReM[j][r][i+1];
        if( i == (Z.size()-2)) {
          ReR[j][r] = ReM[j][r][Z.size()-1];
          ImR[j][r] = ImM[j][r][Z.size()-1];
        }
     }
  }
}
//---------------------------------------------------------
void  amplitudeR0(const double KappaZ, int t, ofstream &fileR0, t1D gt, t2D ReR, t2D ImR)
{
  //The amplitudeR0() function computes the specular reflected beam amplitude R0

  double D1, D2, Y1, Y2, W, X, Y, R0;

  D1 = ReR[1][1]-ReR[0][0]+ImR[1][0]/KappaZ+KappaZ*ImR[0][1];
  Y1 = -ReR[1][0]/KappaZ-KappaZ*ReR[0][1]-ImR[0][0]+ImR[1][1];
  D2 = ReR[0][0]+ReR[1][1]-ImR[1][0]/KappaZ+KappaZ*ImR[0][1];
  Y2 = ReR[1][0]/KappaZ-KappaZ*ReR[0][1]+ImR[0][0]+ImR[1][1];
  X = D1*D2+Y1*Y2;
  Y = D2*Y1-D1*Y2;
  W = sqr(D2)+sqr(Y2);
  X = X/W;
  Y = Y/W;
  R0 = sqr(X)+sqr(Y);

  cout << gt[t] << '\t' << R0 << endl;
  fileR0 << gt[t] << '\t' << R0 << endl;
}
//---------------------------------------------------------
void  initMemory(int numReturn, int numberOfGrowingLayers, t1D &gt, t1D &Z, t1D &Zn, t1D &Ugz, t1D &UgzAdd, t1D &R0Norm,
                 t2D &coverage, t2D &Ug, t2D &UgAdd, t3D &ReM, t3D &ImM, t3D &ReImM, t3D &ImReM)
{
  //The initMemory() function allocates memory

  const int nLayers = (numberOfGrowingLayers+1.5)*nsl; //number of thin slices in growing layer
  try {
    gt=t1D(numReturn);
    R0Norm=t1D(numReturn);
    Z=t1D(nLayers+nSubstrate);
    Ugz=t1D(nLayers+nSubstrate);
    UgzAdd=t1D(nLayers+nSubstrate);
    Zn=t1D(numberOfSubstrateLayers+numberOfGrowingLayers);

    coverage = t2D(numberOfGrowingLayers, gt);
    Ug=t2D(numberOfSubstrateLayers+numberOfGrowingLayers, Z);
    UgAdd=t2D(numberOfSubstrateLayers+numberOfGrowingLayers, Z);

    ReM=t3D(2, t2D(2, Z));
    ImM=t3D(2, t2D(2, Z));
    ReImM=t3D(2, t2D(2, Z));
    ImReM=t3D(2, t2D(2, Z));
  }
  catch (out_of_range o) {
      cout << o.what() << endl;
  }
}
//---------------------------------------------------------
void loadParam(string RHEEDparam, double &alpha, double &beta, double &T, double &E, double &iAngle) {

   ifstream fileParam(Str(RHEEDparam), ios::in);
   if(!fileParam) {
      cout << "Input RHEED parameters file: " <<RHEEDparam<<" not found";
      cin.get();
      exit(EXIT_FAILURE);
   }
   fileParam >> alpha;                 //imaginary part of the potential
   fileParam >> beta;                  //imaginary part of the potential
   fileParam >> T;                     //temperature of the crystal
   fileParam >> E;                     //energy electron beam
   fileParam >> iAngle;                //glancing angle of the incident beam
   fileParam.close();
}
//---------------------------------------------------------
void loadCoverages(string coverages, int &numReturn, int &numberOfGrowingLayers, t1D &gt, t1D &Z, t1D &Zn, t1D &Ugz,
                   t1D &UgzAdd, t1D &R0Norm, t2D &coverage, t2D &Ug, t2D &UgAdd, t3D &ReM, t3D &ImM, t3D &ReImM, t3D &ImReM)
{
   ifstream fileTh(Str(coverages), ios::in);
   if(!fileTh) {
      cout << "Input file: " <<coverages<<" not found";
      cin.get();
      exit(EXIT_FAILURE);
   }
   fileTh >> numberOfGrowingLayers;
   fileTh >> numReturn;

   initMemory(numReturn, numberOfGrowingLayers, gt, Z, Zn, Ugz, UgzAdd, R0Norm, coverage, Ug, UgAdd, ReM, ImM, ReImM, ImReM);
   for(int n = 0; n < numberOfGrowingLayers; n++)
      for(int t = 0; t < numReturn; t++)
        fileTh >> gt[t] >> coverage[n][t];
   fileTh.close();
}
//---------------------------------------------------------
void saveCrystallinePotentials(string fName, t1D Z, t1D Ugz){

  //The saveCrystallinePotentials() function saves the crystalline potentials (real part)

  ofstream fileUz(Str(fName), ios::out);
  fileUz.precision(5);
  fileUz.setf(ios::fixed);
  for(int i = 0; i < Z.size(); i++)
     fileUz << Z[i] << "\t" << Ugz[i]*sqr(hk)*1.6021e-19/(2*m0*1e-20)<< endl;
     // 1eV = 1.6021e-19 J
     // square Angstrom = 1e-20 square metre
  fileUz.close();
}
//---------------------------------------------------------
void intensNorm(int numReturn, string fName, string nfName, t1D R0Norm, t1D gt){

   //The IntensNorm() function saves normalized values of the specular reflected beam amplitude R0

   ifstream fileIn(Str(fName), ios::in);
   for(int i = 0; i < numReturn; i++)
      fileIn >> gt[i] >> R0Norm[i];
   fileIn.close();

   double mn = *min_element(R0Norm.begin(),R0Norm.end());
   double mx = *max_element(R0Norm.begin(),R0Norm.end());

   ofstream fileN(Str(nfName), ios::out);
   fileN.precision(5);
   fileN.setf(ios::fixed);
   for(int i = 0; i < numReturn; i++){
      R0Norm[i]=(R0Norm[i]-mn)/(mx-mn);
      fileN << gt[i]<< "\t" << R0Norm[i] << endl;
   }
   fileN.close();
}
//----------------------------------------------------------
