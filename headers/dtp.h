#ifndef __DTP_H
#define __DTP_H


const int numberOfSubstrateLayers = 10;                     //number of substrate layers
const int nsl = 100;                                        //number of thin slices in one layer
const int nSubstrate = ((numberOfSubstrateLayers+0)*nsl);   //number of thin slices in substrate


//Doyle and Turner ak, bk parameters of the analytic representation of the electron scattering factors
//P.A. Doyle, P.S. Turner, Acta Crystallogr. A 24 (1968) 390.

//Substrate Si
const double sa0 =2.1293;
const double sa1 =2.5333;
const double sa2 =0.8349;
const double sa3 =0.3216;

const double sb0 =57.7748;
const double sb1 =16.4756;
const double sb2 =2.8796;
const double sb3 =0.3386;

//Growing layers Ag
const double ga0 =2.0355;
const double ga1 =3.2716;
const double ga2 =2.5105;
const double ga3 =0.8372;

const double gb0 =61.4970;
const double gb1 =11.8237;
const double gb2 =2.8456;
const double gb3 =0.3271;



#endif // DTP_H
