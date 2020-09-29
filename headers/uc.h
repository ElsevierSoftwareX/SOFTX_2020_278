#ifndef __UC_H
#define __UC_H

#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265
#endif


//== Universal physical constants ==========================
const double c    = 2.997e+18;   //speed of light [Angstrom/s]
const double hk   = 6.58e-16;    //Planck's constant [eV*s]
const double h    = 1.0545e-34;  //Planck's constant [J*s]
const double m0   = 9.1091e-31;  //rest electron mass [kg]
const double m0c2 = 510.9989e+3; //rest electron mass [eV]
const double kb   = 1.38054e-23; //Boltzmann's constant [J/K]


//=== Constants for the crystals under investigation =======

//---Substrate----------------------------------------------
const double mAtSub  = 28.086*1.656e-27; //atomic mass of Si [kg]
const double TDSub   = 625.0;            //Debye temperature of Si [K]
const double A0Sub = 5.43072;            //lattice constant of Si [Angstrom]
const double npsSub  = 3.84;             //net point spacing in plane (111) of Si [Angstrom]
const double d0Sub = A0Sub*sqrt(3)/3.0;  //distance between atomic planes for Si(111)
const double S0Sub = sqr(npsSub)*sqrt(3)/2.0; //area of the two-dimensional unit cell of Si parallel to the surface (111)

//---Growing Layers-----------------------------------------
const double mAtGl  = 107.87*1.656e-27; //atomic mass of Ag [kg]
const double TDGl   = 215.0;            //Debye temperature of Ag [K]
const double A0Gl = 4.09;               //lattice constant of Ag [Angstrom]
const double npsGl  = 2.892;            //net point spacing in plane (111) of Ag [Angstrom]

const double d0Gl = A0Gl*sqrt(3.0)/3.0; //distance between atomic planes for Ag(111)
const double S0Gl = sqr(npsGl)*sqrt(3.0)/2.0; //area of the two-dimensional unit cell of Ag parallel to the surface (111)



#endif // UC_H
