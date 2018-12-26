#ifndef NGCOMPRESSIONFACTOR_H
#define NGCOMPRESSIONFACTOR_H

#include "AGA892DC.h"
#include "SGERG88.h"

#ifdef BUILDING_DLL
    #define DLL_PUBLIC extern "C" __declspec(dllexport)  __stdcall       
#else
    #define DLL_PUBLIC extern "C" __declspec(dllimport)  __stdcall   
#endif
//__attribute__ ((dllexport)) __stdcall
//__attribute__ ((dllimport)) __stdcall

DLL_PUBLIC  double xipt2Z_AGA(double xi_raw[],int ID_raw[],int NUM_x,double p,double t);
DLL_PUBLIC  double xipt2rhom_AGA(double xi_raw[],int ID_raw[],int NUM_x,double p,double t);
DLL_PUBLIC  double xipt2rho_AGA(double xi_raw[],int ID_raw[],int NUM_x,double p,double t);
DLL_PUBLIC  double Z_SGERG88(double p,double t,double Hs,double d,double xCO2,double xH2);
DLL_PUBLIC  double Rhom_SGERG88(double p,double t,double Hs,double d,double xCO2,double xH2);
DLL_PUBLIC  double Rho_SGERG88(double p,double t,double Hs,double d,double xCO2,double xH2);



#endif