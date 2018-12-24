#ifndef NGCOMPRESSIONFACTOR_H
#define NGCOMPRESSIONFACTOR_H

#include "AGA892DC.h"

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




#endif