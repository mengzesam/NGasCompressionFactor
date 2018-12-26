#define BUILDING_DLL

#include "NGCompressionFactor.h"
#include <windows.h>

/*
    AGA892DC:
    //xi:摩尔分数,p:绝对压力,MPa,t:摄氏度，C
    /*
    p:0<p<=12MPa
    t+273.15:263-338
    x(CH4):0.7-1.0
    x(N2):0-0.7
    x(CO2):0-0.2
    x(C2H6):0-0.1
    x(C3H8):0-0.035
    x(C4H10):0-0.015
    x(C5H12):0-0.005
    x(C6H14):0-0.001
    x(C7H16):0-0.0005
    x(C8+):0-0.0005
    x(H2):0-0.1
    x(CO):0-0.03
    x(He):0-0.005
    x(H2O):0-0.00015
    x(O2):0-0.00015?
    x(H2S):0-0.00015?
    x(Ar+Ne+Kr+Xe):0-0.00015?
    sum(xi)==1.0
    ID:         1	     2	3	4	    5	    6	7	8	9	10	
    Component:  CH4	    N2	CO2	C2H6	C3H8	H2O	H2S	H2	CO	O2	
    ID:             11	    12	    13	    14	    15	    16	    17	1   8	    19	    20	21
    Component:    i-C4H10	n-C4H10	i-C5H12	n-C5H12	C6H14	C7H16	C8H18	C9H20	C10H22	He	Ar
    
*/


DLL_PUBLIC  double xipt2Z_AGA(double xi_raw[],int ID_raw[],int NUM_x,double p,double t){
    double Z,rhom,rho;
    int ret=AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    if(ret==100) 
        return Z;
    else 
        return -999999.0;
}
DLL_PUBLIC  double xipt2rhom_AGA(double xi_raw[],int ID_raw[],int NUM_x,double p,double t){
    double Z,rhom,rho;
    int ret=AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    if(ret==100) 
        return rhom;
    else 
        return -999999.0;
}
DLL_PUBLIC  double xipt2rho_AGA(double xi_raw[],int ID_raw[],int NUM_x,double p,double t){
    double Z,rhom,rho;
    int ret=AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    if(ret==100) 
        return rho;
    else 
        return -999999.0;
}

DLL_PUBLIC  double Z_SGERG88(double p,double t,double Hs,double d,double xCO2,double xH2){
    double Z,rhom,rho;
    int ret=SGERG88::ZRho_SGERG88(p,t,Hs,d,xCO2,xH2,Z,rhom,rho);
    if(ret==0)
        return Z;
    else 
        return -999999.0;
}
DLL_PUBLIC  double Rhom_SGERG88(double p,double t,double Hs,double d,double xCO2,double xH2){
    double Z,rhom,rho;
    int ret=SGERG88::ZRho_SGERG88(p,t,Hs,d,xCO2,xH2,Z,rhom,rho);
    if(ret==0)
        return rhom;
    else 
        return -999999.0;
}
DLL_PUBLIC  double Rho_SGERG88(double p,double t,double Hs,double d,double xCO2,double xH2){
    double Z,rhom,rho;
    int ret=SGERG88::ZRho_SGERG88(p,t,Hs,d,xCO2,xH2,Z,rhom,rho);
    if(ret==0)
        return rho;
    else 
        return -999999.0;
}



BOOL APIENTRY DllMain( HANDLE hModule, DWORD fdwReason, LPVOID lpReserved)
{
    switch (fdwReason)
    {
    case DLL_PROCESS_ATTACH:
        break;
    case DLL_PROCESS_DETACH:
        break;
    case DLL_THREAD_ATTACH:
        break;
    case DLL_THREAD_DETACH: 
        break;
    }
    return TRUE;
}


