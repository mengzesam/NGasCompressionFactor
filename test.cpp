#undef BUILDING_DLL
#include "NGCompressionFactor.h"
#include <iostream>
#include <iomanip>
using namespace std;

/* int main(int argc, char const *argv[])
{
    double rhom,rho,Z;
    double p=6.0;
    double t;
    int NUM_x=14;
    int ID_raw[14]={
        3,2,8,9,1,4,5,11,12,13,14,15,16,17
    };
    double xi_raw[14]={
       0.011,0.117,0,0,0.826,0.035,0.0075,0.0012,0.0012,0.0004,0.0004,0.0002,0.0001,0
    };
    t=-3.15;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    t=6.85;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    t=16.85;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    t=36.85;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    t=56.85;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    p=12.0;
    t=-3.15;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    t=6.85;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    t=16.85;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    t=36.85;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    t=56.85;
    Z=xipt2Z_AGA(xi_raw,ID_raw,NUM_x,p,t);
    cout<<setprecision(10)<<Z<<endl;
    return 0;
} */



int main(){
    double xCO2,xH2,d,Hs,p,t;
    double Z,rho,rhom;
    xCO2=0.006;
    xH2=0.0;
    d=0.581;
    Hs=40.66;
    p=6;
    t=-3.15;
    Z=Z_SGERG88(p,t,Hs,d,xCO2,xH2);
    return 0;
}