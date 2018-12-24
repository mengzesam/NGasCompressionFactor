#ifndef AGA892DC_H
#define AGA892DC_H

#include <math.h>
class AGA892DC{
    //ID:   1	    2	3	4	    5	    6	7	8	9	10	11	    12	    13	    14	    15	    16	    17	1   8	    19	    20	21
    //Comp.:CH4	    N2	CO2	C2H6	C3H8	H2O	H2S	H2	CO	O2	i-C4H10	n-C4H10	i-C5H12	n-C5H12	C6H14	C7H16	C8H18	C9H20	C10H22	He	Ar
public://member method
   static int calcZRho(double xi_raw[],int ID_raw[],int NUM_x,double p,double t,
                        double& Z,double& rhom,double& rho);
private://member method
    static int checkInputs(double xi_raw[],int ID_raw[],int NUM_x,double p,double t);
    static void sortInput(double xi_raw[],int ID_raw[],int ID[],int Index2raw[],int NUM_x);
    static void calcSecondVirial(double xi_raw[],int ID[],int Index2raw[],int NUM_x,double BI[18]);
    static void calcKFQGU(double xi_raw[],int ID[],int Index2raw[],int NUM_x,
                          double& K,double& F,double& Q,double& G,double& U);
    static double rhom2PZ(double rhom,double tt,double BMIX,double CNS[46],double K,
                        double F,double Q,double G,double U,double& Z);
private:
    static int status;
private://state parameters ref. to B.1 of GB/T 17747.2-2011 (or ISO 12213-2 2006)
    const static int bn[58];
    const static int cn[58];
    const static int kn[58];
    const static int gn[58];
    const static int qn[58];
    const static int fn[58];
    const static int sn[58];
    const static int wn[58];
    const static double an[58];
    const static double un[58];
private://characterization parameters ref. to B.1 of GB/T 17747.2-2011 (or ISO 12213-2 2006)
    const static double Mi[21];//kg/kmol
    const static double Ei[21];//K
    const static double Ki[21];//
    const static double Gi[21];//(m3/kmol)^(1/3)
    const static double Qi3,Qi6,Qi7;
    const static double Si6,Si7;
private://binary interaction parameters
    const static double E1j[18];//j:2-19
    const static double U1j[18];//j:2-19
    const static double K1j[18];//j:2-19
    const static double G13;
    const static double G18;
    const static double E2j[12];//j:3-14
    const static double U2j[10];//j:3-12
    const static double K2j[6];//j:3-8
    const static double G23;
    const static double E3j[16];//j:4-19
    const static double U3j[16];//j:4-19
    const static double K3j[16];//j:4-19
    const static double G3j[3];//j:4-6
    const static double E4j[10];//j:5-14
    const static double U4j[10];//j:5-14
    const static double K4j[4];//j:5-8
    const static double E58,E512;
    const static double E7j[5];//j:15-19
    const static double U7j[5];//j:15-19
    const static double K7j[5];//j:15-19
    const static double E89,E811,E812;
};
#endif