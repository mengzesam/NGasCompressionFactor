#ifndef SGERG88_H
#define SGERG88_H

#include <math.h>

class SGERG88{
public://public function
    static int ZRho_SGERG88(double p,double t,double Hs,double d,double xCO2,double xH2,
                      double& Z,double& rhom,double& rho);
private://member function
    static int checkInputs(double p,double t,double Hs,double d,double xCO2,double xH2);
    static int calcIntermediate(double Hs,double d,double xCO2,double xH2,
                          double& xCH,double& xN2,double& xCO,double& HCH,double& MCH,
                          double& Bn,double& rhomn,double& rhon);//d:relative density (comparing to air?)
private://member data:constants for B.1-B.2 Equation
    const static double HH2;
    const static double HCO;
    const static double MN2;
    const static double MCO2;
    const static double MH2;
    const static double MCO;
    const static double R;
    const static double Vmnid;
    const static double Rhonair;
private://member data:constants of numerical values for coefficients b(0)-b(2) 
        //in the temperature expansion of the second virial coefficients for 
        //pure gases and of the unlike-interaction virial coefficients
    const static double b_H0[3];
    const static double b_H1[3];
    const static double b_H2[3];
    const static double b_22[3];
    const static double b_33[3];
    const static double b_44[3];
    const static double b_55[3];
    const static double b_14[3];
    const static double b_15[3];
    const static double b_23[3];
    const static double b_24[3];;
private://member data:constants of numerical values of the coefficients c(0)-c(2) in the temperture 
        //expansion of the third virial coefficient for pure gases and of the unlinke-interaction
        //virial coefficients
    const static double c_H0[3];
    const static double c_H1[3];
    const static double c_H2[3];
    const static double c_222[3];
    const static double c_333[3];
    const static double c_444[3];
    const static double c_115[3];
    const static double c_223[3];
    const static double c_233[3];
};


#endif