#include "SGERG88.h"

//member data:constants for B.1-B.2 Equation
    const double SGERG88::HH2=285.83;
    const double SGERG88::HCO=282.98;
    const double SGERG88::MN2=28.0135;
    const double SGERG88::MCO2=44.01;
    const double SGERG88::MH2=2.0159;
    const double SGERG88::MCO=28.01;
    const double SGERG88::R=0.00831451;
    const double SGERG88::Vmnid=22.414097;
    const double SGERG88::Rhonair=1.292923;
//member data:the second virial coefficients b(0)-b(1)
    const double  SGERG88::b_H0[3]={
        -4.25468E-1,2.86500E-3,-4.62073E-6
    };
    const double  SGERG88::b_H1[3]={
        8.77118E-4,-5.56281E-6,8.81510E-9
    };
    const double  SGERG88::b_H2[3]={
        -8.24747E-7,4.31436E-9,-6.08319E-12
    };
    const double  SGERG88::b_22[3]={
        -1.44600E-1,7.40910E-4,-9.11950E-7
    };
    const double  SGERG88::b_33[3]={
        -8.68340E-1,4.03760E-3,-5.16570E-6
    };
    const double  SGERG88::b_44[3]={
        -1.10596E-3,8.13385E-5,-9.87220E-8
    };
    const double  SGERG88::b_55[3]={
        -1.30820E-1,6.02540E-4,-6.44300E-7
    };
    const double  SGERG88::b_14[3]={
        -5.21280E-2,2.71570E-4,-2.50000E-7
    };
    const double  SGERG88::b_15[3]={
        -6.87290E-2,-2.39381E-6,5.18195E-7
    };
    const double  SGERG88::b_23[3]={
        -3.39693E-1,1.61176E-3,-2.04429E-6
    };
    const double  SGERG88::b_24[3]={
        1.20000E-2,0,0
    };
//member data:the third virial coefficients c(0)-c(1)
    const double  SGERG88::c_H0[3]={
        -3.02488E-1,1.95861E-3,-3.16302E-6
    };
    const double  SGERG88::c_H1[3]={
        6.46422E-4,-4.22876E-6,6.88157E-9
    };
    const double  SGERG88::c_H2[3]={
        -3.32805E-7,2.23160E-9,-3.67713E-12
    };
    const double  SGERG88::c_222[3]={
        7.84980E-3,-3.98950E-5,6.11870E-8
    };
    const double  SGERG88::c_333[3]={
        2.05130E-3,3.48880E-5,-8.37030E-8
    };
    const double  SGERG88::c_444[3]={
        1.04711E-3,-3.64887E-6,4.67095E-9
    };
    const double  SGERG88::c_115[3]={
        7.36748E-3,-2.76578E-5,3.43051E-8
    };
    const double  SGERG88::c_223[3]={
        5.52066E-3,-1.68609E-5,1.57169E-8
    };
    const double  SGERG88::c_233[3]={
        3.58783E-3,8.06674E-6,-3.25798E-8
    };
//end member data
int SGERG88:: checkInputs(double p,double t,double Hs,double d,double xCO2,double xH2){
    //d:相对密度（相对于空气），Hs：高位发热量
    if(p>12.0 || p<=0.0) return -11;
    if(t+273.15>338 || t+273.15<263) return -12;
    if(xCO2>0.3 || xCO2<-1E-8) return -15;
    if(xH2>0.1 || xH2<-1E-8) return -16;
    if(Hs>48 || Hs<20) return -13;
    if(d>0.9 || d<0.55) return -14;    
    return 0;
}

int SGERG88::calcIntermediate(double Hs,double d,double xCO2,double xH2,
                          double& xCH,double& xN2,double& xCO,double& HCH,double& MCH,
                          double& Bn,double& rhomn,double& rhon){
//member function:Calculation of intermediate data
    //x1==xCH,x2==xN2,x3==xCO2,x4==xH2,x5==xCO
    //d:relative density (comparing to air?)
    double err=1E-6;
    rhon=d*Rhonair;
    xCO=0.0964*xH2;
    double xM_right=(xCO2*MCO2+xH2*MH2+xCO*MCO);
    double T=273.15;//T==Tn==273.15
    double B11_C0=b_H0[0]+b_H0[1]*T+b_H0[2]*T*T;
    double B11_C1=b_H1[0]+b_H1[1]*T+b_H1[2]*T*T;
    double B11_C2=b_H2[0]+b_H2[1]*T+b_H2[2]*T*T;
    double B22=b_22[0]+b_22[1]*T+b_22[2]*T*T;
    double B33=b_33[0]+b_33[1]*T+b_33[2]*T*T;
    double B44=b_44[0]+b_44[1]*T+b_44[2]*T*T;
    double B55=b_55[0]+b_55[1]*T+b_55[2]*T*T;
    double B14=b_14[0]+b_14[1]*T+b_14[2]*T*T;
    double B15=b_15[0]+b_15[1]*T+b_15[2]*T*T;
    double B23=b_23[0]+b_23[1]*T+b_23[2]*T*T;
    double B24=b_24[0]+b_24[1]*T+b_24[2]*T*T;
    Bn=-0.065;
    rhomn=1.0/(Vmnid+Bn);
    HCH=1000.0;
    double itera0=0;
    int itera1=0;
    double Hs_v;
    while(1){//outer loop
        MCH=-2.709328+0.021062199*HCH;
        xCH=Hs/(HCH*rhomn)-(xH2*HH2+xCO*HCO)/HCH;
        xN2=1.0-xCH-xCO-xCO2-xH2;
        double rhon_u2=(xCH*MCH+xN2*MN2+xM_right)*rhomn;
        if(abs(rhon_u2-rhon)>err){//inner loop
            double HCH0=HCH;
            double rhon_u0=rhon_u2;
            HCH=HCH0+1.0;
            MCH=-2.709328+0.021062199*HCH;
            xCH=Hs/(HCH*rhomn)-(xH2*HH2+xCO*HCO)/HCH;
            xN2=1.0-xCH-xCO-xCO2-xH2;
            double rhon_u1=(xCH*MCH+xN2*MN2+xM_right)*rhomn;
            double HCH1=HCH;
            HCH=HCH1+(rhon-rhon_u1)/(rhon_u0-rhon_u1)*(HCH0-HCH1);
            MCH=-2.709328+0.021062199*HCH;
            xCH=Hs/(HCH*rhomn)-(xH2*HH2+xCO*HCO)/HCH;
            xN2=1.0-xCH-xCO-xCO2-xH2;
            rhon_u2=(xCH*MCH+xN2*MN2+xM_right)*rhomn;
            while(abs(rhon_u2-rhon)>err && itera1<10){
                itera1++;
                HCH0=HCH1;
                rhon_u0=rhon_u1;
                HCH1=HCH;
                rhon_u1=rhon_u2;
                HCH=HCH1+(rhon-rhon_u1)/(rhon_u0-rhon_u1)*(HCH0-HCH1);
                MCH=-2.709328+0.021062199*HCH;
                xCH=Hs/(HCH*rhomn)-(xH2*HH2+xCO*HCO)/HCH;
                xN2=1.0-xCH-xCO-xCO2-xH2;
                rhon_u2=(xCH*MCH+xN2*MN2+xM_right)*rhomn;
            }
            if(abs(rhon_u2-rhon)>err && itera1>=10) return -1;//inner loop not convergency
        }
        double B11=B11_C0+B11_C1*HCH+B11_C2*HCH*HCH;
        double B12=(0.72+1.875E-5*(320.0-T)*(320.0-T))*(B11+B22)/2.0;
        double B13=-0.865*sqrt(B11*B33);
        Bn=B11*xCH*xCH+2*B12*xCH*xN2+2*B13*xCH*xCO2+2*B14*xCH*xH2+2*B15*xCH*xCO
            +B22*xN2*xN2+2*B23*xN2*xCO2+2*B24*xN2*xH2
            +B33*xCO2*xCO2+B44*xH2*xH2+B55*xCO*xCO;
        rhomn=1.0/(Vmnid+Bn);
        Hs_v=(xCH*HCH+xH2*HH2+xCO*HCO)*rhomn;
        itera0++;
        if(abs(Hs-Hs_v)<=err || itera0>=10) break;//outer loop end
    }
    if(abs(Hs-Hs_v)>err && itera0>=10) return -2;//outer loop not convergency
    return 0;
}

int SGERG88::ZRho_SGERG88(double p,double t,double Hs,double d,double xCO2,double xH2,
                      double& Z,double& rhom,double& rho){
//public function:calculation Z,rhom,rho from inputs p,t,H2,d,xCO2,xH2
    double err=1E-6;
    double xCH,xN2,xCO;
    double HCH,MCH;
    double Bn;
    double rhomn,rhon;
    int ret=checkInputs(p,t,Hs,d,xCO2,xH2);
    if(ret<0) return ret;//check inputs invalid
    ret=calcIntermediate(Hs,d,xCO2,xH2,xCH,xN2,xCO,HCH,MCH,Bn,rhomn,rhon);
    if(ret<0) return ret;//calculation of intermediate is not convergency,failure and exit
    double T=t+273.15;    
    //the second virial coefficients:
        double B11_C0=b_H0[0]+b_H0[1]*T+b_H0[2]*T*T;
        double B11_C1=b_H1[0]+b_H1[1]*T+b_H1[2]*T*T;
        double B11_C2=b_H2[0]+b_H2[1]*T+b_H2[2]*T*T;
        double B22=b_22[0]+b_22[1]*T+b_22[2]*T*T;
        double B33=b_33[0]+b_33[1]*T+b_33[2]*T*T;
        double B44=b_44[0]+b_44[1]*T+b_44[2]*T*T;
        double B55=b_55[0]+b_55[1]*T+b_55[2]*T*T;
        double B14=b_14[0]+b_14[1]*T+b_14[2]*T*T;
        double B15=b_15[0]+b_15[1]*T+b_15[2]*T*T;
        double B23=b_23[0]+b_23[1]*T+b_23[2]*T*T;
        double B24=b_24[0]+b_24[1]*T+b_24[2]*T*T;
        double B11=B11_C0+B11_C1*HCH+B11_C2*HCH*HCH;
        double B12=(0.72+1.875E-5*(320.0-T)*(320.0-T))*(B11+B22)/2.0;
        double B13=-0.865*sqrt(B11*B33);
    double B=B11*xCH*xCH+2*B12*xCH*xN2+2*B13*xCH*xCO2+2*B14*xCH*xH2+2*B15*xCH*xCO
            +B22*xN2*xN2+2*B23*xN2*xCO2+2*B24*xN2*xH2
            +B33*xCO2*xCO2+B44*xH2*xH2+B55*xCO*xCO;
    //the third virial coefficients:
        double C111_C0=c_H0[0]+c_H0[1]*T+c_H0[2]*T*T;
        double C111_C1=c_H1[0]+c_H1[1]*T+c_H1[2]*T*T;
        double C111_C2=c_H2[0]+c_H2[1]*T+c_H2[2]*T*T;
        double C222=c_222[0]+c_222[1]*T+c_222[2]*T*T;
        double C333=c_333[0]+c_333[1]*T+c_333[2]*T*T;
        double C444=c_444[0]+c_444[1]*T+c_444[2]*T*T;
        double C115=c_115[0]+c_115[1]*T+c_115[2]*T*T;
        double C223=c_223[0]+c_223[1]*T+c_223[2]*T*T;
        double C233=c_233[0]+c_233[1]*T+c_233[2]*T*T;
        double y112=0.92+0.0013*(T-270.0);//y112==y122
        double y113=0.92;//y113==y133
        double y114=1.20;
        double y123=1.10;
        double C111=C111_C0+C111_C1*HCH+C111_C2*HCH*HCH;
        double C112=y112*pow(C111*C111*C222,1.0/3.0);
        double C113=y113*pow(C111*C111*C333,1.0/3.0);
        double C114=y114*pow(C111*C111*C444,1.0/3.0);
        double C122=y112*pow(C111*C222*C222,1.0/3.0);//y112==y122
        double C123=y123*pow(C111*C222*C333,1.0/3.0);
        double C133=y113*pow(C111*C333*C333,1.0/3.0);//y113==y133
    double C=C111*xCH*xCH*xCH+3*C112*xCH*xCH*xN2+3*C113*xCH*xCH*xCO2+3*C114*xCH*xCH*xH2
            +3*C115*xCH*xCH*xCO+3*C122*xCH*xN2*xN2+6*C123*xCH*xN2*xCO2+3*C133*xCH*xCO2*xCO2
            +C222*xN2*xN2*xN2+3*C223*xN2*xN2*xCO2+3*C233*xN2*xCO2*xCO2
            +C333*xCO2*xCO2*xCO2+C444*xH2*xH2*xH2;
    rhom=1.0/(R*T/p+B);
    double p_w=R*T*rhom*(1.0+B*rhom+C*rhom*rhom);
    int itera=0;
    while(abs(p-p_w)>err && itera<10){
        itera++;
        rhom=1.0/((R*T/p)*(1+B*rhom+C*rhom*rhom));
        p_w=R*T*rhom*(1.0+B*rhom+C*rhom*rhom);
    }
    if(abs(p-p_w)>err && itera>=10) return -3;//the third loop not convergency,failure and exit
    Z=1.0+B*rhom+C*rhom*rhom;//outputs:Z and rhom
    //the third virial coefficients at normal conditions:Tn=273.15K,pn=0.101325Mpa
    double Tn=273.15;
    double pn=0.101325;
        C111_C0=c_H0[0]+c_H0[1]*Tn+c_H0[2]*Tn*Tn;
        C111_C1=c_H1[0]+c_H1[1]*Tn+c_H1[2]*Tn*Tn;
        C111_C2=c_H2[0]+c_H2[1]*Tn+c_H2[2]*Tn*Tn;
        C222=c_222[0]+c_222[1]*Tn+c_222[2]*Tn*Tn;
        C333=c_333[0]+c_333[1]*Tn+c_333[2]*Tn*Tn;
        C444=c_444[0]+c_444[1]*Tn+c_444[2]*Tn*Tn;
        C115=c_115[0]+c_115[1]*Tn+c_115[2]*Tn*Tn;
        C223=c_223[0]+c_223[1]*Tn+c_223[2]*Tn*Tn;
        C233=c_233[0]+c_233[1]*Tn+c_233[2]*Tn*Tn;
        y112=0.92+0.0013*(Tn-270.0);//y112==y122
        y113=0.92;//y113==y133
        y114=1.20;
        y123=1.10;
        C111=C111_C0+C111_C1*HCH+C111_C2*HCH*HCH;
        C112=y112*pow(C111*C111*C222,1.0/3.0);
        C113=y113*pow(C111*C111*C333,1.0/3.0);
        C114=y114*pow(C111*C111*C444,1.0/3.0);
        C122=y112*pow(C111*C222*C222,1.0/3.0);//y112==y122
        C123=y123*pow(C111*C222*C333,1.0/3.0);
        C133=y113*pow(C111*C333*C333,1.0/3.0);//y113==y133
    double Cn=C111*xCH*xCH*xCH+3*C112*xCH*xCH*xN2+3*C113*xCH*xCH*xCO2+3*C114*xCH*xCH*xH2
             +3*C115*xCH*xCH*xCO+3*C122*xCH*xN2*xN2+6*C123*xCH*xN2*xCO2+3*C133*xCH*xCO2*xCO2
             +C222*xN2*xN2*xN2+3*C223*xN2*xN2*xCO2+3*C233*xN2*xCO2*xCO2
             +C333*xCO2*xCO2*xCO2+C444*xH2*xH2*xH2;
    double Zn=1.0+Bn*rhomn+Cn*rhomn*rhomn;
    rho=d*Rhonair*p*Zn*Tn/(pn*Z*T);//output:rho    
    return 0;//calculation success
}

int main(){
    double xCO2,xH2,d,Hs,p,t;
    double Z,rho,rhom;
    xCO2=0.006;
    xH2=0.0;
    d=0.581;
    Hs=40.66;
    p=6;
    t=-3.15;
    SGERG88::ZRho_SGERG88(p,t,Hs,d,xCO2,xH2,Z,rhom,rho);
    t=6.85;
    SGERG88::ZRho_SGERG88(p,t,Hs,d,xCO2,xH2,Z,rhom,rho);
    t=16.85;
    SGERG88::ZRho_SGERG88(p,t,Hs,d,xCO2,xH2,Z,rhom,rho);
    t=36.85;
    SGERG88::ZRho_SGERG88(p,t,Hs,d,xCO2,xH2,Z,rhom,rho);
    t=56.85;
    SGERG88::ZRho_SGERG88(p,t,Hs,d,xCO2,xH2,Z,rhom,rho);
    return 0;
}