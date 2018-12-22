#include "AGA892DC.h"
#include <iostream>
#include <iomanip>

using namespace std;

int AGA892DC::checkInputs(double xi_raw[],int ID_raw[],int NUM_x,double p,double t){
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
    double maxXi[21]={
        1.00000,0.70000,0.20000,0.10000,0.03500,0.00015,0.00015,0.10000,0.03000,0.00015,0.01500,
        0.01500,0.00500,0.00500,0.00100,0.00050,0.00050,0.00050,0.00050,0.00500,0.00015
    };
    if(p<=0 || p>12) return -100;//压力超出范围
    if((t+273.15)<263 || (t+273.15)>338) return -200;//温度超出范围
    double sum_0=0.0;
    for(int i=0;i<NUM_x;i++){
        int ii=ID_raw[i];
        if(xi_raw[i]>maxXi[ii-1] || xi_raw[i]<-1E-8) return -ii;//哪个组分输入超限，返回该组分标识的负值
        if(ii==1 && xi_raw[i]<0.7) return -ii;//CH4<0.7,超下限 
        sum_0+=xi_raw[i];
    }
    if(abs(sum_0-1.0)>0.0001)
        return 0;//所有组分摩尔分数之和必须1+-0.0001
    return 1;
}
void AGA892DC::sortInput(double xi_raw[],int ID_raw[],int ID[],int Index2raw[],int NUM_x){
    if(status<1) return;
    ID[0]=ID_raw[0];
    Index2raw[0]=0;
    for(int i=1;i<NUM_x;i++){
        int key=ID_raw[i];
        int j=i-1;
        while(j>=0 && key<ID[j]){
            ID[j+1]=ID[j];
            Index2raw[j+1]=Index2raw[j];
            j--;
        }
        ID[j+1]=key;
        Index2raw[j+1]=i;
    }
    status=2;
}
void AGA892DC::calcSecondVirial(double xi_raw[],int ID[],int Index2raw[],int NUM_x,double BI[18]){
    if(status<2) return;
    for(int n=0;n<18;n++){
        BI[n]=0.0;
    }
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        int Fi_sqrt=0;
        int W_i=0;
        double Q_i=0.0;
        double S_i=0.0;
        double K_i=Ki[ii-1];
        if(ii==3)
            Q_i=Qi3;
        else if(ii==6){
            Q_i=Qi6;
            S_i=Si6;
            W_i=1;                
        }else if(ii==7){
            Q_i=Qi7;
            S_i=Si7;  
        }else if(ii==8)
            Fi_sqrt=1;
        double G_i=Gi[ii-1];
        double E_i=Ei[ii-1];
        for(int j=i;j<NUM_x;j++){
            double Gij;
            double Eij;
            int jj=ID[j];            
            double xj=xi_raw[Index2raw[j]];
            if(i!=j) xj=2*xj;//i，j对称
            int Fj_sqrt=0;
            int W_j=0;
            double Q_j=0.0;
            double S_j=0.0;
            double K_j=Ki[jj-1];
            if(jj==3)
                Q_j=Qi3;
            else if(jj==6){
                Q_j=Qi6;
                S_j=Si6;
                W_j=1;                
            }else if(jj==7){
                Q_j=Qi7;
                S_j=Si7;  
            }else if(jj==8)
                Fj_sqrt=1;

            double G_j=Gi[jj-1];
            double E_j=Ei[jj-1];
            double E_ij=1;     
            double G_ij=1;    
            if(ii==1 && jj>=2 && jj<=19){
                E_ij=E1j[jj-2];
                if(jj==3)
                    G_ij=G13;
                else if(jj==8)
                    G_ij=G18;
            }else if(ii==2 && jj>=3 && jj<=14){
                E_ij=E2j[jj-3];
                if(jj==3)
                    G_ij=G23;
            }else if(ii==3 && jj>=4 && jj<=19){
                E_ij=E3j[jj-4];
                if(jj<=6)
                    G_ij=G3j[jj-4];
            }else if(ii==4 && jj>=5 && jj<=14){
                E_ij=E4j[jj-5];
            }else if(ii==5){
                if(jj==8)
                    E_ij=E58;
                else if(jj==12)
                    E_ij=E512;
            }else if(ii==7 && jj>=15 && jj<=19){
                E_ij=E7j[jj-15];
            }else if(ii==8){
                if(jj==9)
                    E_ij=E89;
                else if(jj==11)
                    E_ij=E811;
                else if(jj==12)
                    E_ij=E812; 
            }
            Gij=G_ij*(G_i+G_j)/2.0;
            Eij=E_ij*sqrt(E_i*E_j);
            for(int n=0;n<18;n++){
                double Bnij=pow(Gij+1.0-gn[n],gn[n])
                    *pow(Q_i*Q_j+1.0-qn[n],qn[n])
                    *pow(1.0+Fi_sqrt*Fj_sqrt-fn[n],fn[n])
                    *pow(S_i*S_j+1.0-sn[n],sn[n])
                    *pow(1.0+W_i*W_j-wn[n],wn[n]);
                BI[n]+=an[n]*xi*xj*pow(Eij,un[n])*pow(K_i*K_j,1.5)*Bnij;
            }                
        }
    }
    status=3;
}
void AGA892DC::calcKFQGU(double xi_raw[],int ID[],int Index2raw[],int NUM_x,
                        double& K,double& F,double& Q,double& G,double& U){
    if(status<2) return;
    K=0.0;
    F=0.0;
    Q=0.0;
    G=0.0;
    U=0.0;
    double leftK=0.0;
    double leftU=0.0;
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        leftK+=xi*pow(Ki[ii-1],2.5);
        leftU+=xi*pow(Ei[ii-1],2.5);
        G+=xi*Gi[ii-1];
        if(ii==8) F=xi*xi;
        if(ii==3){
            Q+=xi*Qi3;
        }else if(ii==6){
            Q+=xi*Qi6;
        }else if(ii==7){
            Q+=xi*Qi7;
        }
    }
    double rightK=0.0;
    double rightU=0.0;
    for(int i=0;i<NUM_x-1;i++){
        int ii=ID[i];
        if(ii>=1 && ii<=3){
            double xi=xi_raw[Index2raw[i]];
            double G_i=Gi[ii-1];            
            for(int j=i+1;j<NUM_x;j++){
                int jj=ID[j];
                double G_ij=1.0;
                if(ii==1){
                    if(jj==3)
                        G_ij=G13;
                    else if(jj==8){
                        G_ij=G18;
                    }
                }else if(ii==2 && jj==3){
                    G_ij=G23;
                }else if(ii==3){
                    if(jj==4){
                        G_ij=G3j[0];
                    }else if(jj==6){
                        G_ij=G3j[2];
                    }
                }            
                if(abs(G_ij-1.0)>1E-7){//G_ij!=1.0
                    double xj=xi_raw[Index2raw[j]];
                    double G_j=Gi[jj-1];
                    G+=xi*xj*(G_ij-1.0)*(G_i+G_j);
                }//else G+=0;    
            }
        }//else G+=0;
    }
    for(int i=0;i<NUM_x-1;i++){
        int ii=ID[i];
        if((ii>=1 && ii<=4) || ii==7){
            double xi=xi_raw[Index2raw[i]];
            double K_i=Ki[ii-1];
            double E_i=Ei[ii-1];
            for(int j=i+1;j<NUM_x;j++){
                int jj=ID[j];
                double K_ij=1.0;
                double U_ij=1.0;
                if(ii==1 && jj<=19){
                    K_ij=K1j[jj-2];
                    U_ij=U1j[jj-2];
                }else if(ii==2 && jj<=12){
                    if(jj<=8) K_ij=K2j[jj-3];
                    U_ij=U2j[jj-3];
                }else if(ii==3 && jj<=19){
                    K_ij=K3j[jj-4];
                    U_ij=U3j[jj-4];
                }else if(ii==4 && jj<=14){
                    if(jj<=8) K_ij=K4j[jj-5];
                    U_ij=U4j[jj-5];
                }else if(ii==7 && jj>=15 && jj<=19){
                    K_ij=K7j[jj-15];
                    U_ij=U7j[jj-15];
                }
                double xj=xi_raw[Index2raw[j]];
                if(abs(K_ij-1.0)>1E-8){//K_ij!=1.0
                    double K_j=Ki[jj-1];
                    rightK+=xi*xj*(K_ij*K_ij*K_ij*K_ij*K_ij-1.0)*pow(K_i*K_j,2.5);                    
                }//else rightK+=0; 
                if(abs(U_ij-1.0)>1E-8){//U_ij!=1.0
                    double E_j=Ei[jj-1];
                    rightU+=xi*xj*(U_ij*U_ij*U_ij*U_ij*U_ij-1.0)*pow(E_i*E_j,2.5);                   
                }//else rightU+=0; 
            }
        }//else right+=0;
    }
    K=pow(leftK*leftK+2*rightK,0.2);
    U=pow(leftU*leftU+2*rightU,0.2);
    status=4;    
}
double AGA892DC::rhom2PZ(double rhom,double tt,double BMIX,double CNS[46],double K,
                        double F,double Q,double G,double U,double& Z){    
    double rhor=K*K*K*rhom;
    double T=tt+273.15;
    double R=8.314510E-3;
    double TUN;
    double mid=0;
    double right=0;
    for(int i=12;i<18;i++){
        TUN=pow(T,un[i]);
        mid+=rhor*CNS[i-12]/TUN;
        right+=CNS[i-12]/TUN*(bn[i]-cn[i]*kn[i]*pow(rhor,kn[i]))*pow(rhor,bn[i])
                *exp(-cn[i]*pow(rhor,kn[i]));
    }
    for(int i=18;i<58;i++){
        TUN=pow(T,un[i]);
        right+=CNS[i-12]/TUN*(bn[i]-cn[i]*kn[i]*pow(rhor,kn[i]))*pow(rhor,bn[i])
                *exp(-cn[i]*pow(rhor,kn[i]));
    }
    Z=1.0+BMIX*rhom-mid+right;
    return rhom*R*T*Z;
}
int AGA892DC::calcZRho(double xi_raw[],int ID_raw[],int NUM_x,double p,double t,
                        double& Z,double& rhom,double& rho){
    //xi:摩尔分数,p:绝对压力,MPa,t:摄氏度，C
    status=checkInputs(xi_raw,ID_raw,NUM_x,p,t);
    if(status<1) return status;//inputs error
    int ID[21]={};
    int Index2raw[21]={};
    if(status<2)
        sortInput(xi_raw,ID_raw,ID,Index2raw,NUM_x);
    double err=1E-8;
    double BI[18];
    calcSecondVirial(xi_raw,ID,Index2raw,NUM_x,BI);
    double K,F,Q,G,U;
    calcKFQGU(xi_raw,ID,Index2raw,NUM_x,K,F,Q,G,U);
    double T=t+273.15;
    double R=8.314510E-3;
    double TUN;
    double BMIX=0;
    for(int i=0;i<18;i++){
        TUN=pow(T,un[i]);
        BMIX+=BI[i]/TUN;
    }
    double CNS[46];
    for(int i=12;i<58;i++){
        CNS[i-12]=an[i]*pow(G+1.0-gn[i],gn[i])*pow(Q*Q+1.0-qn[i],qn[i])*pow(F+1.0-fn[i],fn[i])*pow(U,un[i]);
    }
    rhom=2.0;//p/(0.91501*(t+273.15)*(8.314510E-3));
    double pp=rhom2PZ(rhom,t,BMIX,CNS,K,F,Q,G,U,Z);
    int itera=0;
    if(abs(pp-p)>err){
        double rhom0=rhom;
        double p0=pp;
        double rhom1=rhom+1.0;
        double p1=rhom2PZ(rhom1,t,BMIX,CNS,K,F,Q,G,U,Z);
        rhom=rhom1+(p-p1)/(p0-p1)*(rhom0-rhom1);
        pp=rhom2PZ(rhom,t,BMIX,CNS,K,F,Q,G,U,Z);
        while(abs(pp-p)>err && itera<=10){
            itera++;
            rhom0=rhom1;
            p0=p1;
            rhom1=rhom;
            p1=pp;
            rhom=rhom1+(p-p1)/(p0-p1)*(rhom0-rhom1);
            pp=rhom2PZ(rhom,t,BMIX,CNS,K,F,Q,G,U,Z);
        }
    }
    if(itera>10) return -3;//iteration exceed 10,not convergency
    double M=0.0;
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        M+=Mi[i]*xi_raw[Index2raw[i]];
    }
    rho=rhom*M;
    status=100; 
    return status;//calculation success
}


int main(int argc, char const *argv[])
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
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    t=6.85;
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    t=16.85;
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    t=36.85;
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    t=56.85;
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    p=12.0;
    t=-3.15;
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    t=6.85;
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    t=16.85;
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    t=36.85;
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    t=56.85;
    AGA892DC::calcZRho(xi_raw,ID_raw,NUM_x,p,t,Z,rhom,rho);
    cout<<setprecision(10)<<Z<<endl;
    return 0;
}
