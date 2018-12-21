#include "AGA892DC.h"
#include <iostream>
#include <iomanip>

using namespace std;

void AGA892DC::sortInput(){
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
void AGA892DC::calcSecondVirial(double BI[18]){
    if(status<1) return;
    if(status<2)
        sortInput();
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
double AGA892DC::calcK(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double K=0.0;
    double left=0.0;
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        double K_i=Ki[ii-1];
        left+=xi*pow(K_i,2.5);
    }
    double right=0.0;
    for(int i=0;i<NUM_x-1;i++){
        int ii=ID[i];
        if((ii>=1 && ii<=4) || ii==7){
            double xi=xi_raw[Index2raw[i]];
            double K_i=Ki[ii-1];
            for(int j=i+1;j<NUM_x;j++){
                int jj=ID[j];
                double K_ij=1.0;
                if(ii==1 && jj<=19){
                    K_ij=K1j[jj-2];
                }else if(ii==2 && jj<=8){
                    K_ij=K2j[jj-3];
                }else if(ii==3 && jj<=19){
                    K_ij=K3j[jj-4];
                }else if(ii==4 && jj<=8){
                    K_ij=K4j[jj-5];
                }else if(ii==7 && jj>=15 && jj<=19){
                    K_ij=K7j[jj-15];
                }
                //if(abs(K_ij-1.0)>1E-7){//K_ij!=1.0
                    double xj=xi_raw[Index2raw[j]];
                    double K_j=Ki[jj-1];
                    right+=xi*xj*(K_ij*K_ij*K_ij*K_ij*K_ij-1.0)*pow(K_i*K_j,2.5);                    
                //}//else right+=0; 
            }
        }//else right+=0;
    }
    K=pow(left*left+2*right,0.2);
    status=4;
    return K;
}
double AGA892DC::calcF(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double F=0.0;
    status=5;
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        if(ii==8){
            F=xi_raw[Index2raw[i]];
            return F*F;
        }
    }
    return F;
}
double AGA892DC::calcQ(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double Q=0.0;
    status=6;
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        if(ii==3){
            Q+=xi_raw[Index2raw[i]]*Qi3;
        }else if(ii==6){
            Q+=xi_raw[Index2raw[i]]*Qi6;
        }else if(ii==7){
            Q+=xi_raw[Index2raw[i]]*Qi7;
        }
    }
    return Q;
}
double AGA892DC::calcG(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double G=0.0;    
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        double G_i=Gi[ii-1];
        G+=xi*G_i;
    }
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
    status=7;
    return G;
}
double AGA892DC::calcU(){
    if(status<1) return -999999.0;
    if(status<2)
        sortInput();
    double U=0.0;
    double left=0.0;        
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        double E_i=Ei[ii-1];
        left+=xi*pow(E_i,2.5);
    }    
    double right=0.0;
    for(int i=0;i<NUM_x-1;i++){
        int ii=ID[i];
        if((ii>=1 && ii<=4) || ii==7){
            double xi=xi_raw[Index2raw[i]];
            double E_i=Ei[ii-1];
            for(int j=i+1;j<NUM_x;j++){
                int jj=ID[j];
                double U_ij=1.0;
                if(ii==1 && jj<=19){
                    U_ij=U1j[jj-2];
                }else if(ii==2 && jj<=12){
                    U_ij=U2j[jj-3];
                }else if(ii==3 && jj<=19){
                    U_ij=U3j[jj-4];
                }else if(ii==4 && jj<=14){
                    U_ij=U4j[jj-5];
                }else if(ii==7 && jj>=15 && jj<=19){
                    U_ij=U7j[jj-15];
                }
                if(abs(U_ij-1.0)>1E-7){//U_ij!=1.0
                    double xj=xi_raw[Index2raw[j]];
                    double E_j=Ei[jj-1];
                    right+=xi*xj*(U_ij*U_ij*U_ij*U_ij*U_ij-1.0)*pow(E_i*E_j,2.5);
                }//else right+=0; 
            }
        }//else right+=0;        
    }
    U=pow(left*left+2.0*right,0.2);
    status=8;
    return U;
}
double AGA892DC::rhom2PZ(double rhom,double tt,double BI[18],double CNS[46],double K,
                        double F,double Q,double G,double U,double& Z){    
    double rhor=K*K*K*rhom;
    double T=tt+273.15;
    double R=8.314510E-3;
    double BMIX=0;
    double TUN;
    for(int i=0;i<12;i++){
        TUN=pow(T,un[i]);
        BMIX+=BI[i]/TUN;
    }
    double mid=0;
    double right=0;
    for(int i=12;i<18;i++){
        TUN=pow(T,un[i]);
        BMIX+=BI[i]/TUN;
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
void AGA892DC::calcZRhom(double& Z,double& rhom){
    if(status<1) return;
    if(status<2)
        sortInput();
    double err=1E-6;
    double BI[18];
    calcSecondVirial(BI);
    double K=calcK();
    double F=calcF();
    double Q=calcQ(); 
    double G=calcG();  
    double U=calcU();
    double CNS[46];
    for(int i=12;i<58;i++){
        CNS[i-12]=an[i]*pow(G+1.0-gn[i],gn[i])*pow(Q*Q+1.0-qn[i],qn[i])*pow(F+1.0-fn[i],fn[i])*pow(U,un[i]);
    }
    rhom=2.0;//p/(0.91501*(t+273.15)*(8.314510E-3));
    double pp=rhom2PZ(rhom,t,BI,CNS,K,F,Q,G,U,Z);
    if(abs(pp-p)>err){
        double rhom0=rhom;
        double p0=pp;
        double rhom1=rhom+1.0;
        double p1=rhom2PZ(rhom1,t,BI,CNS,K,F,Q,G,U,Z);
        rhom=rhom1+(p-p1)/(p0-p1)*(rhom0-rhom1);
        pp=rhom2PZ(rhom,t,BI,CNS,K,F,Q,G,U,Z);
        while(abs(pp-p)>err){
            rhom0=rhom1;
            p0=p1;
            rhom1=rhom;
            p1=pp;
            rhom=rhom1+(p-p1)/(p0-p1)*(rhom0-rhom1);
            pp=rhom2PZ(rhom,t,BI,CNS,K,F,Q,G,U,Z);
        }
    }
    status=10;
}


int main(int argc, char const *argv[])
{
    double rhom,Z;
    double pp=6.0;
    AGA892DC::p=pp;
    AGA892DC::t=-3.15;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    AGA892DC::p=pp;
    AGA892DC::t=6.85;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    AGA892DC::p=pp;
    AGA892DC::t=16.85;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    AGA892DC::p=pp;
    AGA892DC::t=36.85;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    AGA892DC::p=pp;
    AGA892DC::t=56.85;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    pp=12.0;
    AGA892DC::p=pp;
    AGA892DC::t=-3.15;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    AGA892DC::p=pp;
    AGA892DC::t=6.85;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    AGA892DC::p=pp;
    AGA892DC::t=16.85;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    AGA892DC::p=pp;
    AGA892DC::t=36.85;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    AGA892DC::p=pp;
    AGA892DC::t=56.85;
    AGA892DC::calcZRhom(Z,rhom);
    cout<<setprecision(10)<<Z<<endl;
    return 0;
}
