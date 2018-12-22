#include "AGA892DC2.h"
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
    double err=1E-9;
    double BI[18];
    double K,F,Q,G,U,E;
    double CNS[46];
    K=0;F=0;Q=0;G=0;U=0;E=0;
    for(int i=0;i<18;i++)
        BI[i]=0;
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        K+=xi*pow(Ki[ii-1],2.5);
        U+=xi*pow(Ei[ii-1],2.5);
        G+=xi*Gi[ii-1];
        Q+=xi*Qi[ii-1];
        F+=xi*xi*Fi[ii-1];
        E+=xi*Ei[ii-1];
    }
    K=K*K;
    U=U*U;
    for(int i=0;i<NUM_x-1;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        if(abs(xi)>1E-8 && ii>=1 && ii<=8){
            for(int j=i+1;j<NUM_x;j++){
                int jj=ID[j];
                double xij=xi*xi_raw[Index2raw[j]];
                if(abs(xij)>1E-8 && jj<=19){
                    K+=2.0*xij*(pow(Kij0[ii-1][jj-1],5.0)-1.0)*pow(Ki[ii-1]*Ki[jj-1],2.5);
                    U+=2.0*xij*(pow(Uij0[ii-1][jj-1],5.0)-1.0)*pow(Ei[ii-1]*Ei[jj-1],2.5);
                    G+=xij*(Gij0[ii-1][jj-1]-1.0)*(Gi[ii-1]+Gi[jj-1]);
                }
            }
        }
    }
    for(int i=0;i<NUM_x;i++){
        int ii=ID[i];
        double xi=xi_raw[Index2raw[i]];
        if(abs(xi)>1E-8){
            for(int j=i;j<NUM_x;j++){
                int jj=ID[j];
                double xij=xi*xi_raw[Index2raw[j]];
                if(abs(xij)>1E-14){
                    if(ii!=jj) xij*=2;
                    double Eij=Eij0[ii-1][jj-1]*sqrt(Ei[ii-1]*Ei[jj-1]);
                    double Gij=Gij0[ii-1][jj-1]*(Gi[ii-1]+Gi[jj-1])/2.0;
                    for(int k=0;k<18;k++){
                        double BN=pow(Gij+1.0-gn[k],1.0*gn[k])
                            *pow(Qi[ii-1]*Qi[jj-1]+1.0-qn[k],qn[k])
                            *pow(sqrt(Fi[ii-1]*Fi[jj-1])+1.0-fn[k],fn[k])
                            *pow(Si[ii-1]*Si[jj-1]+1.0-sn[k],sn[k])
                            *pow(Wi[ii-1]*Wi[jj-1]+1.0-wn[k],wn[k]);
                        BI[k]+=an[k]*xij*pow(Eij,un[k])*pow(Ki[ii-1]*Ki[jj-1],1.5)*BN;
                    }
                }
            }
        }
    }
    K=pow(K,0.2);
    U=pow(U,0.2); 
    for(int i=12;i<58;i++){
        CNS[i-12]=an[i]*pow(G+1.0-gn[i],gn[i])*pow(Q*Q+1.0-qn[i],qn[i])*pow(F+1.0-fn[i],fn[i])*pow(U,un[i]);   
    }
    rhom=2.0;
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
