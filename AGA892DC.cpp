#include "AGA892DC.h"
//const static member datas:
const int AGA892DC::bn[58]={
    1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,
    3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5,5,5,6,6,7,7,8,8,8,9,9
};
const int AGA892DC::cn[58]={
    0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,1,0,1,
    1,1,1,1,1,1,1,1,0,0,1,1,1,1,1,0,1,1,1,1,0,1,0,1,1,1,1,1,1
};
const int AGA892DC::kn[58]={
    0,0,0,0,0,0,0,0,0,0,0,0,3,2,2,2,4,4,0,0,2,2,2,4,4,4,4,0,1,
    1,2,2,3,3,4,4,4,0,0,2,2,2,4,4,0,2,2,4,4,0,2,0,2,1,2,2,2,2
};
const int AGA892DC::gn[58]={
    0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,1,
    0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,0,0
};
const int AGA892DC::qn[58]={
    0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,
    0,0,0,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1,0,1,0,0,1,0,0,0,0,0,1
};
const int AGA892DC::fn[58]={
    0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,
    1,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
const int AGA892DC::sn[58]={
    0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
const int AGA892DC::wn[58]={
    0,0,0,0,0,0,0,0,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
};
const double AGA892DC::an[58]={
  0.153832600,
  1.341953000,
 -2.998583000,
 -0.048312280,
  0.375796500,
 -1.589575000,
 -0.053588470,
  0.886594630,
 -0.710237040,
 -1.471722000,
  1.321850350,
 -0.786659250,
  2.291290000E-09,
  0.157672400,
 -0.436386400,
 -0.044081590,
 -0.003433888,
  0.032059050,
  0.024873550,
  0.073322790,
 -0.001600573,
  0.642470600,
 -0.416260100,
 -0.066899570,
  0.279179500,
 -0.696605100,
 -0.002860589,
 -0.008098836,
  3.150547000,
  0.007224479,
 -0.705752900,
  0.534979200,
 -0.079314910,
 -1.418465000,
 -5.999050000E-17,
  0.105840200,
  0.034317290,
 -0.007022847,
  0.024955870,
  0.042968180,
  0.746545300,
 -0.291961300,
  7.294616000,
 -9.936757000,
 -0.005399808,
 -0.243256700,
  0.049870160,
  0.003733797,
  1.874951000,
  0.002168144,
 -0.658716400,
  0.000205518,
  0.009776195,
 -0.020487080,
  0.015573220,
  0.006862415,
 -0.001226752,
  0.002850908
};
const double AGA892DC::un[58]={
    0,0.5,1,3.5,-0.5,4.5,0.5,7.5,9.5,6,12,12.5,-6,2,3,2,2,11,-0.5,0.5,0,4,6,21,23,22,
    -1,-0.5,7,-1,6,4,1,9,-13,21,8,-0.5,0,2,7,9,22,23,1,9,3,8,23,1.5,5,-0.5,4,7,3,0,1,0
};
const double AGA892DC::Mi[21]={
    16.0430,28.0135,44.0100,30.0700,44.0970,18.0153,34.0820,2.0159,28.0100,31.9988,58.1230,
    58.1230,72.1500,72.1500,86.1770,100.2040,114.2310,128.2580,142.2850,4.0026,39.9480
};
const double AGA892DC::Ei[21]={
    151.31830,99.73778,241.96060,244.16670,298.11830,514.01560,296.35500,
    26.95794,105.53480,122.76670,324.06890,337.63890,365.59990,370.68230,
    402.636293,427.72263,450.325022,470.840891,489.558373,2.610111,119.62990
};
const double AGA892DC::Ki[21]={
    0.4619255,0.4479153,0.4557489,0.5279209,0.5837490,0.3825868,0.4618263,
    0.3514916,0.4533894,0.4186954,0.6406937,0.6341423,0.6738577,0.6798307,
    0.7175118,0.7525189,0.7849550,0.8152731,0.8437826,0.3589888,0.4216551
};
const double AGA892DC::Gi[21]={
    0.0000000,0.0278150,0.1890650,0.0793000,0.1412390,0.3325000,0.0885000,
    0.0343690,0.0389530,0.0210000,0.2566920,0.2818350,0.3322670,0.3669110,
    0.2897310,0.3375420,0.3833810,0.4273540,0.4696590,0.0000000,0.0000000
};
const double AGA892DC::Qi3=0.69;
const double AGA892DC::Qi6=1.06775;
const double AGA892DC::Qi7=0.633276;
const double AGA892DC::Si6=1.5822;
const double AGA892DC::Si7=0.39;
const double AGA892DC::E1j[18]={//j:2-19
    0.971640,0.960644,1.000000,0.994635,0.708218,0.931484,1.170520,0.990126,1.000000,
    1.019530,0.989844,1.002350,0.999268,1.107274,0.880880,0.880973,0.881067,0.881161
};
const double AGA892DC::U1j[18]={//j:2-19
    0.8861060,0.9638270,1.0000000,0.9908770,1.0000000,0.7368330,1.1563900,1.0000000,1.0000000,
    1.0000000,0.9922910,1.0000000,1.0036700,1.3025760,1.1919040,1.2057690,1.2196340,1.2334980
};
const double AGA892DC::K1j[18]={//j:2-19
    1.003630,0.995933,1.000000,1.007619,1.000000,1.000080,1.023260,1.000000,1.000000,
    1.000000,0.997596,1.000000,1.002529,0.982962,0.983565,0.982707,0.981849,0.980991
};
const double AGA892DC::G13=0.807653;
const double AGA892DC::G18=1.95731;
const double AGA892DC::E2j[12]={//j:3-14
    1.022740,0.970120,0.945939,0.746954,0.902271,1.086320,
    1.005710,1.021000,0.946914,0.973384,0.959340,0.945520
};
const double AGA892DC::U2j[10]={//j:3-12
    0.8350580,0.8164310,0.9155020,1.0000000,0.9934760,
    0.4088380,1.0000000,1.0000000,1.0000000,0.9935560
};
const double AGA892DC::K2j[6]={//j:3-8
    0.9823610,1.0079600,1.0000000,1.0000000,0.9425960,1.0322700
};
const double AGA892DC::G23=0.982746;
const double AGA892DC::E3j[16]={//j:4-19
    0.9250530,0.9602370,0.8494080,0.9550520,1.2817900,1.5000000,1.0000000,0.9068490,
    0.8973620,0.7262550,0.8597640,0.8551340,0.8312290,0.8083100,0.7863230,0.7651710
};
const double AGA892DC::U3j[16]={//j:4-19
    0.969870,1.000000,1.000000,1.045290,1.000000,0.900000,1.000000,1.000000,
    1.000000,1.000000,1.000000,1.066638,1.077634,1.088178,1.098291,1.108021
};
const double AGA892DC::K3j[16]={//j:4-19
    1.008510,1.000000,1.000000,1.007790,1.000000,1.000000,1.000000,1.000000,
    1.000000,1.000000,1.000000,0.910183,0.895362,0.881152,0.867520,0.854406
};
const double AGA892DC::G3j[3]={//j:4-6
    0.3702960,1.0,1.6730900
};
const double AGA892DC::E4j[10]={//j:5-14
    1.022560,0.693168,0.946871,1.164460,1.000000,
    1.000000,1.000000,1.013060,1.000000,1.005320
};
const double AGA892DC::U4j[10]={//j:5-14
    1.0651730,1.0000000,0.9719260,1.6166600,1.0000000,
    1.0000000,1.2500000,1.2500000,1.2500000,1.2500000
};
const double AGA892DC::K4j[4]={//j:5-8
    0.9868930,1.0000000,0.9999690,1.0203400
};
const double AGA892DC::E58=1.034787;
const double AGA892DC::E512=1.0049;
const double AGA892DC::E7j[5]={//j:15-19
    1.008692,1.010126,1.011501,1.012821,1.014089
};
const double AGA892DC::U7j[5]={//j:15-19
    1.028973,1.033754,1.038338,1.042735,1.046966

};
const double AGA892DC::K7j[5]={//j:15-19
    0.968130,0.962870,0.957828,0.952441,0.948338
};
const double AGA892DC::E89=1.1;
const double AGA892DC::E811=1.3;
const double AGA892DC::E812=1.3;
int AGA892DC::status=0;
//end const static member datas
//static member function:
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
    if(NUM_x>21 || NUM_x<1) return -22;
    double maxXi[21]={
        1.00000,0.70000,0.20000,0.10000,0.03500,0.00015,0.00015,0.10000,0.03000,0.00015,0.01500,
        0.01500,0.00500,0.00500,0.00100,0.00050,0.00050,0.00050,0.00050,0.00500,0.00015
    };
    if(p<=0 || p>12) return -100;//压力超出范围
    if((t+273.15)<263 || (t+273.15)>338) return -200;//温度超出范围
    double sum_0=0.0;
    for(int i=0;i<NUM_x;i++){
        int ii=ID_raw[i];
        if(ii<1 || ii>21) return -22;
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
        M+=Mi[ii-1]*xi_raw[Index2raw[i]];
    }
    rho=rhom*M;
    status=100; 
    return status;//calculation success
}

