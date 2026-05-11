//  /BookVRTE/prog/earth/earthNoT/earth1.cpp
//
//  Created by Olivier Pironneau on 08/06/2025.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
using namespace std;
#define sqr(x) ((x)*(x))

const double Ce= 2., Cs=0.2e-5 ;
const double drho=-0.7; // density gradient
const double Zatmo=1.2, Z = Zatmo*(1+drho*Zatmo/2); // optical thickness
const double B0=1.4744e-8, T0=4798;
const double pi = atan(1.)*4, stefan=sqr(sqr(pi))/15;

const double Te=(273+18)/T0, Ts=5700/T0;
double kappanu=0.5,  q0=-0.3;
const double mus=0.5;

const int Nz=100; // nb points in z
const int quadra=0;  //0,1,2  quadrature precision order; degree 2 fails

const int kmax=15;  // nb fixed point iterations
const double dz=Z/(Nz-1);
double J00[Nz], J0Z[Nz], J0[Nz], S[Nz], T[Nz];

string basedir("/Users/pironneau/Dropbox/aranger/TeX2026/BookVRTE/prog2/earthandsun/");



double expint_E1(const double t=1){ // will not work if if t>2.5
     double abst=fabs(t);
    const int Kexpint = 9+(abst-1)*4;; // precision in exponential integral function E1
    const double  gaNtaua =0.577215664901533; // special integration for log(t)
    if(abst<1e-5) return 0;
    if(abst>2.5) {
        cout << "value of E_1 is incorrect with "<<t<<">2.5"<<endl; return 0;}
    double ak=abst, soNtaue=-gaNtaua - log(abst)+ak;
    for(int k=2;k<Kexpint;k++){
        ak *= -abst*(k-1)/sqr(k);
        soNtaue += ak;
    }
    return soNtaue;
}

double expint_E1b(const double t=1){ // same as expint_E1 without the log
     double abst=fabs(t);
    const int Kexpint = 19+(abst-1)*4;; // precision in exponential integral function E1
    const double  gaNtaua =0.577215664901533; // special integration for log(t)
    if(abst<1e-5) return 0;
    if(abst>2.5) {
        cout << "value of E_1 is incorrect with "<<t<<">2.5"<<endl; return 0;}
    double ak=abst, soNtaue=-gaNtaua +ak;
    for(int k=2;k<Kexpint;k++){
        ak *= -abst*(k-1)/sqr(k);
        soNtaue += ak;
    }
    return soNtaue;
}
double expint_E2(const double t=1){
    double t1=fabs(t);
    return exp(-t1) - t1*expint_E1(t1);
}
double expint_E3(const double t=1){
    double t1=fabs(t);
    return (exp(-t1) - t1*expint_E2(t1))/2;
}
double expint_E4(const double t=1){
    double t1=fabs(t);
    return (exp(-t1) - t1*expint_E3(t1))/3;
}
double expint_E5(const double t=1){
    double t1=fabs(t);
    return (exp(-t1) - t1*expint_E4(t1))/4;
}

//string mykappafile(basedir+"transmitance25.txt"); // 1 - kappa
string mykappafile(basedir+"_kappa.txt"); // 1 - kappa

double BB(const double nu, const double T)
{  return sqr(nu)*nu/(exp(fmin(nu/T,50)) -1);} // Boltzmann

int getI(){
    const double dmu=0.051;
    
    for(double mu=-1; mu<1;mu+=dmu){
        for(int j=0;j<Nz;j++){
            double sums=0;//(1-b)/2*exp(-j*(Nz-1)/mus)*delta; // can't handle a Dirac
            for(int i=0;i<Nz;i++)
                sums += ((i<j)*(mu<0)
                         + (i>j)*(mu>0))*dz*exp(-fabs((i-j)*dz/mu)) * S[i]/fabs(mu);
            cout<<j*dz<<" "<<mu<<" "<< sums << endl;
        }
        cout<<endl;
    }
    return 0;
}

double backtoz(const double z){
    return 1/drho * (sqrt(1+2*drho*z)-1);
}

int main(int argc, const char * argv[]) {
    const int kappaM=2,iE1max=kappaM*Nz; // to host kappanuMax*Z
    double E1[iE1max];
    
    for(int i=0;i<iE1max;i++) E1[i] = expint_E1(i*dz);
    
    const double I0 = Ce*stefan*sqr(sqr(Te)), IZ=Cs*stefan*sqr(sqr(Ts));
        
    for(int i=0;i<Nz;i++){
        J00[i] = I0*expint_E3(i*dz*kappanu)/2;
        J0Z[i] = IZ*exp(-(Z-i*dz)*kappanu/mus)/2;
        J0[i] = J00[i]+J0Z[i];
        S[i] = kappanu*J0[i] ;   // initialization
    }
     for(int k=0;k<kmax;k++){ // iterations on the source
        double Q0= - q0*kappanu*mus*IZ*exp(-kappanu*Z/mus);
        for(int i=0;i<Nz;i++){
            Q0 -= q0*kappanu* expint_E2(i*dz*kappanu)*J0[i]*dz;
        }
        double J0log=0, J0z=0;
        for(int i=0;i<Nz;i++) {
            double z=i*dz;
            double J0aux = J00[i] + J0Z[i] + Q0*expint_E2(z*kappanu)/2;
            J0z=J0aux;J0log=0;
            for(int j=1;j<Nz;j++){
                double zp=j*dz, Sj=(S[j]+S[j-1])/4;
                if(quadra==0){
                    J0z += expint_E1((z-zp+dz/2)*kappanu)*Sj*dz;
 /*                   if(z<zp-dz) // not better
                        J0z -= (expint_E2(kappanu*(zp-z))
                                -expint_E2(kappanu*(zp-z-dz)))*Sj/kappanu;
                     else if(z>zp)
                        J0z += (expint_E2(kappanu*(zp-z))
                                -expint_E2(kappanu*(zp-z-dz)))*Sj/kappanu;
                    else
                        J0z += (2*expint_E2(0.)-expint_E2(kappanu*(zp-z))
                                -expint_E2(kappanu*(zp-z-dz)))*Sj/kappanu;
  */
                }else if (quadra==1){
                    J0z += expint_E1b((z-zp+dz/2)*kappanu)*Sj*dz;
                    if(z!=zp)
                        J0log +=(zp-z)*(log(kappanu*fabs(zp-z))-1)*Sj;
                    if(z!=zp-dz)
                        J0log -=(zp-z-dz)*(log(kappanu*fabs(zp-z-dz))-1)*Sj;
                } else if (quadra==2){  // bug??
                    if(z!=zp)
                        J0log +=(zp-z)*log(kappanu*fabs(zp-z))*(S[j]+(zp-z)/dz*S[j-1])/2
                        -S[j-1]/dz/4*sqr(zp-z)*(log(kappanu*fabs(zp-z))-0.5);
                    if(z!=zp-dz)
            J0log +=-((zp-z-dz)*log(kappanu*fabs(zp-z-dz))+dz)*(S[j]+(zp-z)/dz*S[j-1])/2
                        +S[j-1]/dz/4*sqr(zp-dz-z)*(log(kappanu*fabs(zp-dz-z))-0.5);
               }
                
            }
            J0[i]=J0z;
            if(quadra>0) {
                J0[i] -= J0log;
            }
        }
        for(int i=0;i<Nz;i++){
            T[i] = sqrt(sqrt(J0[i]/stefan));
            S[i] = kappanu*J0[i];
        }
            cout<<k<<" "<< T[1]*4798-273<<endl;
        }
    cout<<endl;
    for(int i=0;i<Nz;i+=1){
        cout<< backtoz(i*dz)<<" "<< B0*J0[i]<<" "<<T[i]*T0-273<<endl;
    }
    ofstream  resultfile(basedir+"J0"+to_string(Nz)+to_string(quadra)+".txt");
    for(int i=0;i<Nz-1;i++)
        resultfile<< backtoz(i*dz)<<" "<<  B0*J0[i]<<" "<<T[i]*T0-273<<endl;
 //   getI();
    return 0;
}
