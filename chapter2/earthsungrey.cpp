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
#define sqr(x) (x*x)

const double Ce= 2., Cs=2.0e-6 ;
const double rho0=1.,drho=-0.7; // density gradient
const double TOA=1.2; // Troposphere
const double Z=TOA*(1+drho*TOA/2);
const double B0=1.4744e-8, T0=4798;
const double pi = atan(1.)*4, stefan=sqr(sqr(pi))/15;

const double Te=(273+18)/T0, Ts=5798/T0;
double kappanu=0.5,  q0=-0.3;
const double mus=0.5;

const int Nz=60; // nb points in z
const int kmax=14;  // nb fixed point iterations
const double dz=Z/(Nz-1);
double  J0[Nz], S[Nz], T[Nz];

string basedir("");//("/Users/pironneau/Dropbox/aranger/TeX2026/BookVRTE/prog2/earthsungrey/");

double expint_E1(const double t=1){ // will not work if if t>2.5
    const double tmin=1.e-10; // min t in ExpInt(t)
    double t1=fabs(t);
    const int Kexpint = 9+(t1-1)*4;; // precision in exponential integral function E1
    const double  gaNtaua =0.577215664901533; // special integration for log(t)
    if(t1<tmin) return 0.;// because integral_0^t(logx)dx ~0
    if(t1>4) {
        cout << "value of E_1 is incorrect with "<<t<<">4"<<endl; return 0;}
    double ak=t1, soNtaue=-gaNtaua - log(t1)+ak;
    for(int k=2;k<Kexpint;k++){
        ak *= -t1*(k-1)/sqr(k);
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
    return (drho==0)?z: (sqrt(1+2*drho*z/rho0)-1)/drho;
}

int main(int argc, const char * argv[]) {
    const double I0 = Ce*stefan*sqr(sqr(Te)),
                 IZ = Cs*stefan*sqr(sqr(Ts));
        
    for(int i=0;i<Nz;i++)  J0[i] = 0;   // initialization
    
    for(int k=0;k<kmax;k++){ // iterations on the source
        double c0= -q0*mus*IZ*exp(-kappanu*Z/mus);
        for(int i=0;i<Nz;i++){
            c0 -= q0*expint_E2(i*dz*kappanu)*S[i]*dz;
            S[i] = kappanu*J0[i];
        }
        for(int i=0;i<Nz;i++) {
            double z=i*dz,
            J0z =  I0*expint_E3(z*kappanu)/2 + IZ*exp(-(Z-z)*kappanu/mus)/2 +
                    c0*expint_E2(z*kappanu)/2;
            for(int j=1;j<Nz;j++){
                double zp=j*dz;
                J0z += (S[j]+S[j-1])*expint_E1(kappanu*(fabs(z-zp+dz/2)))*dz/4;
            }
            J0[i] = J0z;
            T[i] = sqrt(sqrt(J0z/stefan));
         }
        cout<<k<<" "<< T[2]*T0-273<<endl;
    }
    cout<<endl;
    for(int i=0;i<Nz;i+=1){
        cout<< backtoz(i*dz)<<" "<<T[i]*T0-273<<endl;
    }
//    ofstream  resultfile(basedir+"T"+to_string(Nz)+".txt");
    ofstream  resultfile(basedir+"Txx.txt");
    for(int i=0;i<Nz;i++) resultfile<< backtoz(i*dz)<<" "<<T[i]*T0-273<<" "<<  B0*J0[i]<<endl;
 //   getI();
    return 0;
}
