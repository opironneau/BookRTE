//
//  Created by Olivier Pironneau on 08/06/2025.
//
// compile with g++ earthsungreyquadra.cpp
// run by ./a.out
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
using namespace std;

const double Ce= 2., Cs=0.2e-5 ;
const double drho=-0.7; // density gradient
const double Zatmo=1.2, Z = Zatmo*(1+drho*Zatmo/2); // optical thickness
const double B0=1.4744e-8, T0=4798;
const double pi = atan(1.)*4, stefan=pi*pi*pi*pi/15;

const double Te=(273+18)/T0, Ts=5800/T0;
double kappanu=0.5,  q0=-0.3;
const double mus=0.5;

const int Nz=60; // nb points in z

const int kmax=15;  // nb fixed point iterations
const double dz=Z/(Nz-1);
double J00Z[Nz], J0[Nz], S[Nz], T[Nz];

double expint_E1(const double t){ // will not work if if t>14
     double abst=fabs(t)+1e-15;
    const int Kexpint = 10; // precision in exponential integral function E1
    const double  gaNtaua =0.577215664901533; // special integration for log(t)
    double ak=abst, soNtaue=-gaNtaua - log(abst)+ak;
    for(int k=2;k<Kexpint;k++){
        ak *= -abst*(k-1)/(k*k);
        soNtaue += ak;
    }
    return soNtaue;
}
double expint_E2(const double t){
    double t1=fabs(t);
    return exp(-t1) - t1*expint_E1(t1);
}
double expint_E3(const double t){
    double t1=fabs(t);
    return (exp(-t1) - t1*expint_E2(t1))/2;
}
double expint_E1b(const double t){ // will not work if if t>14
     double abst=fabs(t)+1e-15;
    const int Kexpint = 10; // precision in exponential integral function E1
    const double  gaNtaua =0.577215664901533; // special integration for log(t)
    double ak=abst, soNtaue=-gaNtaua +ak; // no log term
    for(int k=2;k<Kexpint;k++){
        ak *= -abst*(k-1)/(k*k);
        soNtaue += ak;
    }
    return soNtaue;
}

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

int main(int argc, const char * argv[])
{
    const double I0 = Ce*stefan*pow(Te,4), IZ=Cs*stefan*pow(Ts,4);
    
    for(int quadra=0;quadra<2;quadra++){
        
        for(int i=0;i<Nz;i++){
            J00Z[i] = I0*expint_E3(i*dz*kappanu)/2 + IZ*exp(-(Z-i*dz)*kappanu/mus)/2;
            J0[i] = J00Z[i];
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
                double J0aux = J00Z[i] + Q0*expint_E2(z*kappanu)/2;
                J0z=J0aux;J0log=0;
                for(int j=1;j<Nz;j++){
                    double zp=j*dz, Sj2=(S[j]+S[j-1])/4;
                    if(quadra==0){
                        J0z += expint_E1((fabs(z-zp)+dz/2)*kappanu)*Sj2*dz;
                    }else if (quadra==1){
                        J0z += expint_E1b((fabs(z-zp)+dz/2)*kappanu)*Sj2*dz;
                        if(z!=zp )
                            J0log +=(zp-z)*(log(kappanu*fabs(zp-z))-1)*Sj2;
                        if(z!=zp-dz)
                            J0log -=(zp-z-dz)*(log(kappanu*fabs(zp-z-dz))-1)*Sj2;
                    }
                }
                J0[i]=J0z-J0log;
            }
            for(int i=0;i<Nz;i++){
                T[i] = sqrt(sqrt(J0[i]/stefan));
                S[i] = kappanu*J0[i];
            }
                       cout<<k<<" "<< T[1]*4798-273<<endl;
        }
        cout<<endl;
        for(int i=0;i<Nz;i+=4){
            cout<< backtoz(i*dz)<<"     "<<T[i]*T0-273<<"   "<< J0[i]/B0<<endl;
        }
        //   getI();
    }
    return 0;
}
