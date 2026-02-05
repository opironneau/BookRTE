//  /Users/pironneau/Dropbox/aranger/TeX2025/BookVRTE/prog1/Milne1/Milne1/milneProblem.cpp
//
//  Created by Olivier Pironneau on 31/05/2025.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
using namespace std;
#define sqr(x) ((x)*(x))

const int Nz=200; // nb points in z and tau in (0,1)
const int kmax=15;  // nb fixed point iterations
const double tau1=1; // tropopause
const double dtau=tau1/(Nz-1);
const double mus=1/sqrt(3.); // Sun inclination is 60°
const double Is= 1; // Sunlight intensity
double J00[Nz], J0[Nz], S[Nz];
bool verbose = false;

double expint_E1(const double t=1){
    double t1=fabs(t);
    const int Kexpint =14+(t1-1)*4;; // precision in exponential integral function E1
    const double  gaNtaua =0.577215664901533; // special integration for log(t)
    if(t1<1e-5) return 0;
    if(t1>4) {
        if(verbose) cout << "argument of E_1>2.5"<<endl;
        double tx=1./t1;
        return exp(-t1)*tx * (1 +(-1+(2+(-6+(24+(-120+720*tx)*tx)*tx)*tx)*tx)*tx );
    }
    double ak=t1, soNtaue=-gaNtaua - log(t1)+ak;
    for(int k=2;k<Kexpint;k++){
        ak *= -t1*(k-1)/sqr(k);
        soNtaue += ak;
    }
    return soNtaue;
}

int getI(){
    const double pi=4*atan(1.), dmu=0.031;
    const int b=1 ; // 1 for Chandrasekhar, 0 for Dirac
    const double lam=12,  w = sqrt(pi)/2/lam * ( erf((mus+1)*lam) - erf((mus-1)*lam) ) ;
    for(double mu=-1; mu<1;mu+=dmu){
        for(int j=0;j<Nz;j++){
            double Iss= (mu<0)*Is*exp(-sqr((mu+mus)*lam))*exp(-j*dtau/mus)/w, I=(1-b)*Iss;
            for(int i=0;i<Nz;i++)
                I +=( (i<j)*(mu<0) + (i>j)*(mu>0) )* dtau*exp(-fabs((i-j)*dtau/fabs(mu)))
                        * ( J0[i] + b*Is*exp(-i*dtau/mus)/2 )/fabs(mu);
            cout<<j*dtau<<"\t"<<mu<<"\t"<< I+ b*Iss<< endl;
        }
        cout<<endl;
    }
    return 0;
}

int main(int argc, const char * argv[]) {
 
    for(int i=0;i<Nz;i++){
        J00[i] =  exp(-i*dtau/mus)*Is/2;
        J0[i] = J00[i];
    }
    
    for(int k=0;k<kmax;k++){ // iterations on the source
        for(int i=0;i<Nz;i++)
            S[i] = J0[i] ;
        for(int i=0;i<Nz;i++) {
            J0[i] = J00[i];
            for(int j=0;j<Nz;j++)
                 J0[i] += dtau/2 * expint_E1((j-i+0.5)*dtau)*S[j];
        }
        cout<<k<<" "<< J0[Nz/2]<<endl;
    }
    cout<<endl;
    for(int i=0;i<Nz;i++) cout<< i*dtau<<" "<< J0[i]<<endl;
    
//    getI();
    return 0;
}
