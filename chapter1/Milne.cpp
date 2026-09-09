// Solves Milne's Problem
// compile with g++ -O3 -std=c++17  -o Milne MilneProblem.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
using namespace std;
#define sqr(x) ((x)*(x))

const int Nz=60, Nmu=100; // nb points in z and tau in (0,1)
const int kmax=15;  // nb fixed point iterations
const double tau1=2;
const double dtau=tau1/(Nz-1);
const double mus=1./sqrt(3.); // an inclination is 60°
const double Is= 1; // Sunlight intensity
double J0[Nz], S[Nz];

double expint_E1(const double t=1){
     double t1=fabs(t); // allow calls with negative argument
    const int Kexpint = 9+(t1-1)*4; // precision in exponential integral function E1
    const double  gaNtaua =0.577215664901533; // special integration for log(t)
    if(t1<1e-5) return 0; // E1(t), t small is always called with t*E1(t)
    if(t1>4) {
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
    const double pi=4*atan(1.), dmu=2./Nmu;
    const int b=1 ; // 1 for Chandrasekhar, 0 for Dirac
    const double lam=12,  w = sqrt(pi)/2/lam * ( erf(mus*lam) + erf((1-mus)*lam) ) ;
    cout<<" tau,mu --> I(tau,mu) \n";
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

double schwarzschild(const double tau){
    const double J00=Is/2;
    return J00*(1-(1+tau)/(1+tau1)) +J00*(1-tau/(1+tau1));
}
int main(int argc, const char * argv[]) {
    
    for(int i=0;i<Nz;i++) // initialization
        J0[i] =  exp(-i*dtau/mus)*Is/2;
    
    cout<<"k     J0[1]"<<endl;
    for(int k=0;k<kmax;k++){ // iterations on the source
        for(int i=0;i<Nz;i++)
            S[i] = J0[i] ;
        for(int i=0;i<Nz;i++) {
            J0[i] = exp(-i*dtau/mus)*Is/2;
            for(int j=1;j<Nz;j++)
                J0[i] += dtau/4 * expint_E1((j-i-0.5)*dtau)*(S[j]+S[j-1]);
        }
        cout<<k<<" "<<J0[1]<<endl;
    }
    cout<<endl<<"tau      J0[tau]   Schwarzschild[tau]"<<endl;
    for(int i=0;i<Nz;i++) cout<< i*dtau<<" "<< J0[i]<<" "<< schwarzschild(i*dtau) << endl;
 //   getI(); //uncomment to display I(tau,mu)
    return 0;
}
