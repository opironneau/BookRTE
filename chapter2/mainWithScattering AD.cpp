
//  Solve RTE with scattering and z-dependent kappa
//  Created by Olivier Pironneau on 10/05/2025.
// kappa is fixed. For variable add "readkappa" function
//  and put kappa file  next to basedir
// file names for results: K=0, no CO2,  K=1 with CO2

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
#include "ddouble.hpp"
using namespace std;
#define sqr(x) ((x)*(x))
const double Ce= 2.5, Cs=0.2e-5 ; // Earth and Sun power
const double rho0=1.,drho=-0.7; // density and  gradient
const double cloud = 1.; // augmented kappa in cloud
const double z1 = 0.4, z2=0.7, z3=0.8;  // altitudes for scattering
const double Zatmo=1.2, Z = rho0*(Zatmo*(1+drho*Zatmo/2)
+(1-cloud)*(z1-z2)*(1+drho*(z1+z2)/2)); // optical thickness
const double B0=1.4744e-8, T0=4798; // scaling for I and T
const double pi = atan(1.)*4, stefan=sqr(sqr(pi))/15;

const double Te=(273+18.)/T0, Ts=5778/T0;
const double q0=-0.3;
const double mus=0.5;

const bool verbose=false; // use it to debug
const int Nz=80; // nb points in tau
const double dz=Z/(Nz-1);
const int kmax=9;  // nb fixed point iterations
const int jmaxmax=600; // max of max nb of points for integration in nu range
const int newton = 50; // to compute the temperature from int k*Botlzmann=int k*Imean
const double epsdycho=1e-4, epsnewton=1.e-10;  // precision for dychotomy & Newton
const double kappamin=0.001;  // if kappa read is too small max it with kappamin
const int nt = 5; // min nb of integration step in anal formula
//const double z1 = Z*0.4, z2=Z*0.8, z3=Z*0.8;  // altitudes for scattering
//const double z1 = Z, z2=Z, z3=Z;  // comment line above for no cloud
const double nu1=0.5,nu2=1; // range of Rayleigh scattering

double  nu[jmaxmax];
ddouble kappanu[jmaxmax];      // uneven discretization of [numin,numax]
ddouble  J0[jmaxmax][Nz], J0old[jmaxmax][Nz],T[Nz];// mean radiation and temperature

string basedir("/Users/pironneau/Dropbox/aranger/TeX2026/BookVRTE/prog2/greenhouse4/");
string mykappafile(basedir+"_kappa.txt"); // 1 - kappa
string myresulttemperature(basedir+"temperature2xxx");
string myresultmeanintensity(basedir+"imean0");
int jmax=500;

const int ntablev=10;
const double kappatrunc=2.1;  // needed for the tabulations



const double as(double z, double nu){ return 0.3*(1+0.5*(z2-z)*(z-z1)*(z>z1)*(z<z2)*4/sqr(z1-z2)
                           + 0.5*(nu<nu2)*(nu>nu1)*sqr(sqr(nu/nu2))*(z>z3));} // scattering coefficient


ddouble expint_E1(const ddouble t=1){ // will not work if if t>2.5
     ddouble abst=abs(t);
    const int Kexpint = 9+abst.val[0]*4;; // precision in exponential integral function E1
    const ddouble  gaNtaua =0.577215664901533; // special integration for log(t)
    if(abst<1e-5) return 0;
    if(abst>2.5) {
        cout << "value of E_1 is incorrect with "<<t<<">2.5"<<endl; return 0;}
    ddouble ak=abst, soNtaue=-gaNtaua - log(abst)+ak;
    for(int k=2;k<Kexpint;k++){
        ak *= -abst*(k-1)/sqr(k);
        soNtaue += ak;
    }
    return soNtaue;
}
ddouble expint_E2(const ddouble t=1){
    ddouble t1=abs(t);
    return exp(-t1) - t1*expint_E1(t1);
}
ddouble expint_E3(const ddouble t=1){
    ddouble t1=abs(t);
    return (exp(-t1) - t1*expint_E2(t1))/2;
}
double Bsun(const double nu){ return sqr(nu)*nu/(exp(nu/Ts) -1);} // Boltzmann
double Bearth(const double nu){ return  sqr(nu)*nu/(exp(nu/Te) -1);} // Boltzmann
ddouble BB(const double nu, const ddouble T){  return sqr(nu)*nu/(exp(nu/T) -1);} // Boltzmann
ddouble dBB(const double nu, const ddouble T){
        ddouble a = exp(nu/T);
        if(a.val[0]>1e100) return sqr(nu*nu/T)*1e-30; // needed but could endanger the Newton iterations
        return a*sqr(nu*nu/(a -1)/T);
}


int updateJ(){
    // returns the convolution t-integral for the mean in mu of kappa*I_nu and mu^2 kappa*I_nu
    for(int jnu=0; jnu<jmax; jnu++){  // frequency loop
        ddouble S[Nz];
        ddouble kappanuj = kappanu[jnu];
        double nuj=nu[jnu];
        ddouble Q0= - q0*kappanuj*mus*Cs*Bsun(nuj)*exp(-kappanuj*Z/mus);
        
        for(int i=0;i<Nz;i++){ // initial altitude loop;
            double z=i*dz;
            Q0 -= q0*kappanuj* expint_E2(i*dz*kappanuj)*J0old[jnu][i]*dz;
            S[i] = kappanuj*(as(z,nuj)*J0old[jnu][i]+(1-as(z,nuj))*BB(nuj,T[i])) ;
        }
        
        for(int i=0;i<Nz;i++) { // main altitude loop
            double z=i*dz;
            ddouble J0z  = Ce*Bearth(nuj)*expint_E3(i*dz*kappanuj)/2
                        + Cs*Bsun(nuj)*exp(-(Z-z)*kappanuj/mus)/2
                        + Q0*expint_E2(z*kappanuj)/2;
            for(int j=1;j<Nz;j++){ // convolution loop
                ddouble expE1 = expint_E1(kappanuj*(i-j+0.5)*dz);
                J0z += expE1*(S[j]+S[j-1])*dz/4;
            }
            J0[jnu][i]=J0z;
        }
    }
    return 0;
}

ddouble root(const ddouble rhs, const ddouble T0,int i){
    ddouble myeq=-rhs;
    for(int j=1; j<jmax;j++)
         myeq += (1-as(i*Z/(Nz-1),nu[j]))*kappanu[j]*BB((nu[j]+nu[j-1])/2,T0)*(nu[j]-nu[j-1]);
    return myeq;
}

ddouble rhsTeq(const int i){
    ddouble rhs=0;
    for(int j=1; j<jmax;j++)
        rhs+=(1-as(i*dz,nu[j]))*kappanu[j]*J0[j][i]*(nu[j]-nu[j-1]);
    return rhs;
}

ddouble getTbydycho(const ddouble rhs, const int i,const ddouble Tstart){
    ddouble Taux, T0=Tstart, T1=2*T0;
    if(rhs==0) return 0;
    int counter=0;
    ddouble myeq0 =root(rhs,T0,i), myeq1=root(rhs,T1,i);
    while (myeq0.val[0]>0 && counter<100){
        T0=T0/2; counter++;
        myeq0=root(rhs,T0,i);
        if(verbose) cout<<T0<<" down "<<myeq0<<endl;
    }
    while (myeq1.val[0]<0 && counter<100){
        T1=2*T1; myeq1=root(rhs,T1,i); counter++;
        if(verbose) cout<<T1<<" up "<<myeq1<<endl;
    }
    if(verbose){
        myeq0=root(rhs,T0,i);
        myeq1=root(rhs,T1,i);
        if(myeq0>myeq1) cout<<"BUG in dychotomy\n";
    }
    while (T1-T0 > epsdycho && counter<100){
        Taux=(T1+T0)/2; counter++;
        myeq0=root(rhs,Taux,i);
        if(myeq0>0) T1=Taux; else T0=Taux;
        if(verbose){
            myeq0=root(rhs,T0,i);
            myeq1=root(rhs,T1,i);
            cout<<T0<<" middle "<<T1<<" "<<myeq0 << " "<< myeq1<<endl;
        }
    }
    if(counter>99) cout<<"Divergence in dichotomy\n";
    return (T1+T0)/2;
}

int genT(){
    for(int i=0;i<Nz;i++){
        ddouble rhs=rhsTeq(i);
        T[i] = getTbydycho(rhs,i,T[i]);
        ddouble presfunc = 1;
        int inewton =0;
        while(inewton++<newton && abs(presfunc) > epsnewton){
            ddouble T0 = T[i];
            ddouble left=0,  deriv=0;
            for(int j=1; j<jmax;j++){
                double dnu=nu[j]-nu[j-1], nu11=(nu[j]+nu[j-1])/2;
                ddouble kappaux=(1-as(i*Z/Nz,nu[j]))*(kappanu[j]+kappanu[j-1])/2;
                 left = left + kappaux*BB(nu11,T0)*dnu;
                deriv = deriv + kappaux*dBB(nu11,T0)*dnu;
            }
            presfunc = rhs-left;
            if(fabs(deriv.val[0])>1e-10) T[i] = T0+presfunc/deriv;
            if(verbose)
                cout<<i<<" "<<" T= "<<T[i]<<" residue= "<<presfunc<<" deriv= "<<deriv<<" Tstefan= "<< sqrt(sqrt(rhs*15*2))/3.1416<<endl; // the 2 is kappanuj
        }
        if(inewton>=newton)cout << "Newton precision doubtful "<<endl;
    }
    return 0;
}


int multiBlock(ddouble initT)
{
    for(int i=0;i<Nz;i++){
        T[i]=initT;  // initialize
        for(int j=0;j<jmax;j++) {
            J0old[j][i]=0;
         }
        }

    for(int k=0;k<kmax; k++){  // fixed point loop: first update F
        updateJ();
        genT();
        for(int i=0;i<Nz;i++)
            for(int j=0;j<jmax;j++) {
                J0old[j][i]=J0[j][i];
             }
        
        cout << "k= "<<k <<" "<<T[2]*4798-273<<endl;
    }
     return 0;
}

double backtoz(const double z){
    return rho0/drho * (sqrt(1+2*drho*z/sqr(rho0))-1);
}

int main(int argc, const char * argv[]) {
    //----------------------------------
    const int lambda=7;
    const double numin=3./lambda, numax=numin+0.02; //defines the region of differentiation
    const double kappa0=0.25;
    
    for(int j=0;j<jmax;j++){
        kappanu[j]=ddouble(kappa0,0.);  // constant kappa
        nu[j]=0.01+j*10./(jmax-1);
        
        if(nu[j]>numin && nu[j]<numax) kappanu[j]=ddouble(kappa0,1.);
    }
        ofstream  resultfile = ofstream(myresulttemperature+to_string(lambda)+".txt");
    
        cout<<"\n iterations \t [T] near earth [T] far ||G|| and ||S||\n";
        double t0 = clock();
        //-----------------------
        multiBlock(Te/2);
        //-----------------------
        printf( " Time CPU = %10.6f\n", (clock() - t0)/CLOCKS_PER_SEC);
        cout<<"\n tau\t [T]:"<<endl;
        for(int i=0;i<Nz;i+=1){
            cout<< backtoz(i*dz)<<" "<<T[i]*T0-273<<endl;
        }
        for(int i=2;i<Nz;i++)
            resultfile<< backtoz(i*dz)<<" "<<T[i].val[0]*T0-273<<" "<<(T[i].val[0]+30*T[i].val[1])*T0-273<<endl;
        resultfile.close();
    return 0;
}

