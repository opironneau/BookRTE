
//  main.cpp : can handle scattering
//  solve Milne problem and
//  solve multi-group for the same Milne for kappa=kappa2 constant (integration in nu)
//  and solve multi-group problem with  kappa=kappa2-dknu*(nu_1<nu<nu_2)
//
//  Created by Olivier Pironneau on 10/05/2021.
//
// kappa file needs be in the same folder
// compile with g++ -O3 -std=c++17 -DIQ4_HAVE_GSL=0 -o earthsunnu earthsunnu.cpp
// run by ./earthsunnu

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <time.h>
using namespace std;
#define sqr(x) ((x)*(x))
const double Ce= 2, Cs=2e-6 ;
const double drho=-0.7; // density gradient
const double Zatmo=1.2, Z = Zatmo*(1+drho*Zatmo/2); // optical thickness
const double B0=1.4744e-8, T0=4798;
const double pi = atan(1.)*4, stefan=sqr(sqr(pi))/15;

const double Te=(273+18)/T0, Ts=5798/T0;
const double q0=-0.3;
const double mus=0.5;
const double z1 = Z*0.5, z2=Z*0.7, z3=Z*0.8;  // altitudes for scattering
const double nu1=0.5,nu2=1; // range of Rayleigh scattering
double lambda=1.5; // for clouds =1.5
double cscat = 0.; // Rayleigh and enhanced scattering in clouds 

const bool verbose=false;
const int Nz=60; // nb points in tau
const double dz=Z/(Nz-1);
const int kmax=15;  // nb fixed point iterations
const int jmaxmax=6000; // max of max nb of points for integration in nu range
const int newton = 50; // to compute the temperature from int k*Botlzmann=int k*Imean
const double epsdycho=1e-4, epsnewton=1.e-10;  // precision for dychotomy & Newton
const double kappamin=0.001;  // if kappa read is too small max it with kappamin

vector<double> nu, kappanu;      // uneven discretization of [numin,numax], sized to jmax
vector<vector<double>> J0, J0old;  // mean radiation, sized jmax x Nz
double T[Nz];// temperature

string basedir(""); // put your own directory if you wish 
string mykappafile(basedir+"_kappaDense.txt");
string myresulttemperature(basedir+"temperaturec");
string myresultmeanintensity(basedir+"imean0");
int jmax;


const double as(double z, double nu){// scattering
    return 0.3                                         // background scattering
       + cscat*(z2-z)*(z-z1)*(z>z1)*(z<z2)*4/sqr(z1-z2)// cloud scattering
        + cscat*(nu<nu2)*(nu>nu1)*sqr(4*(nu-nu1)*(nu-nu2)/sqr(nu2-nu1))*(z>z3);}  // Rayleigh

 
double cloud(double z){ return 1+lambda*((z>z1)-(z>z2));}
double scloud(double z, double zp){
 return zp-z + lambda*((zp-z1)*(zp>z1)-(zp-z2)*(zp>z2) - (z-z1)*(z>z1)+(z-z2)*(z>z2));
}

double expint_E1(const double t=1){
     double t1=fabs(t)+1e-10;
    const int Kexpint = 10; // precision in exponential integral function E1
    double ak=t1, soNtaue=-0.577215664901533 - log(t1)+ak;
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
double Bsun(const double nu){ return sqr(nu)*nu/(exp(nu/Ts) -1);} // Boltzmann
double Bearth(const double nu){ return  sqr(nu)*nu/(exp(nu/Te) -1);} // Boltzmann
double BB(const double nu, const double T){  return sqr(nu)*nu/(exp(fmin(nu/T,100)) -1);
} // Boltzmann
double dBB(const double nu, const double T){
        double a = exp(nu/T);
        if(a>1e100) return sqr(nu*nu/T)*1e-30; // needed but could endanger the Newton iterations
        return a*sqr(nu*nu/fmax(1.e-30,a -1)/T);
}

int readkappa(string mykappafile){
// on each line of kappafile: wavelength, transmitance
    cout<<"reading kappa(nu) from file "<<mykappafile<<endl;
    ifstream kappafile(mykappafile);
    if (!kappafile) {
        throw std::runtime_error("Cannot open file: " + mykappafile);
    }
    double wavel, kappaux;
    while((int)nu.size()<jmaxmax && (kappafile >> wavel >> kappaux)){
        kappanu.push_back(fmax(kappaux,kappamin));
        nu.push_back(3/wavel);
     }
    kappafile.close();
    int j=(int)nu.size();
    cout<<"Number of frequencies "<<j<<endl;
    double auxe=0, auxs=0;
    for(int k=1;k<j;k++){
        auxe+=BB(nu[k],Te)*(nu[k]-nu[k-1]);
        auxs+=BB(nu[k],Ts)*(nu[k]-nu[k-1]);
    }
    if(verbose) cout <<"check Stefan id for the sun   "<< auxs << " "<<stefan*sqr(sqr(Ts))<<endl
                     <<"check Stefan id for the earth "<< auxe << " "<<stefan*sqr(sqr(Te))<<endl;
    return j;
}

int updateJ(){
    // returns the convolution t-integral for the mean in mu of kappa*I_nu and mu^2 kappa*I_nu
    for(int jnu=0; jnu<jmax; jnu++){  // frequency loop
        double S[Nz];
        double kappanuj = kappanu[jnu], nuj=nu[jnu];
        double bearthnu = Ce*Bearth(nu[jnu]),
               bsunnu = Cs*Bsun(nu[jnu]);
        double c0=  -q0*bsunnu*mus*exp(-scloud(0,Z)*kappanuj/mus);
        
        for(int i=0;i<Nz;i++){ // initial altitude loop;
            double z=i*dz, kappanuz = kappanuj*cloud(z);
            double sigs=kappanuz*as(z,nuj), siga=kappanuz-sigs;
            S[i] = sigs*J0old[jnu][i]+siga*BB(nuj,T[i]);
            c0 -= q0*expint_E2(kappanuj*scloud(0,z))*S[i]*dz;
        }
        
        for(int i=0;i<Nz;i++) { // main altitude loop
            double z=i*dz;
            double  J0z1  = bearthnu*expint_E3(kappanuj*scloud(0,z))/2,
                    J0z2  = bsunnu*exp(-kappanuj*scloud(z,Z)/mus)/2,
                    J0z3 =  c0*expint_E2(kappanuj*scloud(0,z))/2,
                    J0z4 = 0;
            for(int j=0;j<Nz;j++){ // convolution loop
                double zp=j*dz, Sj4=S[j]*dz/2;//(S[j]+S[j-1])*dz/4;
                J0z4 += expint_E1(kappanuj*fabs(scloud(zp,z))+0.5*dz) * Sj4;
            }
            J0[jnu][i]=J0z1+J0z2+J0z3+J0z4;
        }
    }
    return 0;
}

double root(const double rhs, const double T0,int i){
    double myeq=-rhs;
    for(int j=1; j<jmax;j++)
         myeq += (1-as(i*Z/(Nz-1),nu[j]))*kappanu[j]*BB((nu[j]+nu[j-1])/2,T0)*(nu[j]-nu[j-1]);
    return myeq;
}

double rhsTeq(const int i){
    double rhs=0;
    for(int j=1; j<jmax;j++)
        rhs+=(1-as(i*dz,nu[j]))*kappanu[j]*J0[j][i]*(nu[j]-nu[j-1]); // second order rule does not better
    return rhs;
}

double getTbydycho(const double rhs, const int i,const double Tstart){
    double Taux, T0=Tstart, T1=2*T0;
    if(rhs==0) return 0;
    int counter=0;
    double myeq0 =root(rhs,T0,i), myeq1=root(rhs,T1,i);
    while (myeq0>0 && counter<100){
        T0=T0/2; counter++;
        myeq0=root(rhs,T0,i);
        if(verbose) cout<<T0<<" down "<<myeq0<<endl;
    }
    while (myeq1<0 && counter<100){
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
        double rhs=rhsTeq(i);
        T[i] = getTbydycho(rhs,i,T[i]);
        double presfunc = 1;
        int inewton =0;
        while(inewton++<newton && fabs(presfunc) > epsnewton){
            double T0 = T[i];
            double left=0,  deriv=0;
            for(int j=1; j<jmax;j++){
                double dnu=nu[j]-nu[j-1], nu11=(nu[j]+nu[j-1])/2,
                 kappaux=(1-as(i*dz,nu[j]))*(kappanu[j]+kappanu[j-1])/2;
                 left += kappaux*BB(nu11,T0)*dnu;
                deriv += kappaux*dBB(nu11,T0)*dnu;
            }
            presfunc = rhs-left;
            if(fabs(deriv)>1e-10) T[i] = T0+presfunc/deriv;
            if(verbose)
                cout<<i<<" "<<" T= "<<T[i]<<" residue= "<<presfunc<<" deriv= "<<deriv<<" Tstefan= "<< sqrt(sqrt(rhs*15*2))/3.1416<<endl; // the 2 is kappanuj
        }
        if(inewton>=newton)cout << "Newton precision doubtful "<<endl;
    }
    return 0;
}


int multiBlock(double initT)
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
        double normG=0;
        for(int j=0; j<jmax;j++)
            for(int i=0;i<Nz;i++){ normG+=sqr(J0[j][i])*dz; }
        cout << "k= "<<k <<" "<<T[2]*4798-273<<"  "<<normG/B0<<endl;
    }
     return 0;
}

double backtoz(const double z){
    return (sqrt(1+2*drho*z)-1)/drho;
}

int main(int argc, const char * argv[]) {
    //----------------------------------
    jmax=readkappa(mykappafile);
    J0.assign(jmax, vector<double>(Nz));
    J0old.assign(jmax, vector<double>(Nz));

    { int K=2;
//    for(int K=0;K<2;K++){
        for(int j=0;j<jmax;j++)
            if (K==1){
                if((nu[j]>3./18)&&(nu[j]<3./14))
                    kappanu[j]=1.2*kappanu[j];//Visible R
            }
            else if (K==2) kappanu[j]=0.5;  // constant kappa
        
        string resstream = myresulttemperature+"5"+to_string(K)+".txt";
        if(lambda==0) resstream =myresulttemperature+to_string(K)+".txt";
 //       string resstream = myresulttemperature+"n0"+to_string(K)+".txt";// night
   //     if(lambda==0) resstream =myresulttemperature+"n"+to_string(K)+".txt";
        ofstream  resultfile = ofstream(resstream);
        cout<<"\n iterations \t [T] near earth [T] far ||G|| and ||S||\n";
        double t0 = clock();
        //-----------------------
        multiBlock(Te/2);
        //-----------------------
        printf( " Time CPU = %10.6f\n", (clock() - t0)/CLOCKS_PER_SEC);
        cout<<"\n tau\t [T]:"<<endl;
        for(int i=0;i<Nz;i+=1){
            double z=i*dz;
            cout<< backtoz(z)<<" "<<T[i]*T0-273<<endl;
        }
        for(int i=1;i<Nz;i++)
            resultfile<< backtoz(i*dz)<<" "<<T[i]*T0-273<<endl;
        resultfile.close();
        double rb=Cs*stefan*sqr(sqr(Ts));
        for(int j=0; j<jmax;j++) rb -=J0[j][Nz-1]*(nu[j]-nu[j-1]);
            cout<<" Radiation Budget "<<rb/B0<<endl;
    }
    return 0;
}
// file names for results.  K==0: 1,10,1Z,   K==1  11,110,11Z,   K==2  21, 210, 21Z
//  K=0, no CO2,  K=1 with CO2,  K=2  kappa=0.5 no CO2
// SB==Up
// epsilon=0  K==0: 2,20,2Z,   K==1  12,120,12Z,   K==2  22, 220, 22SZ
// epsilon!=0 K==0: 1,10,1Z,   K==1  11,110,11Z,   K==2  21, 210, 21Z
// SB==Down
// epsilon=0  K==0: 102,1020,102Z,   K==1  112,1120,112Z,   K==2  122, 1220, 122SZ
// epsilon!=0 K==0: 101,1010,101Z,   K==1  111,1110,111Z,   K==2  121, 1210, 121Z
