
//  Solve VRTE for I and Q (see Chapter 4) when kappa and as
//  are functions of nu and z
//  _beta=0 results in a solver for I only
//  Created by Olivier Pironneau on 10/11/2025.
//
// kappa file needs be in the directory "basedir"
// compile with g++ -O3 -std=c++17 -DIQ4_HAVE_GSL=0 -o IQ4 IQ4.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
using namespace std;
#define sqr(x) ((x)*(x))
const double pi = atan(1.)*4, stefan=sqr(sqr(pi))/15;

const double Ce= 2., Cs=2.0e-6 ;
const double rho0=1.,drho=-0.7; // density gradient
#define zprime(z) (z*rho0*(1+drho*z/2))
const double TOA=1.2, Z = zprime(TOA); // optical thickness
const double B0=1.4744e-8, T0=4798; // scaling

const double Te=(273+18)/T0, Ts=5778/T0;// Earth andSun
const double q0=-0.3; // albedo
const double mus=0.5; // direction of collimated light

const double _beta=0.5;  // % of isotropic to Rayleigh scattering

const bool verbose=false; // more output
const int Nz=50; // nb points in z
const double dz=Z/(Nz-1);
const int kmax=9;  // nb of ISIF iterations
const int jmaxmax=600; // max of max nb of points for integration in nu range
const int newton = 50; // to compute T from the temperature equation
const double epsdycho=1e-4, epsnewton=1.e-10;  // precision for dychotomy & Newton
const double kappamin=0.001;  // if kappa - read - is too small, max it with kappamin
const int nt = 5; // min nb of integration step in anal formula
const double z1 = zprime(0.4), z2=zprime(0.7), z3=zprime(0.8);  // altitudes for cloud and Rayleigh
const double nu1=0.5,nu2=1; // range of Rayleigh scattering

const double cloud = 1.5; // strength of cloud in (z1,z2)

const double kapparhoz(double nu, double z){ //change of kappa not due to density
    return  1 + cloud*((z>z1)-(z>z2)) // altitude dependence
    *((nu>3/10.)*(nu<3/1.) + (nu<3/20.));//frequency dependence from water vapor
}// total kappa is density times this function times kappanu

double  kappanu[jmaxmax], nu[jmaxmax];      // kappanu and uneven discretization of [numin,numax]
double  J0[jmaxmax][Nz], J0old[jmaxmax][Nz],// mu integral of I_nu and
        J2[jmaxmax][Nz], J2old[jmaxmax][Nz], // mu integral of mu^2 I_nu
        K0[jmaxmax][Nz], K0old[jmaxmax][Nz], // mu integral of mu^2 I_nu
        K2[jmaxmax][Nz], K2old[jmaxmax][Nz], // mu integral of mu^2 I_nu
        T[Nz], kappaz[jmaxmax][Nz], kappasum[jmaxmax][Nz];// Temperature, kappa(z) and int_0^z(kappa(z))

string basedir("");
string mykappafile(basedir+"_kappa.txt"); // 1 - kappa
string myresulttemperature(basedir+"temperatureyyx");
string myresultmeanintensity(basedir+"imean0");
int jmax; // defined in readkappa()


const double as(double z, double nu){// scattering
    return 0.3                                                   // background scattering
    + cloud*0.3*4*(z2-z)*(z-z1)*(z>z1)*(z<z2)/sqr(z1-z2)      // cloud scattering
    + 0.3*16*(nu<nu2)*(nu>nu1)*sqr((nu-nu1)*(nu-nu2)/sqr(nu1-nu2))*(z>z3);  // Rayleigh
}
const double expint_E1(const double t=1){
     double t1=fabs(t);
    const int Kexpint = 12+(t1-1)*4;; // precision in exponential integral function E1
    const double  gaNtaua =0.577215664901533; // special integration for log(t)
    if(t1<1e-5) return 0;
    if(t1>12) {
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
const double expint_E2(const double t=1){
    double t1=fabs(t);
    return exp(-t1) - t1*expint_E1(t1);
}
const double expint_E3(const double t=1){
    double t1=fabs(t);
    return (exp(-t1) - t1*expint_E2(t1))/2;
}
const double expint_E4(const double t=1){
    double t1=fabs(t);
    return (exp(-t1) - t1*expint_E3(t1))/3;
}
const double expint_E5(const double t=1){
    double t1=fabs(t);
    return (exp(-t1) - t1*expint_E4(t1))/4;
}

const double Bsun(const double nu){ return sqr(nu)*nu/(exp(nu/Ts) -1);} // Boltzmann
const double Bearth(const double nu){ return  sqr(nu)*nu/(exp(nu/Te) -1);} // Boltzmann
const double BB(const double nu, const double T){  return sqr(nu)*nu/(exp(fmin(nu/T,100)) -1);} // Boltzmann
const double dBB(const double nu, const double T){
        double a = exp(nu/T);
        if(a>1e100) return sqr(nu*nu/T)*1e-30; // needed but could endanger the Newton iterations
        return a*sqr(nu*nu/fmax(1.e-30,a -1)/T);
}


//*******************************
int readkappa(string mykappafile){
//------------------------------
    const double kappamin=0.001;  // if kappa read is too small max it with kappamin
// on each line of kappafile: wavelength, absorption
    cout<<"reading kappa(nu) from file "<<mykappafile<<endl;
    ifstream kappafile(mykappafile);
    if (!kappafile) {
        throw std::runtime_error("Cannot open file: " + mykappafile);
    }
    int j=-1;
    double kappaux, waveno;
    while((j++<=jmaxmax)&&(!kappafile.eof())){
        kappafile >> waveno >>kappaux;
        kappanu[j]  = fmax(kappaux,kappamin);
        nu[j]=3/waveno;
        kappasum[j][0]=0;
        for(int it=0;it<Nz;it++){ // kappa  is kappatau(z)*kappanu(nu);
            double z = it*dz;
            kappaz[j][it] = kapparhoz(nu[j],z);  // density
            if(it>0)
                kappasum[j][it] = kappasum[j][it-1]
                + (kappaz[j][it]+kappaz[j][it-1])*dz/2; // integral of kappaz
        }
    }
    kappafile.close();
    cout<<"Number of frequencies "<<j<<endl;
    return j;
}


int updateJK(){
    // returns the convolution t-integral for the mean in mu of kappa*I_nu and mu^2 kappa*I_nu
    for(int jnu=0; jnu<jmax; jnu++){  // frequency loop
        double H[Nz],S[Nz];
        double kappanuj = kappanu[jnu], nuj=nu[jnu];
        double  cc=kappasum[jnu][Nz-1],
                c0= - q0*mus*Cs*Bsun(nuj)*exp(-kappanuj*kappasum[jnu][Nz-1]/mus);
        
        for(int i=0;i<Nz;i++){ // initial altitude loop;
            double z=i*dz, kappazj=kappaz[jnu][i]*kappanuj,
                    asz=as(z,nuj), sigs=kappazj*asz, siga=kappazj*(1-asz);
            H[i] = 9*_beta*sigs/8 * (J2old[jnu][i] -J0old[jnu][i]/3 - K0old[jnu][i]+K2old[jnu][i]);
            S[i] = sigs*J0old[jnu][i]+siga*BB(nuj,T[i]) -H[i]/3 ;
            c0 -= q0* (expint_E2(kappasum[jnu][i]*kappanuj)*S[i]
                       + expint_E4(kappasum[jnu][i]*kappanuj)*H[i])*dz;
        }
        
        for(int i=0;i<Nz;i++) { // main altitude loop
            double ksz=kappasum[jnu][i]*kappanuj;
            double J0z  = Ce*Bearth(nuj)*expint_E3(ksz)/2;
            double aux = Cs*Bsun(nuj)*exp(-(kappasum[jnu][Nz-1]-kappasum[jnu][i])*kappanuj/mus)/2;
            J0z   += aux;
            J0z+= c0*expint_E2(ksz)/2;
            double J2z  = Ce*Bearth(nuj)*expint_E5(ksz)/2 + sqr(mus)*aux + c0*expint_E4(ksz)/2;
            double K0z=0, K2z=0;
 
            for(int j=1;j<Nz;j++){ // convolution loop
                double Hj=(H[j]+H[j-1])/2;
                double Sj=(S[j]+S[j-1])/2;
                int ip=i, jp=j;
                if(i<j){ip=j; jp=i;}
                double aux = kappanuj*(kappasum[jnu][ip]-kappasum[jnu][jp] -dz/2);
                double expE1ij = expint_E1(aux);
                double expE3ij = expint_E3(aux);
                double expE5ij = expint_E5(aux);
                J0z += (expE1ij*Sj + expE3ij*Hj)*dz/2;
                J2z += (expE3ij*Sj + expE5ij*Hj)*dz/2;
                K0z -= (expE1ij-expE3ij)*Hj*dz/2;
                K2z -= (expE3ij-expE5ij)*Hj*dz/2;
             }
            J0[jnu][i]=J0z;
            J2[jnu][i]=J2z;
            K0[jnu][i]=K0z;
            K2[jnu][i]=K2z;
        }
    }
    return 0;
}

double root(const double rhs, const double T0,int i){
    double myeq=-rhs;
    for(int j=1; j<jmax;j++)
         myeq += (1-as(i*dz,nu[j]))*kappanu[j]*BB((nu[j]+nu[j-1])/2,T0)*(nu[j]-nu[j-1]);
    return myeq;
}
double rhsTeq(const int i){
    double rhs=0;
    for(int j=1; j<jmax;j++)
        rhs+=(1-as(i*dz,nu[j]))*kappanu[j]*(J0[j][i]+J0[j-1][i])/2*(nu[j]-nu[j-1]);
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
                 kappaux=(1-as(i*dz,nu[j]))*kappanu[j];
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
            J0old[j][i]= J2old[j][i]= K0old[j][i]= K2old[j][i]=0;
        }
    }

    for(int k=0;k<kmax; k++){  // fixed point loop: first update F
        updateJK();
        genT();
        for(int i=0;i<Nz;i++)
            for(int j=0;j<jmax;j++) {
                J0old[j][i]=J0[j][i];
                J2old[j][i]=J2[j][i];
                K0old[j][i]=K0[j][i];
                K2old[j][i]=K2[j][i];
            }
        double normG=0;
        for(int j=0; j<jmax;j++)
            for(int i=0;i<Nz;i++){ normG+=fabs(J0[j][i])*(nu[j]-nu[j-1])*dz; }
        cout << "k= "<<k <<" "<<T[2]*4798-273<<"  "<<normG<<endl;
    }
     return 0;
}

double backtoz(const double z){
    return  (drho==0)? z : (sqrt(1+2*drho*z/rho0)-1)/drho;
}

int main(int argc, const char * argv[]) {
    //----------------------------------
    jmax=readkappa(mykappafile);
        
 //  { int K=2;
   for(int K=0;K<2;K++){
        for(int j=0;j<jmax;j++)
            if (K==1){
                double c=1.2; if((nu[j]>3./18)&&(nu[j]<3./14)) //CO2
 //             double c=2; if((nu[j]>3./4.5)&&(nu[j]<3./3.5)) //NO2
 //             double c=2; if((nu[j]>3/8.5)&&(nu[j]<3/7.)) //CH4
                kappanu[j]=c*kappanu[j];//for CO2
            }
            else if (K==2) kappanu[j]=0.5;  // constant kappa
        
        ofstream  resultfile =
        ofstream(myresulttemperature+to_string(K)+".txt");
        cout<<"\n iterations \t [T] near earth [T] far ||G|| and ||S||\n";
        double t0 = clock();
        //-----------------------
        multiBlock(Te/2);
        //-----------------------
        printf( " Time CPU = %10.6f\n", (clock() - t0)/CLOCKS_PER_SEC);
        cout<<"\n tau\t [T]:"<<endl;
        for(int i=0;i<Nz-1;i++){
            for(int j=1;j<jmax;j++) K0[0][i] += K0[j][i]*(nu[j]-nu[j-1]);
                cout<< backtoz(i*dz)<<" "<<T[i]*T0-273<<" "<<K0[0][i]<<endl;
  //          resultfile<< backtoz(i*dz)<<" "<<T[i]*T0-273<<" "<<K0[0][i]<<endl;
        }
        resultfile.close();
    }
    return 0;
}
