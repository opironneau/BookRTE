//  Discrete ordinates solution of the polarized VRTE of Chapter 4: (4.9),(4.10),(4.11)
//
//  Same equations, same data and same cells as MCIQ4.cpp and IQ4.cpp, so that the three
//  methods can be compared. It was written to verify the Monte Carlo, and it is also the
//  third method quoted in the memo.
//
//  As in MCIQ4.cpp the unknowns are I_l=(I+Q)/2 and I_r=(I-Q)/2 of (4.8); the source
//  being a polynomial in mu^2, it is enough to carry the three moments
//      M0l=int I_l dmu,  M0r=int I_r dmu,  M2l=int mu'^2 I_l dmu,
//  the uncollided sunlight being added to them analytically (it is a Dirac).
//  For each frequency:
//    - the angular variable is discretized by a Gauss-Legendre rule on (0,1) and its
//      mirror image on (-1,0);
//    - a downward sweep (mu<0) then an upward sweep (mu>0) of the cells, exact for a
//      source constant per cell (step characteristics):
//          Ibar = S/kappa + (I_in - S/kappa) f,   I_out = S/kappa + (I_in - S/kappa) e,
//          e = exp(-kappa h/|mu|),  f = (|mu|/(kappa h)) (1-e);
//    - between the two sweeps the Lambert condition (4.11) is applied, the downward
//      current on the ground giving c_0 of (4.13);
//    - the two sweeps are repeated until the moments stop moving (source iteration).
//  The temperature equation (4.10) is the outer loop, solved as in IQ4.cpp.
//
//  compile with  g++ -O3 -std=c++17 -o SNIQ4 SNIQ4.cpp
//  run with      ./SNIQ4            or   ./SNIQ4 n=64 iter=25 tol=1e-10

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <ctime>
using namespace std;
#define sqr(x) ((x)*(x))
const double pi = 4*atan(1.);

// ---------------------------------------------------------------- data of IQ4.cpp
const double Ce=2., Cs=2.0e-6;             // c'_e and c_s
const double rho0=1., drho=-0.7;           // density gradient
#define zprime(z) (z*rho0*(1+drho*z/2))
const double TOA=1.2, Z=zprime(TOA);       // optical thickness of the atmosphere
const double T0=4798;                      // scaling (B_0=1.4744e-8)
const double Te=(273+18)/T0, Ts=5778/T0;   // temperature of the ground and of the sun
const double q0=-0.3;                      // Lambert albedo of the ground
const double mus=0.5;                      // direction of the collimated sunlight
const double _beta=0.5;                    // Rayleigh part of the scattering
const int    Nz=50;                        // nodes in z, as in IQ4.cpp
const int    Nc=Nz-1;                      // cells
const double dz=Z/(Nz-1);
const int    jmaxmax=600, nmax=256;
const double z1=zprime(0.4), z2=zprime(0.7), z3=zprime(0.8);
const double nu1=0.5, nu2=1;               // range of the Rayleigh scattering
const double cloud=1.5;                    // strength of the cloud in (z1,z2)
const double kappamin=0.001;

const double kapparhoz(double nu,double z){
    return 1 + cloud*((z>z1)-(z>z2))*((nu>3/10.)*(nu<3/1.) + (nu<3/20.));
}
const double as(double z,double nu){       // scattering albedo a_s
    return 0.3
    + cloud*0.3*4*(z2-z)*(z-z1)*(z>z1)*(z<z2)/sqr(z1-z2)
    + 0.3*16*(nu<nu2)*(nu>nu1)*sqr((nu-nu1)*(nu-nu2)/sqr(nu1-nu2))*(z>z3);
}
const double BB(const double nu,const double T){ return sqr(nu)*nu/(exp(fmin(nu/T,100))-1);}
const double dBB(const double nu,const double T){
    double a=exp(nu/T);
    if(a>1e100) return sqr(nu*nu/T)*1e-30;
    return a*sqr(nu*nu/fmax(1.e-30,a-1)/T);
}
const double Bsun(const double nu){ return BB(nu,Ts);}
const double Bearth(const double nu){ return BB(nu,Te);}
double backtoz(const double z){ return (drho==0)? z : (sqrt(1+2*drho*z/rho0)-1)/drho;}

// ---------------------------------------------------------------- run time parameters
int    ndir   = 32;         // Gauss points on (0,1): 2*ndir directions
int    kmax   = 25;         // outer iterations on the temperature
double tol    = 1e-8;       // stopping criterion of the source iteration
int    itmax  = 500;        // max source iterations

// ---------------------------------------------------------------- unknowns
double nu[jmaxmax], kappanu[jmaxmax];
double kapc[jmaxmax][Nc], asc[jmaxmax][Nc];
double J0[jmaxmax][Nc], J2[jmaxmax][Nc], K0[jmaxmax][Nc], K2[jmaxmax][Nc];
double Tc[Nc], zc[Nc], c0nu[jmaxmax];
double mu[nmax], wg[nmax];                 // Gauss-Legendre on (0,1), sum of wg = 1
double m0l[Nc], m0r[Nc], m2l[Nc], m2r[Nc]; // moments, without the uncollided sunlight
int    jmax;
long   nsweep=0;                           // total number of source iterations

// Gauss-Legendre rule on (0,1) with weights of total mass 1
void gauleg(int n){
    for(int i=0;i<(n+1)/2;i++){
        double z=cos(pi*(i+0.75)/(n+0.5)), z1,pp;
        do{
            double p1=1,p2=0;
            for(int j=0;j<n;j++){ double p3=p2; p2=p1; p1=((2*j+1)*z*p2-j*p3)/(j+1);}
            pp=n*(z*p1-p2)/(z*z-1);
            z1=z; z=z1-p1/pp;
        } while(fabs(z-z1)>1e-14);
        double x=-z, w=2/((1-z*z)*pp*pp);  // node on (-1,1) and its weight
        mu[i]      =(1+x)/2;  wg[i]      =w/2;
        mu[n-1-i]  =(1-x)/2;  wg[n-1-i]  =w/2;
    }
}

//*******************************
int readkappa(string mykappafile){
//------------------------------
    ifstream kappafile(mykappafile);
    if(!kappafile) throw runtime_error("Cannot open file: "+mykappafile);
    cout<<"reading kappa(nu) from file "<<mykappafile<<endl;
    int j=-1; double kappaux,waveno;
    double ks[Nz];
    while((j++<=jmaxmax)&&(!kappafile.eof())){
        kappafile >> waveno >> kappaux;
        kappanu[j]=fmax(kappaux,kappamin);
        nu[j]=3/waveno;
        if((nu[j]>3./18)&&(nu[j]<3./14)) kappanu[j] *=1.2; //CO2
        ks[0]=0;                            // same trapezoidal rule as IQ4.cpp
        for(int i=1;i<Nz;i++)
            ks[i]=ks[i-1]+(kapparhoz(nu[j],i*dz)+kapparhoz(nu[j],(i-1)*dz))*dz/2;
        for(int i=0;i<Nc;i++){
            kapc[j][i]=kappanu[j]*(ks[i+1]-ks[i])/dz;
            asc[j][i] =as((i+0.5)*dz,nu[j]);
        }
    }
    kappafile.close();
    cout<<"Number of frequencies "<<j<<endl;
    return j;
}

// ---------------------------------------------------------------- one frequency
void snOneFrequency(int j){
    const double nuj=nu[j], Bs=Bsun(nuj), Be=Bearth(nuj);
    const double *kap=kapc[j], *asz=asc[j];
    const double b=3*_beta/8, u=(1-_beta)/4;
    double sigs[Nc],siga[Nc],B[Nc],beam[Nc],e[nmax][Nc],f[nmax][Nc],tau=0;
    for(int i=0;i<Nc;i++){
        sigs[i]=kap[i]*asz[i]; siga[i]=kap[i]*(1-asz[i]); B[i]=BB(nuj,Tc[i]);
    }
    for(int i=0;i<Nc;i++) tau+=kap[i]*dz;                  // total optical depth
    double tz=0;
    for(int i=0;i<Nc;i++){                                 // uncollided sunlight
        beam[i]=Cs*Bs*exp(-(tau-tz-kap[i]*dz/2)/mus);
        tz+=kap[i]*dz;
    }
    for(int k=0;k<ndir;k++)                                // the two sweep coefficients
        for(int i=0;i<Nc;i++){
            e[k][i]=exp(-kap[i]*dz/mu[k]);
            f[k][i]=mu[k]/(kap[i]*dz)*(1-e[k][i]);
        }
    const double beam0=mus*Cs*Bs*exp(-tau/mus);            // beam current on the ground

    double Il[Nc],Ir[Nc],I2l[Nc],I2r[Nc];
    for(int it=0;it<itmax;it++){
        nsweep++;
        for(int i=0;i<Nc;i++) Il[i]=Ir[i]=I2l[i]=I2r[i]=0;
        // ---- source: S_l depends on mu^2, S_r does not
        double M0l[Nc],M0r[Nc],M2l[Nc],Sr[Nc],Sa[Nc],Sb[Nc];
        for(int i=0;i<Nc;i++){
            M0l[i]=m0l[i]+beam[i]/2; M0r[i]=m0r[i]+beam[i]/2;
            M2l[i]=m2l[i]+sqr(mus)*beam[i]/2;
            Sa[i]=sigs[i]*(2*b*(M0l[i]-M2l[i])+u*(M0l[i]+M0r[i]))+siga[i]*B[i]/2;
            Sb[i]=sigs[i]*b*(M2l[i]+M0r[i]-2*(M0l[i]-M2l[i]));  // coefficient of mu^2
            Sr[i]=sigs[i]*(b*(M2l[i]+M0r[i])+u*(M0l[i]+M0r[i]))+siga[i]*B[i]/2;
        }
        // ---- downward sweep, mu<0, no diffuse light entering at z=Z
        double down[nmax][2], fdown=beam0;
        for(int k=0;k<ndir;k++){
            double sl=0,sr=0, il=0,ir=0, m2=sqr(mu[k]);
            for(int i=Nc-1;i>=0;i--){
                sl=(Sa[i]+Sb[i]*m2)/kap[i]; sr=Sr[i]/kap[i];
                double bl=sl+(il-sl)*f[k][i], br=sr+(ir-sr)*f[k][i];
                Il[i]+=wg[k]*bl; Ir[i]+=wg[k]*br;
                I2l[i]+=wg[k]*m2*bl; I2r[i]+=wg[k]*m2*br;
                il=sl+(il-sl)*e[k][i]; ir=sr+(ir-sr)*e[k][i];
            }
            down[k][0]=il; down[k][1]=ir;
            fdown+=wg[k]*mu[k]*(il+ir);
        }
        // ---- Lambert ground, then upward sweep
        const double c0=-q0*fdown;
        for(int k=0;k<ndir;k++){
            double m2=sqr(mu[k]), il=Ce*Be*mu[k]/2+c0/2, ir=il, sl,sr;
            for(int i=0;i<Nc;i++){
                sl=(Sa[i]+Sb[i]*m2)/kap[i]; sr=Sr[i]/kap[i];
                double bl=sl+(il-sl)*f[k][i], br=sr+(ir-sr)*f[k][i];
                Il[i]+=wg[k]*bl; Ir[i]+=wg[k]*br;
                I2l[i]+=wg[k]*m2*bl; I2r[i]+=wg[k]*m2*br;
                il=sl+(il-sl)*e[k][i]; ir=sr+(ir-sr)*e[k][i];
            }
        }
        // ---- convergence of the source iteration
        double err=0,nrm=1e-300;
        for(int i=0;i<Nc;i++){
            err=fmax(err,fabs(Il[i]+Ir[i]-m0l[i]-m0r[i]));
            nrm=fmax(nrm,fabs(Il[i]+Ir[i]));
        }
        for(int i=0;i<Nc;i++){ m0l[i]=Il[i]; m0r[i]=Ir[i]; m2l[i]=I2l[i]; m2r[i]=I2r[i];}
        c0nu[j]=c0;
        if(err<tol*nrm) break;
    }
    for(int i=0;i<Nc;i++){     // J_q = (1/2) int mu^q I dmu, K_q the same for Q
        J0[j][i]=(m0l[i]+m0r[i]+beam[i])/2;
        J2[j][i]=(m2l[i]+m2r[i]+sqr(mus)*beam[i])/2;
        K0[j][i]=(m0l[i]-m0r[i])/2;
        K2[j][i]=(m2l[i]-m2r[i])/2;
    }
}

// ---------------------------------------------------------------- temperature equation
double rootT(const double rhs,const double Tx,const int i){
    double s=-rhs;
    for(int j=1;j<jmax;j++)
        s += (1-as(zc[i],nu[j]))*kappanu[j]*BB((nu[j]+nu[j-1])/2,Tx)*(nu[j]-nu[j-1]);
    return s;
}
double rhsT(const int i){
    double r=0;
    for(int j=1;j<jmax;j++)
        r += (1-as(zc[i],nu[j]))*kappanu[j]*(J0[j][i]+J0[j-1][i])/2*(nu[j]-nu[j-1]);
    return r;
}
void genT(){
    for(int i=0;i<Nc;i++){
        double rhs=rhsT(i);
        if(rhs<=0){ Tc[i]=0; continue;}
        double Ta=Tc[i]>0?Tc[i]:Te/2, Tb=2*Ta; int c=0;
        while(rootT(rhs,Ta,i)>0 && c++<100) Ta/=2;
        while(rootT(rhs,Tb,i)<0 && c++<100) Tb*=2;
        while(Tb-Ta>1e-4 && c++<200){
            double Tm=(Ta+Tb)/2; if(rootT(rhs,Tm,i)>0) Tb=Tm; else Ta=Tm;
        }
        Tc[i]=(Ta+Tb)/2;
        for(int it=0;it<50;it++){
            double left=0,deriv=0;
            for(int j=1;j<jmax;j++){
                double dnu=nu[j]-nu[j-1], nm=(nu[j]+nu[j-1])/2,
                       ka=(1-as(zc[i],nu[j]))*kappanu[j];
                left+=ka*BB(nm,Tc[i])*dnu; deriv+=ka*dBB(nm,Tc[i])*dnu;
            }
            double res=rhs-left;
            if(fabs(deriv)>1e-10) Tc[i]+=res/deriv;
            if(fabs(res)<1e-10) break;
        }
    }
}

//*******************************
int main(int argc,const char* argv[]){
//------------------------------
    for(int i=1;i<argc;i++){
        string s(argv[i]); size_t q=s.find('='); if(q==string::npos) continue;
        string k=s.substr(0,q); double v=atof(s.c_str()+q+1);
        if     (k=="n")    ndir=(int)v;
        else if(k=="iter") kmax=(int)v;
        else if(k=="tol")  tol=v;
        else cerr<<"unknown parameter "<<k<<endl;
    }
    if(ndir>nmax) ndir=nmax;
    gauleg(ndir);
    for(int i=0;i<Nc;i++){ zc[i]=(i+0.5)*dz; Tc[i]=Te/2; }
    jmax=readkappa("_kappa.txt");
    cout<<"discrete ordinates: "<<2*ndir<<" directions, "<<Nc<<" cells, tol="<<tol<<endl;

    double t0=clock();
    cout<<"\n iter\t T(ground)\t T(top)\n";
    for(int k=0;k<kmax;k++){
        for(int j=0;j<jmax;j++){
            for(int i=0;i<Nc;i++) m0l[i]=m0r[i]=m2l[i]=m2r[i]=0;
            snOneFrequency(j);
        }
        genT();
        cout<<" "<<k<<"\t "<<Tc[0]*T0-273<<"\t "<<Tc[Nc-1]*T0-273<<endl;
    }
    double cpu=(clock()-t0)/CLOCKS_PER_SEC;
    printf(" Time CPU = %10.3f s   (%ld source iterations, %.1f per frequency and per"
           " temperature iteration)\n",cpu,nsweep,double(nsweep)/jmax/kmax);

    ofstream fT("SNIQ4_T.txt");
    fT<<"# z   T[Celsius]   int J0 dnu   int K0 dnu\n";
    cout<<"\n z\t\t T[C]\t\t int J0 dnu\t int K0 dnu\n";
    for(int i=0;i<Nc;i++){
        double sJ=0,sK=0;
        for(int j=1;j<jmax;j++){
            sJ+=(J0[j][i]+J0[j-1][i])/2*(nu[j]-nu[j-1]);
            sK+=(K0[j][i]+K0[j-1][i])/2*(nu[j]-nu[j-1]);
        }
        fT<<backtoz(zc[i])<<" "<<Tc[i]*T0-273<<" "<<sJ<<" "<<sK<<"\n";
        if(i%6==0) cout<<backtoz(zc[i])<<"\t "<<Tc[i]*T0-273<<"\t "<<sJ<<"\t "<<sK<<endl;
    }
    fT.close();
    cout<<"\n c_0 of (4.13) at nu="<<nu[jmax/2]<<" : "<<c0nu[jmax/2]<<endl;
    cout<<" file SNIQ4_T.txt written"<<endl;
    return 0;
}
