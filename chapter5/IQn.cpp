
//  Created by Olivier Pironneau on 01/03/2026.
//

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
using namespace std;
#define sqr(x) ((x)*(x))
const double pi = atan(1.)*4;

const double Ce= 2.4, Cs=00*2e-6 ;
const double rho0=1.,drho=-0.7; // density gradient
const double Z=1.2;
const double B0=1.4744e-8, T0=4798;

const double Te=(273+18)/T0, Ts=5778/T0;
const double q0=-0.3;
const double mus=0.5;
const double beta=0.5;

const bool verbose=false;
const int Nmu=251, Nz=80; // nb points in mu and in z
const double dmu=2./(Nmu-1), dz=Z/(Nz-1);
const int kmax=10;  // nb fixed point iterations
const int jmaxmax=600; // max of max nb of lines in file _kappa.txt
const int newton = 50; // max iters to solve inte k*Botlzmann=inte k*Imean
const double epsdycho=1e-4, epsnewton=1.e-10;  // precision for dychotomy & Newton

const bool cloud=false;
const double z1 = 0.8, z2=0.9, z3=0.8;  // altitudes for cloud and Rayleigh scattering
const double nu1=0.5,nu2=1; // range of Rayleigh scattering

const double Y=0.5, nm=1., np=nm-0.3; // must match jump of n(). See also kapparho
bool approx=1; // use approximation to compute jump contribution to
const double n(double z) {  return  nm + (np-nm)*(z>Y);}  // diffraction index

const double kapparhoz(double z){ return (1+drho*z)*(1*(z<Y) + 1*(z>=Y));}// density

struct doublet{ double d0=0,d1=0;};

double  kappanu[jmaxmax], nu[jmaxmax];      // kappanu and uneven discretization of [numin,numax]
double  J0[jmaxmax][Nz], J0old[jmaxmax][Nz],// mu integral of I_nu and
        J2[jmaxmax][Nz], J2old[jmaxmax][Nz], // mu integral of mu^2 I_nu
        K0[jmaxmax][Nz], K0old[jmaxmax][Nz], // mu integral of mu^2 I_nu
        K2[jmaxmax][Nz], K2old[jmaxmax][Nz], // mu integral of mu^2 I_nu
        T[Nz], nn[Nz],  // Delta for jump and Temeprature
        kappaz[Nz], kappasum[Nz+1];         // kappa(z) and int_0^z(kappa(z))

string basedir("/Users/pironneau/Dropbox/aranger/TeX2026/BookVRTE/prog5/");
string mykappafile(basedir+"_kappa.txt"); // 1 - kappa
string myresulttemperature(basedir+"temperatureff"); // use z with np<nm, w otherwise
string myresultmeanintensity(basedir+"imean0");
int jmax;


const double as(double z, double nu){// scattering
    return 0.3                                          // background scattering
    + 0.3*4*(z2-z)*(z-z1)*(z>z1)*(z<z2)/sqr(z1-z2)      // cloud scattering
    + 0.3*16*(nu<nu2)*(nu>nu1)*sqr((nu-nu1)*(nu-nu2)/sqr(nu1-nu2))*(z>z3);}  // Rayleigh


/********************************** END of DECLARATIONS ****************************************/

double expint_E1(const double t=1){
     double t1=fabs(t);
    const int Kexpint = 9+(t1-1)*4;; // precision in exponential integral function E1
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
    return (exp(-t1) - t*expint_E3(t1))/3;
}
double expint_E5(const double t=1){
    double t1=fabs(t);
    return (exp(-t1) - t*expint_E4(t1))/4;
}

double Bsun(const double nu){ return sqr(nu)*nu/(exp(nu/Ts) -1);} // Boltzmann
double Bearth(const double nu){ return  sqr(nu)*nu/(exp(nu/Te) -1);} // Boltzmann
double BB(const double nu, const double T){  return sqr(nu)*nu/(exp(fmin(nu/T,100)) -1);} // Boltzmann
double dBB(const double nu, const double T){
        double a = exp(nu/T);
        if(a>1e100) return sqr(nu*nu/T)*1e-30; // needed but could endanger the Newton iterations
        return a*sqr(nu*nu/fmax(1.e-30,a -1)/T);
}
double I0(const double mu, const double nu){ return Ce *mu* Bearth(nu);}
const double xlam=7;
const double lamsc=sqrt(pi)/2/xlam*(erf((1+mus)*xlam)+erf((1-mus)*xlam));
double IZ(const double mu, const double nu){ return Cs/lamsc * exp(-sqr((xlam*(mu+mus)))) * Bsun(nu);}

double kappamin=3, kappamax=0;  // updated in readKappa
//*******************************
int readkappa(string mykappafile){
//------------------------------
    kappasum[0]=0;
    for(int it=0;it<Nz;it++){ // kappa  is kappatau(z)*kappanu(nu);
        double z = it*dz;
        nn[it] = n(z);
        double aux= (cloud && z>z1 && z<z2) ? 2 : 1;
        kappaz[it] = aux*kapparhoz(z);  // density
        kappasum[it+1] = kappasum[it] + kappaz[it]*dz; // integral of kappaz
    }
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
        kappanu[j]  = fmin(3.,fmax(kappaux,0.01));
        nu[j]=3/waveno;
        if(kappaux>kappamax) kappamax=kappaux;
        if(kappaux<kappamin) kappamin=kappaux;
    }
    kappafile.close();
    cout<<"Number of frequencies "<<j
        << " kappamax = "<<kappamax<< " kappamin = "<<kappamin<<endl;
    return j;
}

// in Gemini kappamax is 2.5
const int Nkappamax=100; // kappanuj will be rounded to
//      kappanuj=kappamin+(kappamax-kappamin)*Nkappa/Nkappamax
double PHI[Nkappamax][Nz][Nz][Nmu]; // A trick not to compute phi twice
inline double phi(const double kappanuj,const double zp,const double zpp,const double mu)
{
    // compute exp(-kappanuj\int_zp^zpp kapparho(z")/mu(z")dz") with mu(z") as in mu2 below
    if(zpp-zp < 1e-6) return 1;
    if(zp>zpp || mu<=0){
        cout << "bug in calling phi\n";
        return 1;
    }
    int izp=zp/dz, izpp=zpp/dz, jmu=(Nmu-1)*(1+mu)/2,
    knuj=(kappanuj-kappamin)*Nkappamax/(kappamax-kappamin);
    if(PHI[knuj][izp][izpp][jmu]>=0)
        return PHI[knuj][izp][izpp][jmu]; // this saves computing exponentials
    else{
        int izp=zp/dz, izpp=zpp/dz;
        double aux = fabs((kappasum[izpp]-kappasum[izp])*kappanuj/mu);
        if(aux>100) return 0.;
        PHI[knuj][izp][izpp][jmu] = exp(-aux);
        return PHI[knuj][izp][izpp][jmu];
    }
}
// first index is a tabulation. If Xpm[Nmu][2][2]=-200, then getXY is called
double Xpm[Nmu][2][2], Ypm[Nmu][2][2], Xmp[Nmu][2][2],Ymp[Nmu][2][2];

int getXY(const int mppm, const double mu){
    // Compute Xpm,Ypm if mppm=1 else Xmp,Ymp
    int jmu=(Nmu-1)*(1+mu)/2;
    if(Xpm[jmu][0][0]>-100)
        return 0;
    double mu2=sqr(mu), fmu=fabs(mu);
    for(int ii=0;ii<2;ii++)
        for(int jj=0;jj<2;jj++){
            Xmp[jmu][ii][jj]=0; Ymp[jmu][ii][jj]=(ii==jj);
            Xpm[jmu][ii][jj]=0; Ypm[jmu][ii][jj]=(ii==jj);
        }
    if(mppm==1){
        double npm=np/nm, mucn2=1-1/sqr(npm), etapm2=1-sqr(npm)*(1-sqr(mu));
        if(etapm2>0 ){//&& mucn2>0){
            double etapm = sqrt(etapm2),
            aux1= sqr((fmu-npm*etapm)/(fmu+npm*etapm))/2,
            aux2= sqr((npm*fmu-etapm)/(npm*fmu+etapm))/2,
            auy1= sqr(npm)*2*npm*fmu*etapm/sqr(fmu+npm*etapm),
            auy2= sqr(npm)*2*npm*fmu*etapm/sqr(npm*fmu+etapm);
            if(npm>1){
                Xpm[jmu][0][0] = (mu2>mucn2)* (aux1+aux2) + (mu2<mucn2) ;
                Xpm[jmu][0][1] = (mu2>mucn2)* (aux1-aux2);
                Ypm[jmu][0][0] = (mu2>mucn2)* (auy1+auy2);
                Ypm[jmu][0][1] = (mu2>mucn2)* (auy1-auy2);
            } else {
                Xpm[jmu][0][0] = aux1+aux2;
                Xpm[jmu][0][1] = aux1-aux2;
                Ypm[jmu][0][0] = auy1+auy2;
                Ypm[jmu][0][1] = auy1-auy2;
            }
            Xpm[jmu][1][0]= Xpm[jmu][0][1];
            Xpm[jmu][1][1]= Xpm[jmu][0][0];
            Ypm[jmu][1][1]= Ypm[jmu][0][0];
            Ypm[jmu][1][0]= Ypm[jmu][0][1];
        }
    } else {
        double nmp=nm/np, mucn2=1-sqr(1/nmp), etamp2 = 1-sqr(nmp)*(1-mu2);
        if(etamp2>0){// && mucn2>0){
            double etamp=sqrt(etamp2),
            aux1= sqr((fmu-nmp*etamp)/(fmu+nmp*etamp))/2,
            aux2= sqr((nmp*fmu-etamp)/(nmp*fmu+etamp))/2,
            auy1= sqr(nmp)*2*nmp*fmu*etamp/sqr(fmu+nmp*etamp),
            auy2= sqr(nmp)*2*nmp*fmu*etamp/sqr(nmp*fmu+etamp);
            if(nmp>1){
                Xmp[jmu][0][0] = (mu2>mucn2)* (aux1+aux2) + (mu2<mucn2) ;
                Xmp[jmu][0][1] = (mu2>mucn2)* (aux1-aux2);
                Ymp[jmu][0][0] = (mu2>mucn2)* (auy1+auy2);
                Ymp[jmu][0][1] = (mu2>mucn2)* (auy1-auy2);
            } else {
                Xmp[jmu][0][0] = aux1+aux2;
                Xmp[jmu][0][1] = aux1-aux2;
                Ymp[jmu][0][0] = auy1+auy2;
                Ymp[jmu][0][1] = auy1-auy2;
            }
            Xmp[jmu][1][0]= Xmp[jmu][0][1];
            Xmp[jmu][1][1]= Xmp[jmu][0][0];
            Ymp[jmu][1][1]= Ymp[jmu][0][0];
            Ymp[jmu][1][0]= Ymp[jmu][0][1];
        }
    }
    return 0;
}

double Delta0Z[2][jmaxmax][Nmu];
double deltaS[2][2][2][2];//[matrix i][matrix j][power of mu][power of mu or eta]

doublet getDelta0Z(const int jnu, const int jmu){
    // compute the I0 and IZ contribution to Delat(mu)
    doublet delta;
    if(Y<0 || np==nm) return delta;
    double deltaa[2];
    double  kappanuj=kappanu[jnu],
            mu=-1.+jmu*dmu,
            nuj=nu[jnu];
    double  etapm2 = 1-sqr(np/nm)*(1-sqr(mu)),
            etamp2 = 1-sqr(nm/np)*(1-sqr(mu));
    for(int k=0;k<2;k++){
        if(mu>0 && etapm2>0){
            getXY(1,mu);
            deltaa[k] =- phi(kappanuj,0,Y,mu)*I0(mu,nuj)
                       + phi(kappanuj,Y,Z,mu)*Xpm[jmu][k][0]*IZ(-mu,nuj);
            double etapm=sqrt(etapm2);
            deltaa[k] += phi(kappanuj,0,Y,etapm)*Ypm[jmu][k][0]*I0(etapm,nuj);
        }
        if(mu<0 && etamp2>0){
            getXY(0,-mu);
            deltaa[k] = phi(kappanuj,Y,Z,-mu)*IZ(mu,nuj)
                      + phi(kappanuj,0,Y,-mu)*Xmp[jmu][k][0]*I0(-mu,nuj);
            double etamp=sqrt(etamp2);
             deltaa[k] -= phi(kappanuj,Y,Z,etamp)*Ymp[jmu][k][0]*IZ(-etamp,nuj);
        }
    }
    delta.d0=deltaa[0];
    delta.d1=deltaa[1];
    return delta;
}

doublet getDelta0Zapprox(const int i, const int jnu){
    double z=i*dz, nuj=nu[jnu], kappanuj=kappanu[jnu];
    int iY=Y/dz;
    doublet Aq;
    double F0=0, F2=0 ;
    if(i>iY){
        double nn=np/nm, Xp = 0.5*sqr((1-nn)/(1+nn)),
        Yp = 2*sqr(nn/(1+sqr(nn)))+0.5;
        double kaYz=(kappasum[i]-kappasum[iY])*kappanuj,
        ka0Y=kappasum[iY]*kappanuj, ka0z=kappasum[i]*kappanuj;
        F0 = Xp*Cs*Bsun(nuj)*phi(kappanuj,Y,Z,mus)*expint_E2(kaYz)
        + Ce*Bearth(nuj)*(Yp*expint_E3(ka0Y/nn+kaYz)-expint_E3(ka0z));  //removed from R
        F2 = Xp*Cs*Bsun(nuj)*phi(kappanuj,Y,Z,mus)*expint_E4(kaYz)
        + Ce*Bearth(nuj)*(Yp*expint_E5(ka0Y/nn+kaYz));//-expint_E5(ka0z));
        Aq.d0=F0; Aq.d1=F2;
    } else{
        double nn=nm/np, Xm = 0.5*sqr((1-nn)/(1+nn)),
        Ym = 2*sqr(nn/(1+sqr(nn)))+0.5;
        double kazY=(kappasum[iY]-kappasum[i])*kappanuj,ka0Y=kappasum[iY]*kappanuj;
        F0 = Cs*Bsun(nuj)*(phi(kappanuj,z,Y,mus)*expint_E2(kazY)
                           -Ym*phi(kappanuj,Y,Z,nn*mus)*expint_E2(kazY))
        - Ce*Bearth(nuj)*Xm*expint_E3(ka0Y+kazY);
        F2 = Cs*Bsun(nuj)*(phi(kappanuj,z,Y,mus)*expint_E4(kazY)
                           -Ym*phi(kappanuj,Y,Z,nn*mus)*expint_E4(kazY))
        -Ce*Bearth(nuj)*Xm*expint_E5(ka0Y+kazY);
        Aq.d0=-F0; Aq.d1=-F2;
    }
    return Aq;
}

double DeltaS[Nz][Nz][Nkappamax][2][2][2][2]; // will be used to tabulate deltaS

int getDeltaS(const double z, const int izp, const int jnu ){
    //   double deltaS[2][2][2][2]; [matrix i][matrix j][power q][power k]
    for(int q=0;q<2;q++)
        for(int i=0;i<2;i++)
            for(int j=0;j<2;j++)
                for(int k=0;k<2;k++)
                    for(int l=0;l<2;l++) deltaS[i][j][q][k]=0;
    if(Y<0 || np==nm) return 0;
   
    double zp=izp*dz, kappanuj=kappanu[jnu] ;
    if(z>Y)
       for (double mu=dmu/2;mu<1;mu+=dmu){
            getXY(1,mu);
            int jmu=(Nmu-1)*(1+mu)/2;
            double etapm2 = 1-sqr(np/nm)*(1-sqr(mu));
            for(int q=0;q<2;q++)
                for(int k=0;k<2;k++){
                    double auxkq = (q*sqr(mu)+1-q)*phi(kappanuj,Y,z,mu);
                    for(int i=0;i<2;i++)
                        for(int j=0;j<2;j++)
                            if(zp>Y && etapm2>0){
                                double etapm=sqrt(etapm2);
                               deltaS[i][j][q][k] +=
                                  auxkq*phi(kappanuj,Y,zp,mu)/mu
                                            *(k*sqr(mu)+1-k)*Xpm[jmu][i][j]*dmu
                                  - auxkq*(i==j)*phi(kappanuj,zp,Y,mu)/mu
                                                    *(k*sqr(mu)+1-k)*dmu
                                  + auxkq*phi(kappanuj,zp,Y,etapm)/etapm
                                            *(k*etapm2+1-k)*Ypm[jmu][i][j]*dmu;
                           }
                }
    }
    if(z<Y)
        for(double mu=-1+dmu/2;mu<0;mu+=dmu){
        double etamp2 = 1-sqr(nm/np)*(1-sqr(mu));
            getXY(0,mu);
            int jmu=(Nmu-1)*(1+mu)/2;
            for(int q=0;q<2;q++)
                for(int k=0;k<2;k++){
                    double auxkq = (q*sqr(mu)+1-q)*phi(kappanuj,z,Y,-mu);
                    for(int i=0;i<2;i++)
                        for(int j=0;j<2;j++)
                            if(zp<Y && etamp2>0){
                                double etamp=sqrt(etamp2);
                            deltaS[i][j][q][k] -=
                                auxkq*((i==j)*phi(kappanuj,Y,zp,-mu)
                                       *(k*sqr(mu)+1-k)/fabs(mu)
                                       + auxkq *phi(kappanuj,zp,Y,-mu) *(k*sqr(mu)+1-k)*Xmp[jmu][i][j]/fabs(mu))*dmu
                                    + auxkq*phi(kappanuj,Y,zp,etamp)
                                    *(k*etamp2+1-k)/etamp*Ymp[jmu][i][j]*dmu;
                    }
            }
        }
   return 0;
}

int updateJK(){
    // returns the convolution t-integral for the mean in mu of kappa*I_nu and mu^2 kappa*I_nu
    for(int jnu=0; jnu<jmax; jnu++){  // frequency loop
        double H[Nz],S[Nz];
        double kappanuj = kappanu[jnu], nuj=nu[jnu];
        double c0= - q0*mus*Cs*Bsun(nuj)*exp(-kappanuj*kappasum[Nz-1]/mus);
        
        for(int i=0;i<Nz;i++){ // initial altitude loop;
            double z=i*dz, kappazj=kappaz[i]*kappanuj,
            asz=as(z,nuj), sigs=kappazj*asz, siga=kappazj*(1-asz);
            S[i] = sigs*J0old[jnu][i]+siga*BB(nuj,T[i]) ;
            H[i] = 9*beta*sigs/8 * (J2old[jnu][i] -J0old[jnu][i]/3
                                    - K0old[jnu][i]+K2old[jnu][i]);
            c0 -= q0* expint_E2(kappasum[i]*kappanuj)*S[i]*dz;
        }
        
        for(int i=0;i<Nz;i++) { // main altitude loop
            double  ksz=kappasum[i]*kappanuj,
                    kszZ= kappasum[Nz-1]*kappanuj-ksz;
            // compute Rqnu
            double J0z  =Ce*Bearth(nuj)*expint_E3(ksz)/2
                        + c0*expint_E2(ksz)/2
                        + Cs*Bsun(nuj)*exp(-kszZ/mus)/2;
            double J2z  = Ce*Bearth(nuj)*expint_E5(ksz)/2
                           + c0*expint_E4(ksz)/2
                           + Cs*sqr(mus)*Bsun(nuj)*exp(-kszZ/mus)/2;
            double K0z=0, K2z=0;
            
            // compute Aqj, contribution to Delta of bdy conditions
            double z=i*dz, A[2][2];
            if(Y>0 && np!=nm){
                if(!approx){
                    A[0][0]=0;  A[1][0]=0; A[0][1]=0; A[1][1]=0;
                    for(int jmu = 0; jmu<Nmu; jmu++){
                        double mu=-1+jmu*dmu, aux=0;
                        if(z>Y && mu>0) aux = phi(kappanuj,Y,z,mu);
                        if(z<Y && mu<0) aux = -phi(kappanuj,z,Y,-mu);
                        for(int k=0;k<2;k++){
                            double bux = aux* Delta0Z[k][jnu][jmu];
                            for(int jj=0;jj<2;jj++)
                                A[k][jj] += bux*(1-jj + jj*sqr(mu))*dmu/2;
                        }
                    }
                    J0z += A[0][0];
                    J2z += A[0][1];
                    K0z += A[1][0];
                    K2z += A[1][1];
                }else{
                    doublet aux= getDelta0Zapprox(i,jnu);
                    J0z += aux.d0/2; J2z += aux.d1/2;
                }
                // compute contribution to Delta of source terms
                if(J0z<-1e-6 ){
                    cout<<jnu<<" "<<i<<" "<<J0z<<" "<<J2z<<" bug J0<0\n"; J0z=0; J2z=0;}
            }
                
                double A00=0, A20=0, A01=0, A21=0;
                for(int j=1;j<Nz;j++){ // convolution loop // add delta terms
                    double Hj=(H[j]+H[j-1])/2;
                    double Sj=(S[j]+S[j-1])/2-Hj/3;
                    double auxx = kappanuj*(kappasum[i]-(kappasum[j]+kappasum[j-1])/2);
                    double expE1ij = expint_E1(auxx);
                    double expE3ij = expint_E3(auxx);
                    double expE5ij = expint_E5(auxx);
                    J0z += (expE1ij*Sj + expE3ij*Hj)*dz/2;
                    J2z += (expE3ij*Sj + expE5ij*Hj)*dz/2;
                    K0z -= (expE1ij-expE3ij)*Hj*dz/2;
                    K2z -= (expE3ij-expE5ij)*Hj*dz/2;
                    // contribution of Delta
                    //getDeltaS(z,j,jnu); replaced by tabulation
                    if(Y>0 && np!=nm){
                        int Nkappaj= (kappanuj-kappamin)*Nkappamax/(kappamax-kappamin);
                        A00 += DeltaS[i][j][Nkappaj][0][0][0][0]*Sj // contrib to J0 of S10
                        + DeltaS[i][j][Nkappaj][0][0][0][1]*Hj // of S12
                        - DeltaS[i][j][Nkappaj][0][1][0][0]*Hj // of S20
                        + DeltaS[i][j][Nkappaj][0][1][0][1]*Hj; // of S22
                        A20 += DeltaS[i][j][Nkappaj][0][0][1][0]*Sj // contrib to J2 of S10
                        + DeltaS[i][j][Nkappaj][0][0][1][1]*Hj // of S12
                        - DeltaS[i][j][Nkappaj][0][1][1][0]*Hj // of S20
                        + DeltaS[i][j][Nkappaj][0][1][1][1]*Hj; // of S22
                        A01 += DeltaS[i][j][Nkappaj][1][1][0][0]*Sj // contrib to K0 of S10
                        + DeltaS[i][j][Nkappaj][1][0][0][1]*Hj // of S12
                        - DeltaS[i][j][Nkappaj][1][1][0][0]*Hj // of S20
                        + DeltaS[i][j][Nkappaj][1][1][0][1]*Hj; // of S22
                        A21 += DeltaS[i][j][Nkappaj][1][1][1][0]*Sj // contrib to K2 of S10
                        + DeltaS[i][j][Nkappaj][1][0][1][0]*Hj // of S12
                        - DeltaS[i][j][Nkappaj][1][0][1][0]*Hj // of S20
                        + DeltaS[i][j][Nkappaj][1][0][1][1]*Hj; // of S22
                        J0z += A00*dz/2;
                        J2z += A20*dz/2;
                        K0z += A01*dz/2;
                        K2z += A21*dz/2;
                    }
                }
//            if(J0z<-1e-6){
  //              cout<<jnu<<" "<<i<<" "<<J0z<<" "<<J2z<<" bug J0<0\n"; J0z=0; J2z=0;}
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
         myeq += (1-as(i*Z/(Nz-1),nu[j]))*kappanu[j]*BB((nu[j]+nu[j-1])/2,T0)*(nu[j]-nu[j-1]);
    return myeq;
}

double rhsTeq(const int i){
    double rhs=0;
    for(int j=1; j<jmax;j++)
        rhs+=(1-as(i*dz,nu[j]))*kappanu[j]*J0[j][i]*(nu[j]-nu[j-1]);
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
            J2old[j][i]=0;
            K0old[j][i]=0;
            K2old[j][i]=0;
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
        double normG=0, normS=0;
        for(int j=0; j<jmax;j++)
            for(int i=0;i<Nz;i++){ normG+=J0[j][i]*(nu[j]-nu[j-1])*dz;
                normS+=K2[j][i]*(nu[j]-nu[j-1])*dz; }
        cout << "k= "<<k <<"\t"<<T[2]*4798-273<<"\t"<<normG<<"\t"<<normS<<endl;
    }
     return 0;
}

doublet getIKzmu(const double z, const int jmu, const int jnu){  // UP only
    // Get I(z,mu) and K(z,mu) from J
    int jz=Nz*z/Z, jY=Nz*Y/Z;
    double mu=-1+jmu*dmu, kappanuj=kappanu[jnu], nuj=nu[jnu];
    int Nkappa = (kappanuj-kappamin)/(kappamax-kappamin)*Nkappamax;
    double kappazj=kappanuj*kappaz[jz],
            asz=as(z,nuj),
            sigs=kappazj*asz,
            siga=kappazj*(1-asz);
    double Izmu = 0, Kzmu=0;
    double H=J2[jnu][jz]-J0[jnu][jz]-K0[jnu][jz]+K2[jnu][jz];
    double Sj = siga*BB(nuj,T[jz]) + sigs*J0[jnu][jz]
                        +3*beta*sigs/8*(3*sqr(mu)-1)*H;
    double Kj =-9*beta*sigs/8*(1-sqr(mu))*H;
    
    if(mu>0){
        Izmu = phi(kappanuj,0.,z,mu)*I0(mu,nuj);
        if(z>Y){
            Izmu += phi(kappanuj,Y,z,mu)*(
                        Delta0Z[0][jnu][jmu]
                        +DeltaS[jY][jz][Nkappa][0][0][0][0]*Sj
                        +DeltaS[jY][jz][Nkappa][0][1][0][0]*Kj);
            Kzmu += phi(kappanuj,Y,z,mu)*(
                        Delta0Z[1][jnu][jmu]
                        +DeltaS[jY][jz][Nkappa][1][0][0][0]*Sj
                        +DeltaS[jY][jz][Nkappa][1][1][0][0]*Kj);
        }
        for(int jzp=0;jzp<jz;jzp++){
            H=J2[jnu][jzp]-J0[jnu][jzp]-K0[jnu][jzp]+K2[jnu][jzp];
            Sj = siga*BB(nuj,T[jzp]) + sigs*J0[jnu][jzp]
                                +3*beta*sigs/8*(3*sqr(mu)-1)*H;
            Izmu += phi(kappanuj,jzp*dz,z,mu)*Sj*dz;
            Kzmu += phi(kappanuj,jzp*dz,z,mu)*Kj*dz;
        }
    } else { // mu<0
        Izmu = phi(kappanuj,z,Z,-mu)*IZ(mu,nuj);
        if(z<Y){
            Izmu -= phi(kappanuj,z,Y,-mu)*(
                                Delta0Z[0][jnu][jmu]
                               +DeltaS[jY][jz][Nkappa][0][0][0][0]*Sj
                               +DeltaS[jY][jz][Nkappa][0][1][0][0]*Kj);
            Kzmu -= phi(kappanuj,z,Y,-mu)*(
                                Delta0Z[1][jnu][jmu]
                                +DeltaS[jY][jz][Nkappa][1][0][0][0]*Sj
                                +DeltaS[jY][jz][Nkappa][1][1][0][0]*Kj);
        }
        for(int jzp=jz;jzp<Nz;jzp++){
                H=J2[jnu][jzp]-J0[jnu][jzp]-K0[jnu][jzp]+K2[jnu][jzp];
                Sj = siga*BB(nuj,T[jzp]) + sigs*J0[jnu][jzp]
                +3*beta*sigs/8*(3*sqr(mu)-1)*H;
                Izmu += phi(kappanuj,z,jzp*dz,-mu)*Sj*dz;
                Kzmu += phi(kappanuj,jzp*dz,z,-mu)*Kj*dz;
            }
}
//    if(mu<0 && z<Y && Izmu<-1e-7)
  //      cout<<"bug";
    doublet dd; dd.d0=Izmu; dd.d1=Kzmu;
    return dd;
}

int main(int argc, const char * argv[]) {
    //----------------------------------
    
    jmax=readkappa(mykappafile);
    cout<<"Prepare the tabulations \n";
    for(int jmu=0;jmu<Nmu;jmu++) Xpm[jmu][0][0]=-200; // tabulation of Xpm, Ypm

    for(int iz=0;iz<Nz;iz++)
        for(int izp=0;izp<Nz;izp++)
            for(int Nkappaj=0;Nkappaj<Nkappamax;Nkappaj++){
                int jnu=jmax*double(Nkappaj)/Nkappamax;
                getDeltaS(iz*dz,izp,jnu);
                for(int q=0;q<2;q++)
                    for(int i=0;i<2;i++)
                        for(int j=0;j<2;j++)
                            for(int k=0;k<2;k++)
                                for(int l=0;l<2;l++)
                                    DeltaS[iz][izp][Nkappaj][i][j][q][k]
                                        =deltaS[i][j][q][k];
            }
    for(int j=0;j<Nkappamax;j++)
        for(int k=0;k<Nz;k++)
            for(int l=0;l<Nz;l++)
                for(int i=0;i<Nmu;i++)
                    PHI[j][k][l][i]=-200; // arbitrary value never reached
    
    
    for(int jnu=0; jnu<jmax; jnu++) // Tabluation
        for(int jmu=0;jmu<Nmu;jmu++){
            doublet aux =getDelta0Z(jnu,jmu);
            Delta0Z[0][jnu][jmu] = aux.d0;
            Delta0Z[1][jnu][jmu] = aux.d1;
        }
    cout<<"\n Tabulations done\n";

//    { int K=2;
    for(int K=0;K<2;K++){
        for(int j=0;j<jmax;j++)
            if (K==1){
                if((nu[j]>3./18)&&(nu[j]<3./14))
                    kappanu[j]=fmax(0.5,1.5*kappanu[j]);//Visible R
                if( kappanu[j]>kappamax) kappamax= kappanu[j];
                if( kappanu[j]<kappamin) kappamin= kappanu[j];
            }
            else if (K==2) kappanu[j]=0.5;  // constant kappa
        
        ofstream  resultfile =
        ofstream(myresulttemperature+to_string(cloud)+to_string(K)+".txt");
        cout<<"\n iterations  T[2*dz]  norm of J0  and of K2\n";
        double t0 = clock();
        //-----------------------
        multiBlock(Te/2);
        //-----------------------
        printf( " Time CPU = %10.6f\n", (clock() - t0)/CLOCKS_PER_SEC);
        cout<<"\n tau\t [T]:"<<endl;
        for(int i=1;i<Nz;i++){
            cout<< i*dz <<" "<<T[i]*T0-273<<endl;
            resultfile<< i*dz <<" "<<T[i]*T0-273<<endl;
        }
        resultfile.close();
        double J0total=0, aux = Cs*sqr(sqr(pi*Ts))/30;
        for(int j=0;j<jmax;j++)
            J0total += J0[j][Nz-1]*(nu[j+1]-nu[j]);
        cout<<"Radiative Forcing ="<< J0total-aux<<" "<<J0total <<" "<<aux<< endl;

    }
    // display a surface I(z,mu) at given kappa[j], use splot "test.txt" matrix w l
        int j=262; // wavelength 15.9175
        ofstream  ff =  ofstream(basedir+"gnuplotI00.txt");
        ofstream  gg =  ofstream(basedir+"gnuplotQ00.txt");
        for(double z=dz;z<Z;z+=2*dz){
            for(double mu=-1+dmu;mu<1-dmu;mu+=2*dmu){
                int jmu=(Nmu-1)*(1+mu)/2;
                doublet aux = getIKzmu(z,jmu,j);
                ff<<z<<"  "<<mu<<"  "<<
                 aux.d0<<endl;//*((z<Y)+(z>Y)*sqr(np/nm)) <<endl;
                gg<<z<<"  "<<mu<<"  "<<aux.d1<<endl;//*((z<Y)+(z>Y)*sqr(np/nm))<<endl;
            }
            ff<<endl; gg<<endl;
        }
     /*
      Display surfaces by
      gnuplot> splot "gnuplotI.txt" w l
      gnuplot> splot "gnuplotK.txt" w l
      */
return 0;
}
// file names for results:  temperature+cloud+K with cloud=0/1, K=0 K=1 (CO2) K=2 (constant kappa)
