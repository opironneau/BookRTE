
//  Created by Olivier Pironneau on 01/03/2026.
//
// kappa file needs be in the directory "basedir"
// compile with g++ -O3 -std=c++17 -DIQn_HAVE_GSL=0 -o IQn IQn.cpp

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
using namespace std;
#define sqr(x) ((x)*(x))
const double pi = atan(1.)*4;

const double Ce= 2.4, Cs=2e-6 ;
const double rho0=1.,drho=-0.7; // density gradient
const double Z=1.2;
const double B0=1.4744e-8, T0=4798;

const double Te=(273+18)/T0, Ts=5778/T0;
const double q0=-0.3;
const double mus=0.5;
const double _beta=0.5;

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

const double Y=0.5*Z, nm=1, np=0.9; // must match jump of n(). See also kapparho
const double n(double z) {  return  nm + (np-nm)*(z>Y);}  // diffraction index

const double kapparhoz(double z){ return (1+drho*z)*(1*(z<Y) + 1*(z>=Y));}// density

struct doublet{ double d0=0,d1=0;};

double  kappanu[jmaxmax], nu[jmaxmax];      // kappanu and uneven discretization of [numin,numax]
double  J0[jmaxmax][Nz], J0old[jmaxmax][Nz],// mu integral of I_nu and
        J2[jmaxmax][Nz], J2old[jmaxmax][Nz], // mu integral of mu^2 I_nu
        K0[jmaxmax][Nz], K0old[jmaxmax][Nz], // mu integral of mu^2 I_nu
        K2[jmaxmax][Nz], K2old[jmaxmax][Nz], // mu integral of mu^2 I_nu
        T[Nz], nn[Nz],  // Delta for jump and Temeprature
        c0nu[jmaxmax],                      // c0 of (4.13), kept for the I(0,mu) bdy cond
        kappaz[Nz], kappasum[Nz+1];         // kappa(z) and int_0^z(kappa(z))

string basedir("");
string mykappafile(basedir+"_kappa.txt"); // 1 - kappa
string myresulttemperature(basedir+"temperaturegg"); // use z with np<nm, w otherwise
string myresultmeanintensity(basedir+"imean0");
int jmax;


const double as(double z, double nu){// scattering
    return 0.3                                          // background scattering
    + 0.3*4*(z2-z)*(z-z1)*(z>z1)*(z<z2)/sqr(z1-z2)      // cloud scattering
    + 0.3*16*(nu<nu2)*(nu>nu1)*sqr((nu-nu1)*(nu-nu2)/sqr(nu1-nu2))*(z>z3);}  // Rayleigh


/********************************** END of DECLARATIONS ****************************************/

double expint_E1(const double t=1){
     double t1=fabs(t);
    const int Kexpint = 12+t1*4; // precision in exponential integral function E1
    const double  gaNtaua =0.577215664901533; // special integration for log(t)
    // E1 has an integrable log singularity at 0; returning 0 here loses that mass.
    if(t1<1e-5) return -gaNtaua - log(fmax(t1,1e-12));
    if(t1>12) {   // below 12 the series is more accurate than the asymptotic expansion
        if(verbose) cout << "argument of E_1>12"<<endl;
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
    return (exp(-t1) - t1*expint_E3(t1))/3;  // t1, not t: E_n is called with t<0 in updateJK
}
double expint_E5(const double t=1){
    double t1=fabs(t);
    return (exp(-t1) - t1*expint_E4(t1))/4;  // idem
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
    int j=0;
    double kappaux, waveno;
    // test the read itself (eof() is only set *after* a failed read) and stay in bounds
    while((j<jmaxmax) && (kappafile >> waveno >> kappaux)){
        kappanu[j]  = fmin(3.,fmax(kappaux,0.01));
        nu[j]=3/waveno;
        // bracket the *clamped* values: kappanu[] is what indexes PHI and DeltaS
        if(kappanu[j]>kappamax) kappamax=kappanu[j];
        if(kappanu[j]<kappamin) kappamin=kappanu[j];
        j++;
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

// kappanuj -> bin index, clamped:
inline int knubin(const double kappanuj){
    int k = (kappanuj-kappamin)*Nkappamax/(kappamax-kappamin);
    return k<0 ? 0 : (k>=Nkappamax ? Nkappamax-1 : k);
}

inline double phi(const double kappanuj,const double zp,const double zpp,const double mu)
{
    // compute exp(-kappanuj\int_zp^zpp kapparho(z")/mu(z")dz") with mu(z") as in mu2 below
    if(zpp-zp < 1e-6) return 1;
    if(zp>zpp || mu<=0){
        cout << "bug in calling phi\n";
        return 1;
    }
    int izp=zp/dz, izpp=zpp/dz, jmu=(Nmu-1)*(1+mu)/2, knuj=knubin(kappanuj);
    if(PHI[knuj][izp][izpp][jmu]>=0)
        return PHI[knuj][izp][izpp][jmu]; // this saves computing exponentials
    else{
        double aux = fabs((kappasum[izpp]-kappasum[izp])*kappanuj/mu);
        if(aux>100) return PHI[knuj][izp][izpp][jmu] = 0.;
        PHI[knuj][izp][izpp][jmu] = exp(-aux);
        return PHI[knuj][izp][izpp][jmu];
    }
}
int resetPHI(){ // must be redone whenever kappanu[], kappamin or kappamax change
    for(int j=0;j<Nkappamax;j++)
        for(int k=0;k<Nz;k++)
            for(int l=0;l<Nz;l++)
                for(int i=0;i<Nmu;i++)
                    PHI[j][k][l][i]=-200; // arbitrary value never reached
    return 0;
}
// first index is a tabulation; done[jmu][mppm] says whether that family is already computed.
// One flag per family: the old test on Xpm[jmu][0][0] made the *first* call at a given jmu
// mask every later one, so only one of the two families was ever really computed.
double Xpm[Nmu][2][2], Ypm[Nmu][2][2], Xmp[Nmu][2][2],Ymp[Nmu][2][2];
bool doneXY[Nmu][2];

int getXY(const int mppm, const double mu){
    // Compute Xpm,Ypm if mppm=1 else Xmp,Ymp. Indexed by mu itself (may be <0).
    int jmu=(Nmu-1)*(1+mu)/2;
    if(doneXY[jmu][mppm])
        return 0;
    doneXY[jmu][mppm]=true;
    double mu2=sqr(mu), fmu=fabs(mu);
    for(int ii=0;ii<2;ii++)
        for(int jj=0;jj<2;jj++){
            if(mppm==1){ Xpm[jmu][ii][jj]=0; Ypm[jmu][ii][jj]=(ii==jj); }
            else       { Xmp[jmu][ii][jj]=0; Ymp[jmu][ii][jj]=(ii==jj); }
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
        } else {
            // Total internal reflection: etapm2<=0 <=> mu2<mucn2. X=Gamma=Id, Y=0.
            // (the (mu2<mucn2) term above is unreachable, both conditions are the same)
            for(int ii=0;ii<2;ii++)
                for(int jj=0;jj<2;jj++){
                    Xpm[jmu][ii][jj]=(ii==jj); Ypm[jmu][ii][jj]=0;
                }
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
        } else {
            // Total internal reflection (here nmp=nm/np=1.43 so |mu|<0.714): X=Id, Y=0
            for(int ii=0;ii<2;ii++)
                for(int jj=0;jj<2;jj++){
                    Xmp[jmu][ii][jj]=(ii==jj); Ymp[jmu][ii][jj]=0;
                }
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
    double deltaa[2]={0,0};  // must be initialized: mu==0, and total internal reflection,
                             // leave both branches below untaken
    double  kappanuj=kappanu[jnu],
            mu=-1.+jmu*dmu,
            nuj=nu[jnu];
    double  etapm2 = 1-sqr(np/nm)*(1-sqr(mu)),
            etamp2 = 1-sqr(nm/np)*(1-sqr(mu));
    // getXY is called with the *same* mu used to index Xpm/Xmp below; it used to be called
    // with -mu for mu<0, which keys the tabulation on the mirrored index.
    if(mu>0) getXY(1,mu);
    if(mu<0) getXY(0,mu);
    for(int k=0;k<2;k++){
        if(mu>0){  // Delta(mu)|_{mu>0}, eq. (5.15)
            deltaa[k] =- phi(kappanuj,0,Y,mu)*I0(mu,nuj)
                       + phi(kappanuj,Y,Z,mu)*Xpm[jmu][k][0]*IZ(-mu,nuj);
            if(etapm2>0){ // under total internal reflection Ypm==0 and eta does not exist
                double etapm=sqrt(etapm2);
                deltaa[k] += phi(kappanuj,0,Y,etapm)*Ypm[jmu][k][0]*I0(etapm,nuj);
            }
        }
        if(mu<0){  // Delta(mu)|_{mu<0}: the X_mp term carries a MINUS sign in (5.15)
            deltaa[k] = phi(kappanuj,Y,Z,-mu)*IZ(mu,nuj)
                      - phi(kappanuj,0,Y,-mu)*Xmp[jmu][k][0]*I0(-mu,nuj);
            if(etamp2>0){
                double etamp=sqrt(etamp2);
                deltaa[k] -= phi(kappanuj,Y,Z,etamp)*Ymp[jmu][k][0]*IZ(-etamp,nuj);
            }
        }
    }
    delta.d0=deltaa[0];
    delta.d1=deltaa[1];
    return delta;
}

double DeltaS[Nz][Nz][Nkappamax][2][2][2][2]; // will be used to tabulate deltaS

int getDeltaS(const double z, const int izp, const double kappanuj ){
    //   double deltaS[2][2][2][2]; [matrix i][matrix j][power q][power k]
    // Each of the 3 terms of Delta lives on its own side of Y: the X_pm integral is on
    // (Y,Z), the identity and Y_pm integrals on (0,Y) -- and the mirror image for z<Y.
    // A single "zp>Y" test used to gate all three, which sent phi() outside its interval
    // (where it silently returns 1).
    for(int q=0;q<2;q++)
        for(int i=0;i<2;i++)
            for(int j=0;j<2;j++)
                for(int k=0;k<2;k++) deltaS[i][j][q][k]=0;
    if(Y<0 || np==nm) return 0;

    double zp=izp*dz;
    if(z>Y)
       for (double mu=dmu/2;mu<1;mu+=dmu){
            getXY(1,mu);
            int jmu=(Nmu-1)*(1+mu)/2;
            double etapm2 = 1-sqr(np/nm)*(1-sqr(mu));
            for(int q=0;q<2;q++)
                for(int k=0;k<2;k++){
                    double auxkq = (q*sqr(mu)+1-q)*phi(kappanuj,Y,z,mu),
                           mu2k  = (k*sqr(mu)+1-k);
                    for(int i=0;i<2;i++)
                        for(int j=0;j<2;j++){
                            if(zp>Y)  // + X_pm \int_Y^Z phi(Y,z') S
                                deltaS[i][j][q][k] +=
                                  auxkq*phi(kappanuj,Y,zp,mu)/mu*mu2k*Xpm[jmu][i][j]*dmu;
                            if(zp<Y){ // - \int_0^Y phi(z',Y) S
                                deltaS[i][j][q][k] -=
                                  auxkq*(i==j)*phi(kappanuj,zp,Y,mu)/mu*mu2k*dmu;
                                if(etapm2>0){ // + Y_pm \int_0^Y phi_eta(z',Y) S_eta
                                    double etapm=sqrt(etapm2);
                                    deltaS[i][j][q][k] +=
                                      auxkq*phi(kappanuj,zp,Y,etapm)/etapm
                                            *(k*etapm2+1-k)*Ypm[jmu][i][j]*dmu;
                                }
                            }
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
                    // A^q|_{z<Y} = -phi(z,Y) Delta|_{mu<0}, so the signs of the X_mp and
                    // Y_mp terms flip relative to the identity term.
                    double auxkq = (q*sqr(mu)+1-q)*phi(kappanuj,z,Y,-mu),
                           mu2k  = (k*sqr(mu)+1-k), fmu=fabs(mu);
                    for(int i=0;i<2;i++)
                        for(int j=0;j<2;j++){
                            if(zp>Y){ // - \int_Y^Z phi(Y,z') S ; + Y_mp \int_Y^Z phi_eta
                                deltaS[i][j][q][k] -=
                                  auxkq*(i==j)*phi(kappanuj,Y,zp,-mu)/fmu*mu2k*dmu;
                                if(etamp2>0){
                                    double etamp=sqrt(etamp2);
                                    deltaS[i][j][q][k] +=
                                      auxkq*phi(kappanuj,Y,zp,etamp)/etamp
                                            *(k*etamp2+1-k)*Ymp[jmu][i][j]*dmu;
                                }
                            }
                            if(zp<Y)  // + X_mp \int_0^Y phi(z',Y) S
                                deltaS[i][j][q][k] += 	
                                auxkq*phi(kappanuj,zp,Y,-mu)/fmu*mu2k*Xmp[jmu][i][j]*dmu;
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
            H[i] = 9*_beta*sigs/8 * (J2old[jnu][i] -J0old[jnu][i]/3
                                    - K0old[jnu][i]+K2old[jnu][i]);
            S[i] = sigs*J0old[jnu][i]+siga*BB(nuj,T[i]) - H[i]/ 3.0 ;
            // (4.13): c0 also carries an E_4(kappa(0,z)) H term, which was missing
            c0 -= q0*( expint_E2(kappasum[i]*kappanuj)*S[i]
                     + expint_E4(kappasum[i]*kappanuj)*H[i] )*dz;
        }
        
        c0nu[jnu]=c0;
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
                // compute contribution to Delta of source terms
                if(J0z<-1e-6 ){
                    cout<<jnu<<" "<<i<<" "<<J0z<<" "<<J2z<<" bug J0<0\n"; J0z=0; J2z=0;}
            }
                
                double A00=0, A20=0, A01=0, A21=0;
                for(int j=1;j<Nz;j++){ // convolution loop // add delta terms
                    double Hj=(H[j]+H[j-1])/2;
                    double Sj=(S[j]+S[j-1])/2;
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
                        int Nkappaj= knubin(kappanuj);
                        // index order is DeltaS[.][.][.][row][col][q][k]: the row selects
                        // J (0) or K (1), the col selects S_1 (0) or S_2 (1).
                        A00 += DeltaS[i][j][Nkappaj][0][0][0][0]*Sj // contrib to J0 of S10
                        + DeltaS[i][j][Nkappaj][0][0][0][1]*Hj // of S12
                        - DeltaS[i][j][Nkappaj][0][1][0][0]*Hj // of S20
                        + DeltaS[i][j][Nkappaj][0][1][0][1]*Hj; // of S22
                        A20 += DeltaS[i][j][Nkappaj][0][0][1][0]*Sj // contrib to J2 of S10
                        + DeltaS[i][j][Nkappaj][0][0][1][1]*Hj // of S12
                        - DeltaS[i][j][Nkappaj][0][1][1][0]*Hj // of S20
                        + DeltaS[i][j][Nkappaj][0][1][1][1]*Hj; // of S22
                        A01 += DeltaS[i][j][Nkappaj][1][0][0][0]*Sj // contrib to K0 of S10
                        + DeltaS[i][j][Nkappaj][1][0][0][1]*Hj // of S12
                        - DeltaS[i][j][Nkappaj][1][1][0][0]*Hj // of S20
                        + DeltaS[i][j][Nkappaj][1][1][0][1]*Hj; // of S22
                        A21 += DeltaS[i][j][Nkappaj][1][0][1][0]*Sj // contrib to K2 of S10
                        + DeltaS[i][j][Nkappaj][1][0][1][1]*Hj // of S12
                        - DeltaS[i][j][Nkappaj][1][1][1][0]*Hj // of S20
                        + DeltaS[i][j][Nkappaj][1][1][1][1]*Hj; // of S22
                    }
                }
                // A00..A21 are the *whole* z'-integral: adding them inside the j-loop
                // above counted the j-th contribution Nz-j times.
                if(Y>0 && np!=nm){
                    J0z += A00*dz/2;
                    J2z += A20*dz/2;
                    K0z += A01*dz/2;
                    K2z += A21*dz/2;
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
        for(int j=1; j<jmax;j++)   // j=0 read nu[-1]
            for(int i=0;i<Nz;i++){ normG+=J0[j][i]*(nu[j]-nu[j-1])*dz;
                normS+=K2[j][i]*(nu[j]-nu[j-1])*dz; }
        cout << "k= "<<k <<"\t"<<T[2]*4798-273<<"\t"<<normG<<"\t"<<normS<<endl;
    }
     return 0;
}

// S(z',m) = [S_1,S_2] of (SS), at node iz and direction cosine m>0. sigma_a, sigma_s and
// H all depend on z', so they must be evaluated at iz, not frozen at the observation point.
doublet source(const int iz, const double m, const int jnu){
    double nuj=nu[jnu], kappazj=kappanu[jnu]*kappaz[iz],
           asz=as(iz*dz,nuj), sigs=kappazj*asz, siga=kappazj*(1-asz);
    double Ht = J2[jnu][iz]-J0[jnu][iz]/3-K0[jnu][iz]+K2[jnu][iz]; // the 1/3 was missing
    doublet s;
    s.d0 = siga*BB(nuj,T[iz]) + sigs*J0[jnu][iz] + 3*_beta*sigs/8*(3*sqr(m)-1)*Ht;
    s.d1 = -9*_beta*sigs/8*(1-sqr(m))*Ht;
    return s;
}

// Source part of Delta(mu), eq. (5.15), at one given mu. DeltaS cannot be reused here:
// it is the mu-*integrated* A^q, not Delta(mu) itself.
int getDeltaSourceAtMu(const double mu, const int jnu, double d[2]){
    d[0]=d[1]=0;
    if(Y<0 || np==nm) return 0;
    double kappanuj=kappanu[jnu], fmu=fabs(mu);
    int jmu=(Nmu-1)*(1+mu)/2;
    if(mu>0){
        getXY(1,mu);
        double etapm2=1-sqr(np/nm)*(1-sqr(mu)), etapm=etapm2>0?sqrt(etapm2):0;
        for(int jzp=0;jzp<Nz;jzp++){
            double zp=jzp*dz;
            if(zp>Y){        // + X_pm \int_Y^Z phi(Y,z') S
                doublet s=source(jzp,mu,jnu);
                double w=phi(kappanuj,Y,zp,mu)/mu*dz;
                for(int k=0;k<2;k++) d[k]+= w*(Xpm[jmu][k][0]*s.d0+Xpm[jmu][k][1]*s.d1);
            } else {         // - \int_0^Y phi(z',Y) S
                doublet s=source(jzp,mu,jnu);
                double w=phi(kappanuj,zp,Y,mu)/mu*dz;
                d[0]-= w*s.d0; d[1]-= w*s.d1;
                if(etapm2>0){ // + Y_pm \int_0^Y phi_eta(z',Y) S_eta
                    doublet se=source(jzp,etapm,jnu);
                    double we=phi(kappanuj,zp,Y,etapm)/etapm*dz;
                    for(int k=0;k<2;k++)
                        d[k]+= we*(Ypm[jmu][k][0]*se.d0+Ypm[jmu][k][1]*se.d1);
                }
            }
        }
    } else if(mu<0){
        getXY(0,mu);
        double etamp2=1-sqr(nm/np)*(1-sqr(mu)), etamp=etamp2>0?sqrt(etamp2):0;
        for(int jzp=0;jzp<Nz;jzp++){
            double zp=jzp*dz;
            if(zp>Y){        // + \int_Y^Z phi(Y,z') S
                doublet s=source(jzp,fmu,jnu);
                double w=phi(kappanuj,Y,zp,fmu)/fmu*dz;
                d[0]+= w*s.d0; d[1]+= w*s.d1;
                if(etamp2>0){ // - Y_mp \int_Y^Z phi_eta(Y,z') S_eta
                    doublet se=source(jzp,etamp,jnu);
                    double we=phi(kappanuj,Y,zp,etamp)/etamp*dz;
                    for(int k=0;k<2;k++)
                        d[k]-= we*(Ymp[jmu][k][0]*se.d0+Ymp[jmu][k][1]*se.d1);
                }
            } else {         // - X_mp \int_0^Y phi(z',Y) S
                doublet s=source(jzp,fmu,jnu);
                double w=phi(kappanuj,zp,Y,fmu)/fmu*dz;
                for(int k=0;k<2;k++) d[k]-= w*(Xmp[jmu][k][0]*s.d0+Xmp[jmu][k][1]*s.d1);
            }
        }
    }
    return 0;
}

doublet getIKzmu(const double z, const int jmu, const int jnu){
    // I(z,mu) and Q(z,mu) from (gensolppm)
    doublet dd={0.,0.};
    double mu=-1+jmu*dmu, kappanuj=kappanu[jnu], nuj=nu[jnu], fmu=fabs(mu);
    if(fmu<0.02) return dd;   // mu=0 is a singular direction of the characteristics
    double Izmu=0, Kzmu=0, dsrc[2]={0,0}, dtot[2]={0,0};
    if(Y>0 && np!=nm){
        getDeltaSourceAtMu(mu,jnu,dsrc);
        dtot[0]=Delta0Z[0][jnu][jmu]+dsrc[0];   // boundary part + source part
        dtot[1]=Delta0Z[1][jnu][jmu]+dsrc[1];
    }
    if(mu>0){
        Izmu = phi(kappanuj,0.,z,mu)*(I0(mu,nuj)+c0nu[jnu]/2); // I(0,mu)=c'_e B mu + c0/2
        if(z>Y){
            Izmu += phi(kappanuj,Y,z,mu)*dtot[0];
            Kzmu += phi(kappanuj,Y,z,mu)*dtot[1];
        }
        for(int jzp=0;jzp<Nz;jzp++){       // \int_0^z phi(z',z) S(z') dz', S=S(z',mu)/mu
            double zp=jzp*dz;
            if(zp>=z) break;
            doublet s=source(jzp,mu,jnu);
            double w=phi(kappanuj,zp,z,mu)/mu*dz;
            Izmu += w*s.d0; Kzmu += w*s.d1;
        }
    } else { // mu<0
        Izmu = phi(kappanuj,z,Z,fmu)*IZ(mu,nuj);
        if(z<Y){
            Izmu -= phi(kappanuj,z,Y,fmu)*dtot[0];
            Kzmu -= phi(kappanuj,z,Y,fmu)*dtot[1];
        }
        for(int jzp=0;jzp<Nz;jzp++){       // \int_z^Z phi(z,z') S(z') dz'
            double zp=jzp*dz;
            if(zp<=z) continue;
            doublet s=source(jzp,fmu,jnu);
            double w=phi(kappanuj,z,zp,fmu)/fmu*dz;
            Izmu += w*s.d0; Kzmu += w*s.d1;
        }
    }
    dd.d0=Izmu; dd.d1=Kzmu;
    return dd;
}

int tabulate(){
    // Must be redone whenever kappanu[], kappamin or kappamax change (i.e. at each K).
    // resetPHI() MUST come first: PHI is a static array, so it starts at 0, and phi()
    // treats any value >=0 as a valid cache hit -- every phi() call inside the DeltaS
    // tabulation used to return 0, which made the whole source part of Delta vanish.
    resetPHI();
    for(int Nkappaj=0;Nkappaj<Nkappamax;Nkappaj++){
        // the bin must be tabulated with the kappa it is later looked up with; picking
        // kappanu[jmax*Nkappaj/Nkappamax] filed each bin under an unrelated frequency.
        double kappanuj = kappamin + (kappamax-kappamin)*(Nkappaj+0.5)/Nkappamax;
        for(int iz=0;iz<Nz;iz++)
            for(int izp=0;izp<Nz;izp++){
                getDeltaS(iz*dz,izp,kappanuj);
                for(int q=0;q<2;q++)
                    for(int i=0;i<2;i++)
                        for(int j=0;j<2;j++)
                            for(int k=0;k<2;k++)
                                DeltaS[iz][izp][Nkappaj][i][j][q][k]
                                    =deltaS[i][j][q][k];
            }
    }
    for(int jnu=0; jnu<jmax; jnu++) // Tabluation
        for(int jmu=0;jmu<Nmu;jmu++){
            doublet aux =getDelta0Z(jnu,jmu);
            Delta0Z[0][jnu][jmu] = aux.d0;
            Delta0Z[1][jnu][jmu] = aux.d1;
        }
    return 0;
}

int main(int argc, const char * argv[]) {
    //----------------------------------
    if(argc>1){ // optional output/input directory, so the code runs out of the box
        basedir=string(argv[1]);
        if(basedir.back()!='/') basedir+="/";
        mykappafile=basedir+"_kappa.txt";
        myresultmeanintensity=basedir+"imean0";
    }
    jmax=readkappa(mykappafile);

    { int K=1;
//    for(int K=0;K<2;K++){
        for(int j=0;j<jmax;j++)
            if (K==1){
                if((nu[j]>3./18)&&(nu[j]<3./14))
                    kappanu[j]=fmax(0.5,1.5*kappanu[j]);//Visible R
                if( kappanu[j]>kappamax) kappamax= kappanu[j];
                if( kappanu[j]<kappamin) kappamin= kappanu[j];
            }
            else if (K==2) kappanu[j]=0.5;  // constant kappa
        if(K==2){ kappamin=0.5; kappamax=0.5+1e-12; }
        cout<<"Prepare the tabulations (K="<<K<<") \n";
        tabulate();   // kappanu[] has just changed, so PHI and DeltaS must be rebuilt
        cout<<" Tabulations done\n";

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
        for(int j=0;j<jmax-1;j++)   // j=jmax-1 read nu[jmax]
            J0total += J0[j][Nz-1]*(nu[j+1]-nu[j]);
        cout<<"Radiative Forcing ="<< J0total-aux<<" "<<J0total <<" "<<aux<< endl;

    }
    // display a surface I(z,mu) at given kappa[j], use splot "test.txt" matrix w l
//        int j=515; // wavelength 1.9204
        ofstream  ff =  ofstream(basedir+"gnuplotI00.txt");
        ofstream  gg =  ofstream(basedir+"gnuplotQ00.txt");
    
        for(double z=dz;z<Z;z+=2*dz){
            for(double mu=-1+dmu;mu<1-dmu;mu+=2*dmu){
                int jmu=(Nmu-1)*(1+mu)/2;
                doublet aux1, aux={0.,0.};
                for(int j=1;j<jmax-1;j++){
                    aux1= getIKzmu(z,jmu,j);
                    aux.d0 += aux1.d0; aux.d1 += aux1.d1;
                }
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
