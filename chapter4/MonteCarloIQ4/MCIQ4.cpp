//  Monte Carlo solution of the polarized VRTE of Chapter 4: (4.9),(4.10),(4.11)
//
//    mu dI/dz + kappa I = sigma_a B_nu(T) + (sigma_s/2) int_{-1}^1 I dmu'
//                       + (beta sigma_s/4) P_2(mu) int_{-1}^1 [P_2 I - (1-P_2)Q] dmu'
//    mu dQ/dz + kappa Q = -(beta sigma_s/4)(1-P_2(mu)) int_{-1}^1 [P_2 I - (1-P_2)Q] dmu'
//    int_0^oo sigma_a ( B_nu(T) - (1/2) int_{-1}^1 I dmu ) dnu = 0
//    I(0,mu)|_{mu>0} = c'_e B_nu(T_e) mu + q0 int_{-1}^0 mu I(0,mu) dmu,  Q(0,mu)=0,
//    I(Z,mu)|_{mu<0} = c_s B_nu(T_s) delta(mu+mu_s),                      Q(Z,mu)=0.
//
//  The photons are followed in the equivalent variables (4.8), I_l=(I+Q)/2, I_r=(I-Q)/2,
//  because there the scattering kernel K_{p p'}(mu,mu') is nonnegative and normalized,
//        sum_{p} int_{-1}^1 K_{p p'}(mu,mu') dmu = 1  for p'= l,r ,
//  i.e. it is the probability density of the (direction, polarization state) after a
//  collision: an analog Monte Carlo is then possible, no negative weight is needed.
//  A photon carries the pair of weights (w_l,w_r): it scores w_l+w_r in I and w_l-w_r in Q.
//  Only the direction is sampled, from the density h(mu)=sum_p c_p(mu)/(w_l+w_r) with
//  c_p(mu) = K_{p l}(mu,mu')w_l + K_{p r}(mu,mu')w_r; the new weights are then
//  w_p = (w_l+w_r) c_p(mu)/(c_l(mu)+c_r(mu)), so the total weight is exactly conserved
//  and the polarization carries no sampling noise -- Q is 1 to 2% of I, an l/r random
//  tag would need some 10^4 times more photons for the same accuracy.
//
//  Three sources: the collimated sun at z=Z, the Lambert emission of the ground at z=0
//  and the volume emission 2 sigma_a B_nu(T); the ground reflects with the albedo
//  |q0|/2 and a cosine law. Only the temperature needs outer iterations: for a given T
//  the Monte Carlo solves the transport problem, scattering included, in a single pass.
//
//  Same data and same discretization as IQ4.cpp, so that the two can be compared.
//  compile with  g++ -O3 -std=c++17 -o MCIQ4 MCIQ4.cpp
//  run with      ./MCIQ4          or   ./MCIQ4 N=20000 iter=20 nuplot=0.75

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <string>
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
const int    Nc=Nz-1;                      // cells for the Monte Carlo tallies
const double dz=Z/(Nz-1);
const int    jmaxmax=600;
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
long   Nphot  = 20000;      // photons per frequency and per iteration
int    kmax   = 20;         // outer iterations on the temperature
double nuplot = 0.68;      // frequency at which I(z,mu) and Q(z,mu) are kept
int    Nmu    = 40;         // angular bins for the plot
unsigned long seed0 = 20260908UL;
int    doplot = 1;
long   Nplot  = 0;          // photons at nu=nuplot on the last iteration (default 40*N)

// ---------------------------------------------------------------- unknowns
double nu[jmaxmax], kappanu[jmaxmax];
double kapc[jmaxmax][Nc], asc[jmaxmax][Nc];        // kappa and a_s in each cell
double J0[jmaxmax][Nc], J2[jmaxmax][Nc], K0[jmaxmax][Nc], K2[jmaxmax][Nc];
double Tc[Nc], zc[Nc], c0nu[jmaxmax];              // temperature, cell centres, c_0 of (4.13)
double tauN[Nc+1];                                 // optical depth of the cell edges
int    jmax, jplot;
double psiI[Nc*128], psiQ[Nc*128], exI[128], exQ[128], dnI[128], dnQ[128];  // for the plot

mt19937_64 gen;
inline double rnd(){ return (gen()>>11)*0x1.0p-53; }

// ---------------------------------------------------------------- scattering of (4.8)
// K(p<-p', mu<-mu'): probability density of the state and direction after a collision
inline double Kscat(int p,int pp,double mu,double mup){
    const double b=3*_beta/8, u=(1-_beta)/4;
    if(pp==0) return (p==0)? b*(2*(1-mu*mu)*(1-mup*mup)+mu*mu*mup*mup)+u : b*mup*mup+u;
    else      return (p==0)? b*mu*mu+u : b+u;
}
inline double sample1mu2(double U){ return 2*cos((acos(1-2*U)+4*pi)/3);}   // ~ (1-mu^2)
inline double samplemu2 (double U,double V){ double m=cbrt(U); return (V<0.5)?-m:m;} // ~ mu^2
// intensity scattered into the direction mu, by state: c_p = K_{p l} w_l + K_{p r} w_r
inline void cvec(double wl,double wr,double mu,double mup,double& cl,double& cr){
    cl=Kscat(0,0,mu,mup)*wl+Kscat(0,1,mu,mup)*wr;
    cr=Kscat(1,0,mu,mup)*wl+Kscat(1,1,mu,mup)*wr;
}
void scatter(double& wl,double& wr,double& mu){
    const double b=3*_beta/8, u=(1-_beta)/4, mup=mu, W=wl+wr; (void)u;
    // c_l+c_r = A(1-mu^2) + B mu^2 + C, of total mass W: sample mu from it
    double A=2*b*(1-mup*mup)*wl, B=b*mup*mup*wl+b*wr;   // the rest of the mass W is
    double V=rnd()*W;                                   // the isotropic part, 2*C

    if     (V<4*A/3)         mu=sample1mu2(rnd());
    else if(V<4*A/3+2*B/3)   mu=samplemu2(rnd(),rnd());
    else                     mu=2*rnd()-1;
    double cl,cr; cvec(wl,wr,mu,mup,cl,cr);    // deterministic polarization update
    wl=W*cl/(cl+cr); wr=W-wl;
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

// ---------------------------------------------------------------- one flight
// moves the photon by an optical depth tau, scores the track lengths in the cells
// crossed; returns the cell of the collision, or -1 (ground) or -2 (top of atmosphere)
int flight(double& z,double mu,double tau,double wl,double wr,bool score,
           const double* kap,double* tj0,double* tj2,double* tk0,double* tk2,
           double* psi,double* psq){
    const double w=wl+wr, s=(wl-wr)/(w>0?w:1), mu2=mu*mu, im=1/fabs(mu);
    int i=int(z/dz); if(i<0)i=0; if(i>=Nc)i=Nc-1;
    int imu=int((mu+1)*Nmu/2); if(imu<0)imu=0; if(imu>=Nmu)imu=Nmu-1;
    double need=tau*fabs(mu);                  // remaining optical depth, counted on z
    while(true){
        double zb=(mu>0)? (i+1)*dz : i*dz, du=kap[i]*fabs(zb-z);
        bool stop=(du>=need);
        double zn= stop? z+(mu>0? need/kap[i] : -need/kap[i]) : zb;
        if(score){
            double len=fabs(zn-z)*im, wl=w*len;
            tj0[i]+=wl; tj2[i]+=wl*mu2; tk0[i]+=s*wl; tk2[i]+=s*wl*mu2;
            if(psi){ psi[i*Nmu+imu]+=wl; psq[i*Nmu+imu]+=s*wl; }
        }
        z=zn;
        if(stop) return i;
        need-=du;
        i+=(mu>0)?1:-1;
        if(i<0){ z=0; return -1;}
        if(i>=Nc){ z=Z; return -2;}
    }
}

// ---------------------------------------------------------------- Monte Carlo, one frequency
void mcOneFrequency(int j,long N,double& bIn,double& bTop,double& bGnd,double& bAbs){
    const double nuj=nu[j], Bs=Bsun(nuj), Be=Bearth(nuj), rho=-q0/2;  // ground reflectance
    const double *kap=kapc[j], *asz=asc[j];
    const bool   plot=(j==jplot);
    double *tj0=J0[j], *tj2=J2[j], *tk0=K0[j], *tk2=K2[j];
    for(int i=0;i<Nc;i++) tj0[i]=tj2[i]=tk0[i]=tk2[i]=0;
    tauN[0]=0; for(int i=0;i<Nc;i++) tauN[i+1]=tauN[i]+kap[i]*dz;

    // ---- strength of the 3 sources, in units of a current (weight per unit area)
    double Wb=mus*Cs*Bs, Wg=Ce*Be/3, Wv[Nc], Wc[Nc+1];
    Wc[0]=0;
    for(int i=0;i<Nc;i++){ Wv[i]=2*kap[i]*(1-asz[i])*BB(nuj,Tc[i])*dz; Wc[i+1]=Wc[i]+Wv[i];}
    double Wtot=Wb+Wg+Wc[Nc];
    bIn=Wtot;
    if(Wtot<=0) return;
    const double w0=Wtot/N, wmin=1e-4*w0;
    double fdown=0;                            // downward current on the ground

    for(long n=0;n<N;n++){
        double z,mu,wl=w0/2,wr=w0/2, U=rnd()*Wtot;   // the 3 sources are unpolarized
        bool beam=false;
        if(U<Wb){ z=Z; mu=-mus; beam=true; }                    // collimated sun
        else if(U<Wb+Wg){ z=0; mu=cbrt(rnd()); }                // Lambert ground, pdf 3mu^2
        else {                                                  // volume emission
            double V=U-Wb-Wg; int lo=0,hi=Nc;
            while(hi-lo>1){ int m=(lo+hi)/2; if(Wc[m]<=V) lo=m; else hi=m;}
            z=(lo+rnd())*dz; mu=2*rnd()-1;
        }
        for(int k=0;k<10000;k++){
            if(fabs(mu)<1e-9) mu=(mu<0?-1:1)*1e-9;
            int i=flight(z,mu,-log(1-rnd()),wl,wr,!beam,kap,tj0,tj2,tk0,tk2,
                          plot?psiI:0,plot?psiQ:0);
            beam=false;
            double w=wl+wr;
            if(i==-2){ bTop+=w; break; }                        // escapes to space
            if(i==-1){                                          // reaches the ground
                fdown+=w;
                if(rnd()>rho){ bGnd+=w; break; }                // absorbed by the ground
                mu=sqrt(rnd()); wl=wr=w/2;                      // Lambert reflection, Q=0
                continue;
            }
            bAbs+=w*(1-asz[i]);  wl*=asz[i]; wr*=asz[i]; w*=asz[i];   // implicit capture
            if(plot){       // next event estimator of the intensities on the boundaries
                for(int jm=0;jm<Nmu;jm++){
                    double m=-1+(jm+0.5)*2./Nmu, tz=tauN[i]+kap[i]*(z-i*dz);
                    double cl,cr; cvec(wl,wr,m,mu,cl,cr);
                    if(m>0){
                        double e=exp(-(tauN[Nc]-tz)/m)/m;
                        exI[jm]+=(cl+cr)*e; exQ[jm]+=(cl-cr)*e;
                    } else {
                        double e=exp(tz/m)/(-m);
                        dnI[jm]+=(cl+cr)*e; dnQ[jm]+=(cl-cr)*e;
                    }
                }
            }
            scatter(wl,wr,mu);
            if(w<wmin){ if(rnd()>0.25) break; wl*=4; wr*=4; }   // russian roulette
        }
    }
    // ---- normalization, J_q = (1/2) int mu^q I dmu,  K_q = (1/2) int mu^q Q dmu
    for(int i=0;i<Nc;i++){
        tj0[i]/=2*dz; tj2[i]/=2*dz; tk0[i]/=2*dz; tk2[i]/=2*dz;
        double e=exp(-(tauN[Nc]-tauN[i]-kap[i]*dz/2)/mus);      // uncollided sunlight
        tj0[i]+=Cs*Bs*e/2;  tj2[i]+=sqr(mus)*Cs*Bs*e/2;
    }
    c0nu[j]=-q0*fdown;                                          // the c_0 of (4.13)
    if(plot)     // I = sum of the track lengths / (dz dmu),  same for Q with its sign
        for(int i=0;i<Nc*Nmu;i++){ psiI[i]*=Nmu/(2*dz); psiQ[i]*=Nmu/(2*dz); }
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
        double Ta=Tc[i]>0?Tc[i]:Te/2, Tb=2*Ta; int c=0;         // dichotomy
        while(rootT(rhs,Ta,i)>0 && c++<100) Ta/=2;
        while(rootT(rhs,Tb,i)<0 && c++<100) Tb*=2;
        while(Tb-Ta>1e-4 && c++<200){
            double Tm=(Ta+Tb)/2; if(rootT(rhs,Tm,i)>0) Tb=Tm; else Ta=Tm;
        }
        Tc[i]=(Ta+Tb)/2;
        for(int it=0;it<50;it++){                               // Newton
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
        if     (k=="N")      Nphot=(long)v;
        else if(k=="iter")   kmax=(int)v;
        else if(k=="nuplot") nuplot=v;
        else if(k=="Nmu")    Nmu=(int)v;
        else if(k=="seed")   seed0=(unsigned long)v;
        else if(k=="plot")   doplot=(int)v;
        else if(k=="Nplot")  Nplot=(long)v;
        else cerr<<"unknown parameter "<<k<<endl;
    }
    if(Nmu>128) Nmu=128;
    if(Nplot<=0) Nplot=40*Nphot;
    for(int i=0;i<Nc;i++){ zc[i]=(i+0.5)*dz; Tc[i]=Te/2; }
    jmax=readkappa("_kappa.txt");
    jplot=0; for(int j=1;j<jmax;j++) if(fabs(nu[j]-nuplot)<fabs(nu[jplot]-nuplot)) jplot=j;
    cout<<"Z="<<Z<<"  beta="<<_beta<<"  mu_s="<<mus<<"  q0="<<q0<<"  c'_e="<<Ce
        <<"  c_s="<<Cs<<"\n"<<Nphot<<" photons per frequency and per iteration, "
        <<Nc<<" cells;  plots at nu="<<nu[jplot]<<" (lambda="<<3/nu[jplot]<<" microns)\n";

    double t0=clock(), bTop,bGnd,bAbs;
    cout<<"\n iter\t T(ground)\t T(top)\t\t balance\n";
    for(int k=0;k<kmax;k++){
        bTop=bGnd=bAbs=0; double win=0;
        for(int i=0;i<Nc*Nmu;i++) psiI[i]=psiQ[i]=0;
        for(int jm=0;jm<Nmu;jm++) exI[jm]=exQ[jm]=dnI[jm]=dnQ[jm]=0;
        for(int j=0;j<jmax;j++){
            gen.seed(seed0+j);            // correlated sampling: same noise at each iteration
            double in;                    // the last sweep is refined at nu=nuplot,
            long np=(j==jplot && k==kmax-1)? Nplot : Nphot;   // Q is a small difference
            mcOneFrequency(j,np,in,bTop,bGnd,bAbs);
            win+=in;
        }
        genT();
        cout<<" "<<k<<"\t "<<Tc[0]*T0-273<<"\t "<<Tc[Nc-1]*T0-273
            <<"\t "<<(bTop+bGnd+bAbs)/win<<endl;
    }
    printf(" Time CPU = %10.3f s\n",(clock()-t0)/CLOCKS_PER_SEC);

    // -------------------------------------------------- temperature and mean intensity
    ofstream fT("MCIQ4_T.txt");
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

    ofstream fJ("MCIQ4_Jnu.txt");     // the 4 moments at the frequency of the plots
    fJ<<"# z   J0   J2   K0   K2   weight of the Dirac of the sunlight   at nu="
      <<nu[jplot]<<", kappa_nu="<<kappanu[jplot]<<"\n";
    for(int i=0;i<Nc;i++)
        fJ<<zc[i]<<" "<<J0[jplot][i]<<" "<<J2[jplot][i]<<" "<<K0[jplot][i]<<" "
          <<K2[jplot][i]<<" "<<2*J0[jplot][i]*0+Cs*Bsun(nu[jplot])
            *exp(-(tauN[Nc]-tauN[i]-kapc[jplot][i]*dz/2)/mus)<<"\n";
    fJ.close();

    // -------------------------------------------------- I(z,mu) and Q(z,mu) at nu=nuplot
    const double dmu=2./Nmu;
    ofstream fI("MCIQ4_I.txt");
    fI<<"# z   mu   I(z,mu)   Q(z,mu)   at nu="<<nu[jplot]<<"\n";
    for(int i=0;i<Nc;i++){
        for(int jm=0;jm<Nmu;jm++)
            fI<<zc[i]<<" "<<-1+(jm+0.5)*dmu<<" "<<psiI[i*Nmu+jm]<<" "<<psiQ[i*Nmu+jm]<<"\n";
        fI<<"\n";
    }
    fI.close();
    ofstream fE("MCIQ4_exit.txt");
    fE<<"# mu   I   Q   -Q/I   (mu>0: leaving at z=Z; mu<0: arriving on the ground)\n";
    for(int jm=0;jm<Nmu;jm++){
        double m=-1+(jm+0.5)*dmu, a=(m>0)?exI[jm]:dnI[jm], b=(m>0)?exQ[jm]:dnQ[jm];
        fE<<m<<" "<<a<<" "<<b<<" "<<((a>0)?-b/a:0)<<"\n";
    }
    fE.close();
    // -------------------------------------------------- the figures
    ofstream fG("MCIQ4.gp");
    fG<<"# gnuplot MCIQ4.gp\n"
        "set terminal pngcairo size 1350,900 enhanced font 'Helvetica,13'\n"
        "set output 'MCIQ4.png'\n"
        "set palette defined (0 '#2b1b6e', 0.35 '#1f8ac0', 0.6 '#7fd06a', 0.8 '#f6d746', 1 '#b3202a')\n"
        "mus="<<mus<<"\n"
        "set multiplot layout 2,2 title \"Monte Carlo for (4.9)-(4.11):  {/Symbol b}="<<_beta
      <<", {/Symbol m}_s="<<mus<<", q_0="<<q0<<", c'_e="<<Ce<<", c_s="<<Cs<<";   surfaces at "
        "{/Symbol n}="<<nu[jplot]<<" ({/Symbol l}="<<3/nu[jplot]<<" {/Symbol m}m, {/Symbol k}_"
        "{/Symbol n}="<<kappanu[jplot]<<") with "<<Nplot<<" photons\"\n"
        "  set xlabel 'z'; set ylabel '{/Symbol m}'; set xrange [0:"<<Z<<"]; set yrange [-1:1]\n"
        "  set pm3d depthorder; set hidden3d; set style fill transparent solid 0.9\n"
        "  set ticslevel 0.05; set view 55,340,1,1.15; unset key; unset colorbox\n"
        "  set title 'diffuse intensity I(z,{/Symbol m}); red: the Dirac of the sunlight'\n"
        "  splot 'MCIQ4_I.txt' using 1:2:3 with pm3d notitle, \\\n"
        "        'MCIQ4_I.txt' using 1:2:3 with lines lc 'black' lw 0.2 nogrid nocontours notitle, \\\n"
        "        'MCIQ4_Jnu.txt' every 2 using 1:(-mus):6 with impulses lw 1.5 lc 'red' \\\n"
        "             nogrid nocontours notitle\n"
        "  set title 'polarization Q(z,{/Symbol m})'\n"
        "  splot 'MCIQ4_I.txt' using 1:2:4 with pm3d notitle, \\\n"
        "        'MCIQ4_I.txt' using 1:2:4 with lines lc 'black' lw 0.2 nogrid nocontours notitle\n"
        "  unset pm3d; unset hidden3d; set view map; set grid; set key top right\n"
        "  set title 'temperature, solution of (4.10)'\n"
        "  set xlabel 'z'; set ylabel 'T  [Celsius]'; set autoscale x; set autoscale y\n"
        "  plot 'MCIQ4_T.txt' using 1:2 with lines lw 2 lc 'blue' title 'T(z)'\n"
        "  set title 'degree of polarization of the light leaving the atmosphere'\n"
        "  set xlabel '{/Symbol m}'; set ylabel '-Q(Z,{/Symbol m})/I(Z,{/Symbol m})'\n"
        "  set xrange [0:1]\n"
        "  plot 'MCIQ4_exit.txt' using 1:4 with linespoints pt 7 ps 0.6 lw 2 lc 'red' \\\n"
        "       title 'at {/Symbol n}="<<nu[jplot]<<"'\n"
        "unset multiplot\n";
    fG.close();
    if(doplot){
        if(system("gnuplot MCIQ4.gp")) cout<<" gnuplot not found: run it later on MCIQ4.gp"<<endl;
        else cout<<" figure written on MCIQ4.png"<<endl;
    }
    cout<<"\n c_0 of (4.13) at nu="<<nu[jplot]<<" : "<<c0nu[jplot]<<endl;
    cout<<" files MCIQ4_T.txt, MCIQ4_Jnu.txt, MCIQ4_I.txt, MCIQ4_exit.txt, MCIQ4.gp written"<<endl;
    return 0;
}
