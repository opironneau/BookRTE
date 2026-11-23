//  Monte Carlo for the VRTE with a variable refractive index: Chapters 5 and 6.
//
//  Stratified, polarized, coupled with the temperature. With n=n(z) the unknowns are
//  the invariants of Chapter 5,  Itilde = I/n^2,  Qtilde = Q/n^2,  which satisfy the
//  system (4.9),(4.10),(4.11) of Chapter 4 with the transport operator of Chapter 6,
//
//      L(I) = mu d_z I + d_z(log n) (1-mu^2) d_mu I + kappa I ,
//
//  and the Fresnel conditions wherever n jumps. A Monte Carlo needs no discretization
//  of the extra term: along a ray
//
//      n(z) sqrt(1-mu^2) = constant          (Snell),
//
//  so the walk follows the exact characteristic. The index is taken piecewise constant
//  on the cells, so that ONE mechanism covers both chapters: inside a cell the ray is
//  straight, at every interface Snell and Fresnel are applied. A graded n(z) is a
//  sequence of small interfaces (the spurious reflection is O((dn/n)^2) each, and the
//  bending converges to the exact ray); a jump at z=Y is one big interface; a ray that
//  turns back inside a gradient is a total internal reflection at one of them.
//
//  A weight is a current of Itilde. The two weights (w_l,w_r) of MCIQ4.cpp are the
//  components parallel and perpendicular to the meridian plane, i.e. exactly the p and
//  s polarizations of Fresnel, so the interface conditions are DIAGONAL: no matrix,
//  two scalars R_p,R_s. As for scattering, the branch (reflected or transmitted) is
//  sampled and the pair of weights is then updated deterministically, which conserves
//  the total weight exactly and leaves no sampling noise on Q.
//  Nothing else changes: scattering kernel, implicit capture, sources, Lambert ground,
//  temperature loop are those of MCIQ4.cpp.
//
//  Data of IQn.cpp (Chapter 5) so that the two can be compared.
//  compile with  g++ -O3 -std=c++17 -o MCIQn MCIQn.cpp
//  run with      ./MCIQn                         (the jump of IQn.cpp: n=1 then 0.9)
//                ./MCIQn graded=1                (the same jump smoothed over the slab)
//                ./MCIQn nm=1 np=1               (no refraction: must give MCIQ4)

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <string>
using namespace std;
#define sqr(x) ((x)*(x))
const double pi = 4*atan(1.);

const bool CO2=false;
// ---------------------------------------------------------------- data of IQn.cpp
const double Ce=2.4, Cs=2.0e-6;
const double drho=-0.7;
const double Z=1.2;
const double T0=4798;
const double Te=(273+18)/T0, Ts=5778/T0;
const double q0=-0.3;
const double mus=0.5;
const double _beta=0.5;
const int    Nz=200;                        // nodes in z, as in IQn.cpp
const int    Nc=Nz-1;                      // cells
const double dz=Z/(Nz-1);
const int    jmaxmax=600;
const double z1=0.8, z2=0.9, z3=0.8;
const double nu1=0.5, nu2=1;
const double kappamin=0.001;

double nm=1.3, np=1, Y=0.5*Z;             // index below and above Y
int    graded=0;                           // 1: the same variation, spread over (0,Z)

const double kapparhoz(double z){ return (1+drho*z)*(1*(z<Y) + 1*(z>=Y)); }
const double as(double z,double nu){
    return 0.3 + 0.3*4*(z2-z)*(z-z1)*(z>z1)*(z<z2)/sqr(z1-z2)
    + 0.3*16*(nu<nu2)*(nu>nu1)*sqr((nu-nu1)*(nu-nu2)/sqr(nu1-nu2))*(z>z3);
}
double nindex(double z){                   // the refractive index
    if(!graded) return nm+(np-nm)*(z>Y);
    return nm+(np-nm)*z/Z;                 // graded: the same variation, linear in z
}
const double BB(const double nu,const double T){ return sqr(nu)*nu/(exp(fmin(nu/T,100))-1);}
const double dBB(const double nu,const double T){
    double a=exp(nu/T);
    if(a>1e100) return sqr(nu*nu/T)*1e-30;
    return a*sqr(nu*nu/fmax(1.e-30,a-1)/T);
}
const double Bsun(const double nu){ return BB(nu,Ts);}
const double Bearth(const double nu){ return BB(nu,Te);}

// ---------------------------------------------------------------- run time parameters
long   Nphot  = 20000;
int    kmax   = 45;
double nuplot = 0.68;
int    Nmu    = 40;
long   Nplot  = 0;
unsigned long seed0 = 20261019UL;
int    doplot = 1;

// ---------------------------------------------------------------- unknowns
double nu[jmaxmax], kappanu[jmaxmax];
double kapc[jmaxmax][Nc], asc[jmaxmax][Nc];
double J0[jmaxmax][Nc], J2[jmaxmax][Nc], K0[jmaxmax][Nc], K2[jmaxmax][Nc];
double Tc[Nc], zc[Nc], nc[Nc], c0nu[jmaxmax];
double psiI[Nc*128], psiQ[Nc*128], exI[128], exQ[128], dnI[128], dnQ[128];
int    jmax, jplot;
mt19937_64 gen;
inline double rnd(){ return (gen()>>11)*0x1.0p-53; }

// ---------------------------------------------------------------- scattering of (4.8)
inline double Kscat(int p,int pp,double mu,double mup){
    const double b=3*_beta/8, u=(1-_beta)/4;
    if(pp==0) return (p==0)? b*(2*(1-mu*mu)*(1-mup*mup)+mu*mu*mup*mup)+u : b*mup*mup+u;
    else      return (p==0)? b*mu*mu+u : b+u;
}
inline double sample1mu2(double U){ return 2*cos((acos(1-2*U)+4*pi)/3);}
inline double samplemu2 (double U,double V){ double m=cbrt(U); return (V<0.5)?-m:m;}
inline void cvec(double wl,double wr,double mu,double mup,double& cl,double& cr){
    cl=Kscat(0,0,mu,mup)*wl+Kscat(0,1,mu,mup)*wr;
    cr=Kscat(1,0,mu,mup)*wl+Kscat(1,1,mu,mup)*wr;
}
void scatter(double& wl,double& wr,double& mu){
    const double b=3*_beta/8, u=(1-_beta)/4, mup=mu, W=wl+wr;
    double A=2*b*(1-mup*mup)*wl, B=b*mup*mup*wl+b*wr;
    double V=rnd()*W;
    if     (V<4*A/3)       mu=sample1mu2(rnd());
    else if(V<4*A/3+2*B/3) mu=samplemu2(rnd(),rnd());
    else                   mu=2*rnd()-1;
    double cl,cr; cvec(wl,wr,mu,mup,cl,cr);
    wl=W*cl/(cl+cr); wr=W-wl; (void)u;
}

// ---------------------------------------------------------------- Fresnel
// medium 1 -> medium 2, incidence cosine c1=|mu1|. Returns 1 if total reflection,
// else 0 and the refracted cosine c2 and the two reflectances (p for I_l, s for I_r).
inline int fresnel(double n1,double n2,double c1,double& c2,double& Rp,double& Rs){
    double s2=sqr(n1/n2)*(1-c1*c1);
    if(s2>=1){ c2=0; Rp=Rs=1; return 1; }                 // total internal reflection
    c2=sqrt(1-s2);
    Rs=sqr((n1*c1-n2*c2)/(n1*c1+n2*c2));
    Rp=sqr((n2*c1-n1*c2)/(n2*c1+n1*c2));
    return 0;
}

//*******************************
int readkappa(string mykappafile){
//------------------------------
    ifstream kappafile(mykappafile);
    if(!kappafile) throw runtime_error("Cannot open file: "+mykappafile);
    cout<<"reading kappa(nu) from file "<<mykappafile<<endl;
    int j=-1; double kappaux,waveno;
    while((j++<=jmaxmax)&&(!kappafile.eof())){
        kappafile >> waveno >> kappaux;
        kappanu[j]=fmax(kappaux,kappamin);
        nu[j]=3/waveno;
        for(int i=0;i<Nc;i++){
            kapc[j][i]=kappanu[j]*kapparhoz((i+0.5)*dz);
            asc[j][i] =as((i+0.5)*dz,nu[j]);
        }
    }
    kappafile.close();
    cout<<"Number of frequencies "<<j<<endl;
    return j;
}

// ---------------------------------------------------------------- one free flight
// Advances the photon by the optical depth tau along the exact ray: straight inside a
// cell, Snell and Fresnel at each interface. Scores the track lengths if score.
// Returns the cell of the collision, -1 on the ground, -2 at the top of the atmosphere.
int flight(double& z,double& mu,double& wl,double& wr,double tau,bool& score,
           const double* kap,double* tj0,double* tj2,double* tk0,double* tk2,
           double* psi,double* psq){
    int i=int(z/dz); if(i<0)i=0; if(i>=Nc)i=Nc-1;
    double need=tau;
    while(true){
        if(fabs(mu)<1e-12) mu=(mu<0?-1:1)*1e-12;
        double zb=(mu>0)? (i+1)*dz : i*dz;
        double path=fabs(zb-z)/fabs(mu), du=kap[i]*path;
        bool stop=(du>=need);
        double len=stop? need/kap[i] : path;               // geometric path in the cell
        double zn=stop? z+mu*len : zb;
        if(score){
            double w=wl+wr, s=(wl-wr)/(w>0?w:1), wl2=w*len, m2=mu*mu;
            tj0[i]+=wl2; tj2[i]+=wl2*m2; tk0[i]+=s*wl2; tk2[i]+=s*wl2*m2;
            if(psi){
                int im=int((mu+1)*Nmu/2); if(im<0)im=0; if(im>=Nmu)im=Nmu-1;
                psi[i*Nmu+im]+=wl2; psq[i*Nmu+im]+=s*wl2;
            }
        }
        z=zn;
        if(stop) return i;
        need-=du;
        int inext=(mu>0)? i+1 : i-1;
        if(inext<0){ z=0; return -1;}
        if(inext>=Nc){ z=Z; return -2;}
        // ---- interface between the cells i and inext
        double c2,Rp,Rs, W=wl+wr;
        if(fresnel(nc[i],nc[inext],fabs(mu),c2,Rp,Rs)){ mu=-mu; continue; } // total refl.
        double WR=Rp*wl+Rs*wr;
        if(rnd()*W<WR){                                    // reflected
            wl=W*Rp*wl/WR; wr=W-wl; mu=-mu;
            score=true;             // a reflected beam is no longer the direct sunlight
        } else {                                           // transmitted
            double WT=W-WR;
            wl=W*(1-Rp)*wl/WT; wr=W-wl;
            mu=(mu>0)? c2 : -c2;  i=inext;
        }
    }
}

// ---------------------------------------------------------------- Monte Carlo, one frequency
void mcOneFrequency(int j,long N,double& bIn,double& bTop,double& bGnd,double& bAbs){
    const double nuj=nu[j], Bs=Bsun(nuj), Be=Bearth(nuj), rho=-q0/2;
    const double *kap=kapc[j], *asz=asc[j];
    const bool   plot=(j==jplot);
    double *tj0=J0[j], *tj2=J2[j], *tk0=K0[j], *tk2=K2[j];
    for(int i=0;i<Nc;i++) tj0[i]=tj2[i]=tk0[i]=tk2[i]=0;

    // ---- the uncollided sunlight: it bends, it is attenuated, and Fresnel polarizes it.
    // Along the ray n^2(1-mu^2) is constant and the physical flux n^2 mu A is conserved,
    // so the amplitude of the Dirac is A = A(Z) n_Z^2 mu_s /(n^2 |mu_b|) x T x exp(-tau).
    double mub[Nc],Al[Nc],Ar[Nc];
    {   double C=sqr(nc[Nc-1])*(1-sqr(mus)), tl=Cs*Bs/2, tr=Cs*Bs/2, tauz=0;
        for(int i=Nc-1;i>=0;i--){
            double m2=1-C/sqr(nc[i]);
            if(m2<=0){ for(int k=i;k>=0;k--){ mub[k]=0; Al[k]=Ar[k]=0;} break;}
            mub[i]=-sqrt(m2);
            if(i<Nc-1){                                    // interface i+1 -> i
                double c2,Rp,Rs;
                if(fresnel(nc[i+1],nc[i],fabs(mub[i+1]),c2,Rp,Rs)){ tl=tr=0; }
                else { tl*=(1-Rp); tr*=(1-Rs); }
            }
            double e=exp(-(tauz+kap[i]*dz/2/fabs(mub[i])));
            double fac=sqr(nc[Nc-1])*mus/(sqr(nc[i])*fabs(mub[i]));
            Al[i]=tl*fac*e; Ar[i]=tr*fac*e;
            tauz+=kap[i]*dz/fabs(mub[i]);
        }
    }
    // ---- next event estimator: the ray that leaves at mu_exit, seen from each cell
    static double Eup[Nc][128],Mup[Nc][128],TPu[Nc][128],TSu[Nc][128];
    static double Edn[Nc][128],Mdn[Nc][128],TPd[Nc][128],TSd[Nc][128];
    if(plot){
        const double dmu=2./Nmu;
        for(int jm=0;jm<Nmu;jm++){
            double me=-1+(jm+0.5)*dmu;
            for(int i=0;i<Nc;i++){ Eup[i][jm]=Edn[i][jm]=0; Mup[i][jm]=Mdn[i][jm]=0;
                                   TPu[i][jm]=TSu[i][jm]=TPd[i][jm]=TSd[i][jm]=0;}
            if(me>0){                                      // upwards, towards space
                double C=sqr(nc[Nc-1])*(1-me*me), tp=1,ts=1,tauz=0;
                for(int i=Nc-1;i>=0;i--){
                    double m2=1-C/sqr(nc[i]); if(m2<=0) break;
                    double m=sqrt(m2);
                    if(i<Nc-1){
                        double c2,Rp,Rs;
                        if(fresnel(nc[i],nc[i+1],m,c2,Rp,Rs)) break;
                        tp*=(1-Rp); ts*=(1-Rs);
                    }
                    Mup[i][jm]=m; TPu[i][jm]=tp; TSu[i][jm]=ts;
                    Eup[i][jm]=exp(-tauz)*sqr(nc[Nc-1])*me/(sqr(sqr(nc[i]))*m2);
                    tauz+=kap[i]*dz/m;
                }
            } else {                                       // downwards, onto the ground
                double C=sqr(nc[0])*(1-me*me), tp=1,ts=1,tauz=0;
                for(int i=0;i<Nc;i++){
                    double m2=1-C/sqr(nc[i]); if(m2<=0) break;
                    double m=sqrt(m2);
                    if(i>0){
                        double c2,Rp,Rs;
                        if(fresnel(nc[i],nc[i-1],m,c2,Rp,Rs)) break;
                        tp*=(1-Rp); ts*=(1-Rs);
                    }
                    Mdn[i][jm]=m; TPd[i][jm]=tp; TSd[i][jm]=ts;
                    Edn[i][jm]=exp(-tauz)*sqr(nc[0])*(-me)/(sqr(sqr(nc[i]))*m2);
                    tauz+=kap[i]*dz/m;
                }
            }
        }
    }
    // ---- the three sources, in units of a current
    double Wb=sqr(nc[Nc-1])*mus*Cs*Bs, Wg=sqr(nc[0])*Ce*Be/3, Wv[Nc], Wc[Nc+1];
    Wc[0]=0;
    for(int i=0;i<Nc;i++){ Wv[i]=2*sqr(nc[i])*kap[i]*(1-asz[i])*BB(nuj,Tc[i])*dz;
                           Wc[i+1]=Wc[i]+Wv[i];}
    double Wtot=Wb+Wg+Wc[Nc];
    bIn=Wtot;
    if(Wtot<=0) return;
    const double w0=Wtot/N, wmin=1e-4*w0;
    double fdown=0;

    for(long n=0;n<N;n++){
        double z,mu,wl=w0/2,wr=w0/2, U=rnd()*Wtot;
        bool score=true;
        if(U<Wb){ z=Z-1e-12; mu=-mus; score=false; }       // collimated sun
        else if(U<Wb+Wg){ z=0; mu=cbrt(rnd()); }           // Lambert ground
        else { double V=U-Wb-Wg; int lo=0,hi=Nc;
               while(hi-lo>1){ int m=(lo+hi)/2; if(Wc[m]<=V) lo=m; else hi=m;}
               z=(lo+rnd())*dz; mu=2*rnd()-1; }
        for(int k=0;k<10000;k++){
            int i=flight(z,mu,wl,wr,-log(1-rnd()),score,kap,tj0,tj2,tk0,tk2,
                         plot?psiI:0,plot?psiQ:0);
            score=true;
            double w=wl+wr;
            if(i==-2){ bTop+=w; break; }
            if(i==-1){ fdown+=w;
                if(rnd()>rho){ bGnd+=w; break;}
                mu=sqrt(rnd()); wl=wr=w/2; continue;
            }
            bAbs+=w*(1-asz[i]); wl*=asz[i]; wr*=asz[i]; w*=asz[i];
            if(plot){
                for(int jm=0;jm<Nmu;jm++){
                    double me=-1+(jm+0.5)*2./Nmu;
                    if(me>0){
                        if(Eup[i][jm]<=0) continue;
                        double m=Mup[i][jm], cl,cr; cvec(wl,wr,m,mu,cl,cr);
                        double e=Eup[i][jm]*exp(-kap[i]*((i+1)*dz-z)/m);
                        cl*=TPu[i][jm]; cr*=TSu[i][jm];   // Fresnel is diagonal in (l,r)
                        exI[jm]+=(cl+cr)*e; exQ[jm]+=(cl-cr)*e;
                    } else {
                        if(Edn[i][jm]<=0) continue;
                        double m=Mdn[i][jm], cl,cr; cvec(wl,wr,-m,mu,cl,cr);
                        double e=Edn[i][jm]*exp(-kap[i]*(z-i*dz)/m);
                        cl*=TPd[i][jm]; cr*=TSd[i][jm];   // Fresnel is diagonal in (l,r)
                        dnI[jm]+=(cl+cr)*e; dnQ[jm]+=(cl-cr)*e;
                    }
                }
            }
            scatter(wl,wr,mu);
            if(w<wmin){ if(rnd()>0.25) break; wl*=4; wr*=4;}
        }
    }
    for(int i=0;i<Nc;i++){           // J_q,K_q of Itilde, + the uncollided sunlight
        double n2=2*dz*sqr(nc[i]);            // Itilde = I/n^2
        tj0[i]/=n2; tj2[i]/=n2; tk0[i]/=n2; tk2[i]/=n2;
        tj0[i]+=(Al[i]+Ar[i])/2;          tj2[i]+=sqr(mub[i])*(Al[i]+Ar[i])/2;
        tk0[i]+=(Al[i]-Ar[i])/2;          tk2[i]+=sqr(mub[i])*(Al[i]-Ar[i])/2;
    }
    c0nu[j]=-q0*fdown/sqr(nc[0]);
    if(plot) for(int i=0;i<Nc;i++) for(int jm=0;jm<Nmu;jm++){
        double f=Nmu/(2*dz*sqr(nc[i])); psiI[i*Nmu+jm]*=f; psiQ[i*Nmu+jm]*=f;}
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
                double dnu=nu[j]-nu[j-1], nmid=(nu[j]+nu[j-1])/2,
                       ka=(1-as(zc[i],nu[j]))*kappanu[j];
                left+=ka*BB(nmid,Tc[i])*dnu; deriv+=ka*dBB(nmid,Tc[i])*dnu;
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
        else if(k=="Nplot")  Nplot=(long)v;
        else if(k=="seed")   seed0=(unsigned long)v;
        else if(k=="plot")   doplot=(int)v;
        else if(k=="nm")     nm=v;
        else if(k=="np")     np=v;
        else if(k=="Y")      Y=v;
        else if(k=="graded") graded=(int)v;
        else cerr<<"unknown parameter "<<k<<endl;
    }
    if(Nmu>128) Nmu=128;
    if(Nplot<=0) Nplot=40*Nphot;
    for(int i=0;i<Nc;i++){ zc[i]=(i+0.5)*dz; Tc[i]=Te/2; nc[i]=nindex(zc[i]); }
    jmax=readkappa("_kappa.txt");
    if(CO2)
        for(int j=0;j<jmax;j++) // CO2
            if((nu[j]>3./18)&&(nu[j]<3./14))
                kappanu[j]=fmax(0.5,1.5*kappanu[j]);//Visible R

    jplot=0; for(int j=1;j<jmax;j++) if(fabs(nu[j]-nuplot)<fabs(nu[jplot]-nuplot)) jplot=j;
    cout<<"Z="<<Z<<"  n: "<<nm<<" -> "<<np<<(graded?" (graded)":" (jump at Y)")
        <<"  Y="<<Y<<"  mu_s="<<mus<<"  q0="<<q0<<"\n"<<Nphot
        <<" photons per frequency and per iteration, "<<Nc<<" cells\n";
    if(!graded){
        double c2,Rp,Rs;
        int tir=fresnel(nm,np,1.,c2,Rp,Rs);
        cout<<"  at normal incidence R="<<(Rp+Rs)/2<<"; upward total reflection for |mu|<"
            <<((np<nm)? sqrt(1-sqr(np/nm)) : 0.)<<(tir?" ":" ")<<endl;
    }

    double t0=clock(), bTop,bGnd,bAbs;
    cout<<"\n iter\t T(ground)\t T(top)\t\t balance\n";
    for(int k=0;k<kmax;k++){
        bTop=bGnd=bAbs=0; double win=0;
        for(int i=0;i<Nc*Nmu;i++) psiI[i]=psiQ[i]=0;
        for(int jm=0;jm<Nmu;jm++) exI[jm]=exQ[jm]=dnI[jm]=dnQ[jm]=0;
        for(int j=0;j<jmax;j++){
            gen.seed(seed0+j);
            double in; long np2=(j==jplot && k==kmax-1)? Nplot : Nphot;
            mcOneFrequency(j,np2,in,bTop,bGnd,bAbs);
            win+=in;
        }
        genT();
        cout<<" "<<k<<"\t "<<Tc[0]*T0-273<<"\t "<<Tc[Nc-1]*T0-273
            <<"\t "<<(bTop+bGnd+bAbs)/win<<endl;
    }
    printf(" Time CPU = %10.3f s\n",(clock()-t0)/CLOCKS_PER_SEC);

    ofstream fT("MCIQn_T.txt");
    fT<<"# z   T[Celsius]   int J0 dnu   int K0 dnu   n(z)\n";
    cout<<"\n z\t\t T[C]\t\t int J0 dnu\t int K0 dnu\n";
    for(int i=0;i<Nc;i++){
        double sJ=0,sK=0;
        for(int j=1;j<jmax;j++){
            sJ+=(J0[j][i]+J0[j-1][i])/2*(nu[j]-nu[j-1]);
            sK+=(K0[j][i]+K0[j-1][i])/2*(nu[j]-nu[j-1]);
        }
        fT<<zc[i]<<" "<<Tc[i]*T0-273<<" "<<sJ<<" "<<sK<<" "<<nc[i]<<"\n";
        if(i%8==0) cout<<zc[i]<<"\t "<<Tc[i]*T0-273<<"\t "<<sJ<<"\t "<<sK<<endl;
    }
    fT.close();

    const double dmu=2./Nmu;
    ofstream fI("MCIQn_I.txt");
    fI<<"# z   mu   Itilde(z,mu)   Qtilde(z,mu)   at nu="<<nu[jplot]<<"\n";
    for(int i=0;i<Nc;i++){
        for(int jm=0;jm<Nmu;jm++)
            fI<<zc[i]<<" "<<-1+(jm+0.5)*dmu<<" "<<psiI[i*Nmu+jm]<<" "<<psiQ[i*Nmu+jm]<<"\n";
        fI<<"\n";
    }
    fI.close();
    ofstream fE("MCIQn_exit.txt");
    fE<<"# mu   I   Q   -Q/I   (mu>0: leaving at z=Z; mu<0: arriving on the ground)\n";
    for(int jm=0;jm<Nmu;jm++){
        double m=-1+(jm+0.5)*dmu, a=(m>0)?exI[jm]:dnI[jm], b=(m>0)?exQ[jm]:dnQ[jm];
        fE<<m<<" "<<a<<" "<<b<<" "<<((a>0)?-b/a:0)<<"\n";
    }
    fE.close();
    cout<<"\n files MCIQn_T.txt, MCIQn_I.txt, MCIQn_exit.txt written"<<endl;
    return 0;
}
