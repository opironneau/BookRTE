//  Discrete ordinates with Fresnel conditions, for the two-medium case of Chapter 5:
//  n = nm below Y, n = np above. Written to verify MCIQn.cpp; same data, same cells,
//  same temperature loop, independent of it in every other respect.
//
//  The difficulty of a discrete ordinates method here is that Snell does not map a
//  quadrature node to a quadrature node, so IQn.cpp must interpolate. It is avoided by
//  labelling a ray by its invariant s = n sqrt(1-mu^2), which is the same on both sides
//  of Y: a node stays a node. Two families of rays:
//
//   - K transmitted rays: a Gauss rule mu_p in (0,1) ABOVE Y; below Y the same rays have
//     mu_m = sqrt(1-(np/nm)^2(1-mu_p^2)), which fills (mu_c,1) with mu_c=sqrt(1-(np/nm)^2);
//     their quadrature weight below Y is w_k np^2 mu_p /(nm^2 mu_m), the Jacobian of Snell;
//   - J trapped rays: a Gauss rule on (0,mu_c), which exist only below Y. They are the
//     rays totally reflected at Y, the ones that warm the lower layer.
//
//  Together they integrate exactly over mu in (0,1) on each side. At Y the invariant
//  Itilde = I/n^2 is simply multiplied by the Fresnel R or T = 1-R, diagonal in (l,r),
//  the reflected flux being taken from the previous iterate: the Fresnel coupling, the
//  Lambert ground and the scattering all converge in the same fixed point.
//
//  compile with  g++ -O3 -std=c++17 -o SNIQn SNIQn.cpp
//  run with      ./SNIQn               or   ./SNIQn K=24 J=12 iter=25 nm=1 np=0.9

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <string>
#include <ctime>
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
const int    Nz=80, Nc=Nz-1;
const double dz=Z/(Nz-1);
const int    jmaxmax=555, Kmax=128;
const double z1=0.8, z2=0.9, z3=0.8;
const double nu1=0.5, nu2=1;
const double kappamin=0.001;
double nm=1.3, np=1., Y=0.5*Z;

const double kapparhoz(double z){ return (1+drho*z)*(1*(z<Y) + 1*(z>=Y)); }
const double as(double z,double nu){
    return 0.3 + 0.3*4*(z2-z)*(z-z1)*(z>z1)*(z<z2)/sqr(z1-z2)
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

int    K=24, J=12, kmax=15, itmax=800;
int    KCO2=0;   // 1: the CO2 case of IQn.cpp, kappa_nu -> max(0.5,1.5 kappa_nu)
double tol=1e-9;
double nu[jmaxmax], kappanu[jmaxmax];
double kapc[jmaxmax][Nc], asc[jmaxmax][Nc];
double J0[jmaxmax][Nc], J2[jmaxmax][Nc], K0[jmaxmax][Nc], K2[jmaxmax][Nc];
double Tc[Nc], zc[Nc], nc[Nc];
int    jmax, iY;                            // iY: first cell above Y
double mub_[Kmax],mum_[Kmax],wq_up[Kmax],wq_dn[Kmax];   // transmitted rays
double mut_[Kmax],wt_[Kmax];                            // trapped rays
double Rp_,Rs_,Tp_,Ts_;                                 // Fresnel at Y, per ray below
double RpK[Kmax],RsK[Kmax],TpK[Kmax],TsK[Kmax];
long   nsweep=0;

void gauleg(int n,double a,double b,double* x,double* w){  // on (a,b)
    for(int i=0;i<(n+1)/2;i++){
        double z=cos(pi*(i+0.75)/(n+0.5)), z1g,pp;
        do{ double p1=1,p2=0;
            for(int j=0;j<n;j++){ double p3=p2; p2=p1; p1=((2*j+1)*z*p2-j*p3)/(j+1);}
            pp=n*(z*p1-p2)/(z*z-1); z1g=z; z=z1g-p1/pp;
        } while(fabs(z-z1g)>1e-14);
        double xm=(a+b)/2, xl=(b-a)/2, ww=xl*2/((1-z*z)*pp*pp);
        x[i]=xm-xl*z;     w[i]=ww;
        x[n-1-i]=xm+xl*z; w[n-1-i]=ww;
    }
}
inline int fresnel(double n1,double n2,double c1,double& c2,double& Rp,double& Rs){
    double s2=sqr(n1/n2)*(1-c1*c1);
    if(s2>=1){ c2=0; Rp=Rs=1; return 1;}
    c2=sqrt(1-s2);
    Rs=sqr((n1*c1-n2*c2)/(n1*c1+n2*c2));
    Rp=sqr((n2*c1-n1*c2)/(n2*c1+n1*c2));
    return 0;
}
int readkappa(string f){
    ifstream kf(f); if(!kf) throw runtime_error("Cannot open file: "+f);
    cout<<"reading kappa(nu) from file "<<f<<endl;
    int j=-1; double ka,wn;
    while((j++<=jmaxmax)&&(!kf.eof())){
        kf >> wn >> ka; kappanu[j]=fmax(ka,kappamin); nu[j]=3/wn;
        if(KCO2==1 && nu[j]>3./18 && nu[j]<3./14)      // same rule as IQn.cpp
            kappanu[j]=fmax(0.5,1.5*kappanu[j]);
        for(int i=0;i<Nc;i++){
            kapc[j][i]=kappanu[j]*kapparhoz((i+0.5)*dz);
            asc[j][i]=as((i+0.5)*dz,nu[j]);
        }
    }
    kf.close(); cout<<"Number of frequencies "<<j<<endl; return j;
}

// ---------------------------------------------------------------- one frequency
void snOneFrequency(int j){
    const double nuj=nu[j], Bs=Bsun(nuj), Be=Bearth(nuj);
    const double *kap=kapc[j], *asz=asc[j];
    const double b=3*_beta/8, u=(1-_beta)/4;
    const int NR=K+J;
    double sigs[Nc],siga[Nc],B[Nc],beamL[Nc],beamR[Nc],mubeam[Nc];
    for(int i=0;i<Nc;i++){ sigs[i]=kap[i]*asz[i]; siga[i]=kap[i]*(1-asz[i]); B[i]=BB(nuj,Tc[i]);}
    // ---- uncollided sunlight: Snell, Fresnel and the Jacobian of the Dirac
    {   double C=sqr(nc[Nc-1])*(1-sqr(mus)), tl=Cs*Bs/2, tr=Cs*Bs/2, tauz=0;
        for(int i=Nc-1;i>=0;i--){
            double m2=1-C/sqr(nc[i]);
            if(m2<=0){ for(int k=i;k>=0;k--){mubeam[k]=0;beamL[k]=beamR[k]=0;} break;}
            mubeam[i]=-sqrt(m2);
            if(i<Nc-1 && nc[i]!=nc[i+1]){
                double c2,Rp,Rs;
                if(fresnel(nc[i+1],nc[i],fabs(mubeam[i+1]),c2,Rp,Rs)) tl=tr=0;
                else { tl*=(1-Rp); tr*=(1-Rs);}
            }
            double e=exp(-(tauz+kap[i]*dz/2/fabs(mubeam[i])));
            double fac=sqr(nc[Nc-1])*mus/(sqr(nc[i])*fabs(mubeam[i]));
            beamL[i]=tl*fac*e; beamR[i]=tr*fac*e;
            tauz+=kap[i]*dz/fabs(mubeam[i]);
        }
    }
    // ---- moments, iterated
    static double m0l[Nc],m0r[Nc],m2l[Nc],m2r[Nc];
    static double RdnL[Kmax],RdnR[Kmax],RupL[Kmax],RupR[Kmax],TrL[Kmax],TrR[Kmax];
    for(int i=0;i<Nc;i++) m0l[i]=m0r[i]=m2l[i]=m2r[i]=0;
    for(int k=0;k<K;k++) RdnL[k]=RdnR[k]=RupL[k]=RupR[k]=0;
    for(int k=0;k<J;k++) TrL[k]=TrR[k]=0;              // trapped, upward value at Y

    for(int it=0;it<itmax;it++){
        nsweep++;
        double nl[Nc],nr[Nc],n2l[Nc],n2r[Nc];
        for(int i=0;i<Nc;i++) nl[i]=nr[i]=n2l[i]=n2r[i]=0;
        double Sa[Nc],Sb[Nc],Sr[Nc];
        for(int i=0;i<Nc;i++){
            double M0l=m0l[i]+beamL[i], M0r=m0r[i]+beamR[i],
                   M2l=m2l[i]+sqr(mubeam[i])*beamL[i];
            Sa[i]=sigs[i]*(2*b*(M0l-M2l)+u*(M0l+M0r))+siga[i]*B[i]/2;
            Sb[i]=sigs[i]*b*(M2l+M0r-2*(M0l-M2l));
            Sr[i]=sigs[i]*(b*(M2l+M0r)+u*(M0l+M0r))+siga[i]*B[i]/2;
        }
        double fdown=0;
        double dnL[Kmax+Kmax],dnR[Kmax+Kmax];              // arriving on the ground
        double dnYL[Kmax],dnYR[Kmax];                      // downward, just ABOVE Y
        // ---------------- downward sweeps
        for(int r=0;r<NR;r++){
            bool trapped=(r>=K);
            int k=trapped? r-K : r;
            double il,ir;
            int istart;
            if(!trapped){                                   // from the top of the slab
                il=ir=0; istart=Nc-1;
                for(int i=Nc-1;i>=iY;i--){
                    double m=mub_[k], e=exp(-kap[i]*dz/m), f=m/(kap[i]*dz)*(1-e);
                    double sl=(Sa[i]+Sb[i]*sqr(m))/kap[i], sr=Sr[i]/kap[i];
                    nl[i]+=wq_up[k]*(sl+(il-sl)*f);  nr[i]+=wq_up[k]*(sr+(ir-sr)*f);
                    n2l[i]+=wq_up[k]*sqr(m)*(sl+(il-sl)*f); n2r[i]+=wq_up[k]*sqr(m)*(sr+(ir-sr)*f);
                    il=sl+(il-sl)*e; ir=sr+(ir-sr)*e;
                }
                dnYL[k]=il; dnYR[k]=ir;                     // keep the value above Y
                il=TpK[k]*il+RpK[k]*RupL[k];                // Fresnel at Y, from both sides
                ir=TsK[k]*ir+RsK[k]*RupR[k];
                istart=iY-1;
            } else { il=TrL[k]; ir=TrR[k]; istart=iY-1; }   // totally reflected at Y
            for(int i=istart;i>=0;i--){
                double m=trapped? mut_[k] : mum_[k], w=trapped? wt_[k] : wq_dn[k];
                double e=exp(-kap[i]*dz/m), f=m/(kap[i]*dz)*(1-e);
                double sl=(Sa[i]+Sb[i]*sqr(m))/kap[i], sr=Sr[i]/kap[i];
                nl[i]+=w*(sl+(il-sl)*f); nr[i]+=w*(sr+(ir-sr)*f);
                n2l[i]+=w*sqr(m)*(sl+(il-sl)*f); n2r[i]+=w*sqr(m)*(sr+(ir-sr)*f);
                il=sl+(il-sl)*e; ir=sr+(ir-sr)*e;
            }
            dnL[r]=il; dnR[r]=ir;
            fdown+=(trapped? wt_[k]*mut_[k] : wq_dn[k]*mum_[k])*(il+ir);
        }
        // the uncollided sunlight also lands on the ground
        if(mubeam[0]!=0)
            fdown+=fabs(mubeam[0])*(beamL[0]+beamR[0])*exp(-kap[0]*dz/2/fabs(mubeam[0]));
        const double c0=-q0*fdown;
        // ---------------- upward sweeps
        for(int r=0;r<NR;r++){
            bool trapped=(r>=K);
            int k=trapped? r-K : r;
            double m0=trapped? mut_[k] : mum_[k];
            double il=Ce*Be*m0/2+c0/2, ir=il;               // Lambert ground
            for(int i=0;i<iY;i++){
                double w=trapped? wt_[k] : wq_dn[k];
                double e=exp(-kap[i]*dz/m0), f=m0/(kap[i]*dz)*(1-e);
                double sl=(Sa[i]+Sb[i]*sqr(m0))/kap[i], sr=Sr[i]/kap[i];
                nl[i]+=w*(sl+(il-sl)*f); nr[i]+=w*(sr+(ir-sr)*f);
                n2l[i]+=w*sqr(m0)*(sl+(il-sl)*f); n2r[i]+=w*sqr(m0)*(sr+(ir-sr)*f);
                il=sl+(il-sl)*e; ir=sr+(ir-sr)*e;
            }
            if(trapped){ TrL[k]=il; TrR[k]=ir; continue; }  // total reflection at Y
            RupL[k]=il; RupR[k]=ir;                         // for the next downward sweep
            double jl=TpK[k]*il+RpK[k]*RdnL[k], jr=TsK[k]*ir+RsK[k]*RdnR[k];
            for(int i=iY;i<Nc;i++){
                double m=mub_[k], e=exp(-kap[i]*dz/m), f=m/(kap[i]*dz)*(1-e);
                double sl=(Sa[i]+Sb[i]*sqr(m))/kap[i], sr=Sr[i]/kap[i];
                nl[i]+=wq_up[k]*(sl+(jl-sl)*f); nr[i]+=wq_up[k]*(sr+(jr-sr)*f);
                n2l[i]+=wq_up[k]*sqr(m)*(sl+(jl-sl)*f); n2r[i]+=wq_up[k]*sqr(m)*(sr+(jr-sr)*f);
                jl=sl+(jl-sl)*e; jr=sr+(jr-sr)*e;
            }
        }
        for(int r=0;r<K;r++){ RdnL[r]=dnYL[r]; RdnR[r]=dnYR[r]; } // downward value at Y+
        double err=0,nrm=1e-300;
        for(int i=0;i<Nc;i++){
            err=fmax(err,fabs(nl[i]+nr[i]-m0l[i]-m0r[i])); nrm=fmax(nrm,fabs(nl[i]+nr[i]));
            m0l[i]=nl[i]; m0r[i]=nr[i]; m2l[i]=n2l[i]; m2r[i]=n2r[i];
        }
        if(err<tol*nrm) break;
    }
    for(int i=0;i<Nc;i++){
        J0[j][i]=(m0l[i]+m0r[i]+beamL[i]+beamR[i])/2;
        J2[j][i]=(m2l[i]+m2r[i]+sqr(mubeam[i])*(beamL[i]+beamR[i]))/2;
        K0[j][i]=(m0l[i]-m0r[i]+beamL[i]-beamR[i])/2;
        K2[j][i]=(m2l[i]-m2r[i]+sqr(mubeam[i])*(beamL[i]-beamR[i]))/2;
    }
}
// ---------------------------------------------------------------- temperature
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
        while(Tb-Ta>1e-4 && c++<200){ double Tm=(Ta+Tb)/2; if(rootT(rhs,Tm,i)>0) Tb=Tm; else Ta=Tm;}
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
        if(k=="K") K=(int)v; else if(k=="KCO2") KCO2=(int)v; else if(k=="J") J=(int)v; else if(k=="iter") kmax=(int)v;
        else if(k=="tol") tol=v; else if(k=="nm") nm=v; else if(k=="np") np=v;
        else if(k=="Y") Y=v; else cerr<<"unknown parameter "<<k<<endl;
    }
    if(K>Kmax) K=Kmax; if(J>Kmax) J=Kmax;
    for(int i=0;i<Nc;i++){ zc[i]=(i+0.5)*dz; Tc[i]=Te/2; nc[i]=nm+(np-nm)*(zc[i]>Y); }
    iY=0; while(iY<Nc && nc[iY]==nc[0]) iY++;              // first cell above Y
    // ---- the two families of rays
    double xp[Kmax],wp[Kmax];
    gauleg(K,0.,1.,xp,wp);                                 // above Y
    double muc=(np<nm)? sqrt(1-sqr(np/nm)) : 0.;
    for(int k=0;k<K;k++){
        mub_[k]=xp[k]; wq_up[k]=wp[k];
        double m2=1-sqr(np/nm)*(1-sqr(xp[k]));             // the same ray below Y
        mum_[k]=sqrt(fmax(m2,1e-12));
        wq_dn[k]=wp[k]*sqr(np)*xp[k]/(sqr(nm)*mum_[k]);    // Jacobian of Snell
        double c2;
        fresnel(nm,np,mum_[k],c2,RpK[k],RsK[k]);
        TpK[k]=1-RpK[k]; TsK[k]=1-RsK[k];
    }
    double xt[Kmax],wt2[Kmax];
    if(muc>0){ gauleg(J,0.,muc,xt,wt2); for(int k=0;k<J;k++){ mut_[k]=xt[k]; wt_[k]=wt2[k];} }
    else J=0;
    double s1=0,s2=0;
    for(int k=0;k<K;k++){ s1+=wq_up[k]; s2+=wq_dn[k]; }
    for(int k=0;k<J;k++) s2+=wt_[k];
    jmax=readkappa("_kappa.txt");
    cout<<"KCO2 = "<<KCO2<<endl;
    if(CO2)
        for(int j=0;j<jmax;j++) // CO2
            if((nu[j]>3./18)&&(nu[j]<3./14))
                kappanu[j]=fmax(0.5,1.5*kappanu[j]);//Visible R
    cout<<"Fresnel discrete ordinates: n="<<nm<<" | "<<np<<" at z="<<iY*dz
        <<",  mu_c="<<muc<<"\n  "<<K<<" transmitted + "<<J<<" trapped rays; the two"
        " quadratures sum to "<<s1<<" (above) and "<<s2<<" (below)\n";

    double t0=clock();
    cout<<"\n iter\t T(ground)\t T(top)\n";
    for(int k=0;k<kmax;k++){
        for(int j=0;j<jmax;j++) snOneFrequency(j);
        genT();
        cout<<" "<<k<<"\t "<<Tc[0]*T0-273<<"\t "<<Tc[Nc-1]*T0-273<<endl;
    }
    printf(" Time CPU = %10.3f s  (%.1f sweeps per frequency and per temperature iteration)\n",
           (clock()-t0)/CLOCKS_PER_SEC,double(nsweep)/jmax/kmax);
    ofstream fT("SNIQn_T.txt");
    fT<<"# z   T[Celsius]   int J0 dnu   int K0 dnu   n(z)\n";
    for(int i=0;i<Nc;i++){
        double sJ=0,sK=0;
        for(int j=1;j<jmax;j++){
            sJ+=(J0[j][i]+J0[j-1][i])/2*(nu[j]-nu[j-1]);
            sK+=(K0[j][i]+K0[j-1][i])/2*(nu[j]-nu[j-1]);
        }
        fT<<zc[i]<<" "<<Tc[i]*T0-273<<" "<<sJ<<" "<<sK<<" "<<nc[i]<<"\n";
    }
    fT.close();
    cout<<" file SNIQn_T.txt written"<<endl;
    return 0;
}
