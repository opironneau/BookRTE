// Monte Carlo solution of Milne's problem
//
//      mu dI/dz + kappa I = (sigma_s/2) int_{-1}^{1} I(z,mu') dmu',   z in (0,Z), mu in (-1,1)
//      I(0,mu) = delta(mu-mu_s)  for mu>0,      I(Z,mu) = 0  for mu<0.
//
// Photons are emitted at z=0 in the direction mu_s and followed until they leave the slab.
// Absorption is treated by weight reduction (implicit capture), with albedo a = sigma_s/kappa;
// scattering is isotropic because the right hand side is (sigma_s/2)*int I dmu'.
//
// Tallies
//   - track length estimator of  phi(z) = int_{-1}^1 I dmu = 2 J_0(z), and of I(z,mu);
//   - next event estimator of the emergent intensities I(Z,mu), mu>0 and I(0,mu), mu<0;
//   - the uncollided (direct) beam is not sampled, it is added analytically:
//         I_dir(z,mu) = delta(mu-mu_s) exp(-kappa z/mu_s),  phi_dir(z) = exp(-kappa z/mu_s).
//   Error bars come from batch means.
//
// compile with  g++ -O3 -std=c++17 -o MilneMC MilneMC.cpp
// run with      ./MilneMC              or  ./MilneMC kappa=1 sigmas=0.9 N=4000000

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <random>
#include <vector>
#include <string>
#include <time.h>
using namespace std;

// ---------------------------------------------------------------- parameters
double Z      = 1.;             // thickness of the slab
double kappa  = 1.;             // extinction coefficient
double sigmas = 1.;             // scattering coefficient (kappa>=sigmas)
double mus    = 1./sqrt(3.);    // cosine of the direction of the incoming beam
long   Nhist  = 2000000;        // number of photons
int    Nz     = 50;             // cells in z
int    Nmu    = 40;             // bins in mu
int    Nbatch = 20;             // batches, for the error bars
unsigned long seed = 20260827UL;
int    doplot = 1;              // 1 = write milneMC_I.gp and call gnuplot for the surface

double dz, dmu, albedo;
const double wmin = 1e-4;       // russian roulette threshold (relative to w0)
const int    kmax = 100000;     // safety bound on the number of scatterings

mt19937_64 gen;
inline double rnd(){ return (gen()>>11)*0x1.0p-53; }   // uniform in [0,1)

// track length scoring of one free flight from z0 to zend, direction mu, weight w
inline void scoreTrack(double z0,double zend,double mu,double w,
                       vector<double>& phi, vector<double>& psi){
    int imu = int((mu+1.)/dmu); if(imu<0) imu=0; if(imu>=Nmu) imu=Nmu-1;
    const double invmu = 1./fabs(mu);
    double za = min(z0,zend), zb = max(z0,zend);
    int ka = int(za/dz), kb = int(zb/dz);
    if(ka<0) ka=0;  if(kb>=Nz) kb=Nz-1;
    for(int k=ka; k<=kb; k++){
        double lo = max(za,k*dz), hi = min(zb,(k+1)*dz);
        if(hi<=lo) continue;
        double len = (hi-lo)*invmu;          // path length (3d) inside cell k
        phi[k]            += w*len/dz;       // int I dmu
        psi[k*Nmu + imu]  += w*len/(dz*dmu); // I(z,mu)
    }
}

int main(int argc, const char* argv[]){
    double t0 = clock();
    for(int i=1;i<argc;i++){                 // command line: key=value
        string s(argv[i]); size_t p=s.find('=');
        if(p==string::npos){ cerr<<"skipping "<<s<<endl; continue; }
        string k=s.substr(0,p); double v=atof(s.c_str()+p+1);
        if     (k=="Z")      Z=v;
        else if(k=="kappa")  kappa=v;
        else if(k=="sigmas") sigmas=v;
        else if(k=="mus")    mus=v;
        else if(k=="N")      Nhist=(long)v;
        else if(k=="Nz")     Nz=(int)v;
        else if(k=="Nmu")    Nmu=(int)v;
        else if(k=="batches")Nbatch=(int)v;
        else if(k=="seed")   seed=(unsigned long)v;
        else if(k=="plot")   doplot=(int)v;
        else cerr<<"unknown parameter "<<k<<endl;
    }
    if(sigmas>kappa){ cerr<<"need sigma_s <= kappa"<<endl; return 1; }
    if(mus<=0 || mus>1){ cerr<<"need 0 < mu_s <= 1"<<endl; return 1; }
    dz = Z/Nz;  dmu = 2./Nmu;  albedo = sigmas/kappa;
    gen.seed(seed);
    const long nper = max(1L, Nhist/Nbatch);   // photons per batch
    Nhist = nper*Nbatch;

    // batch accumulators and their sums, for mean and standard error
    vector<double> phi(Nz), psi(Nz*Nmu), up(Nmu), dn(Nmu);
    vector<double> Sphi(Nz,0), Qphi(Nz,0), Spsi(Nz*Nmu,0), Sup(Nmu,0), Qup(Nmu,0),
                   Sdn(Nmu,0), Qdn(Nmu,0);
    double Sesc[2]={0,0}, Qesc[2]={0,0}, Sabs=0;   // escaping currents, absorbed weight

    for(int ib=0; ib<Nbatch; ib++){
        fill(phi.begin(),phi.end(),0.); fill(psi.begin(),psi.end(),0.);
        fill(up.begin(),up.end(),0.);   fill(dn.begin(),dn.end(),0.);
        double escT=0, escB=0, abso=0;
        // the incoming current is int_0^1 mu I(0,mu) dmu = mu_s, shared by nper photons
        const double w0 = mus/nper;

        for(long n=0; n<nper; n++){
            double z=0., mu=mus, w=w0;
            bool first=true;                    // 1st flight = direct beam, done analytically
            for(int k=0;k<kmax;k++){
                double s = -log(1.-rnd())/kappa;        // free path
                double zend = z + mu*s;
                if(zend>=Z || zend<=0){                 // the photon leaves the slab
                    double zb = (zend>=Z)? Z : 0.;
                    if(!first) scoreTrack(z,zb,mu,w,phi,psi);
                    if(zend>=Z) escT += w; else escB += w;
                    break;
                }
                if(!first) scoreTrack(z,zend,mu,w,phi,psi);
                first=false;
                z = zend;
                abso += w*(1.-albedo);
                w   *= albedo;                          // implicit capture
                if(w < wmin*w0){                        // russian roulette
                    if(rnd()>0.25) break;
                    w *= 4.;
                }
                // next event estimator of the emergent intensities: the scattered source
                // at z is (sigma_s/2)phi(z), it reaches the boundary along mu unscattered
                for(int j=0;j<Nmu;j++){
                    double m = -1.+(j+0.5)*dmu;
                    if(m>0) up[j] += 0.5*w/m*exp(-kappa*(Z-z)/m);
                    else    dn[j] += 0.5*w/(-m)*exp(-kappa*z/(-m));
                }
                mu = 2.*rnd()-1.;                       // isotropic scattering
                if(fabs(mu)<1e-12) mu = (mu<0?-1:1)*1e-12;
            }
        }
        for(int k=0;k<Nz;k++){ Sphi[k]+=phi[k]; Qphi[k]+=phi[k]*phi[k]; }
        for(int k=0;k<Nz*Nmu;k++) Spsi[k]+=psi[k];
        for(int j=0;j<Nmu;j++){ Sup[j]+=up[j]; Qup[j]+=up[j]*up[j];
                                Sdn[j]+=dn[j]; Qdn[j]+=dn[j]*dn[j]; }
        Sesc[0]+=escT; Qesc[0]+=escT*escT; Sesc[1]+=escB; Qesc[1]+=escB*escB; Sabs+=abso;
    }

    // mean and standard error over the batches
    auto stat=[&](double S,double Q,double& m,double& e){
        m = S/Nbatch;
        double v = (Q - Nbatch*m*m)/(Nbatch-1.);
        e = (v>0)? sqrt(v/Nbatch) : 0.;
    };
    printf( " Time CPU = %10.6f\n", (clock() - t0)/CLOCKS_PER_SEC);
    cout<<scientific<<setprecision(4);
    cout<<"# Milne problem by Monte Carlo:  Z="<<Z<<"  kappa="<<kappa<<"  sigma_s="<<sigmas
        <<"  albedo="<<albedo<<"  mu_s="<<mus<<"\n#  "<<Nhist<<" photons in "<<Nbatch
        <<" batches, "<<Nz<<" cells, "<<Nmu<<" angles\n";

    // ---------------------------------------------------- J_0 = (1/2) int I dmu
    ofstream fJ("milneMC_J0.txt");
    fJ<<"# z   J0(z)   error   J0_direct(z)   a*J0 = I(z,0+-)   exp(-kappa z/mu_s)\n";
    cout<<"\n z\t\tJ0(z)\t\terror\t\tJ0 direct\n";
    for(int k=0;k<Nz;k++){
        double z=(k+0.5)*dz, m,e; stat(Sphi[k],Qphi[k],m,e);
        double dir=exp(-kappa*z/mus);                  // uncollided contribution to int I dmu
        double J0=0.5*(m+dir);
        fJ<<z<<" "<<J0<<" "<<0.5*e<<" "<<0.5*dir<<" "<<albedo*J0<<" "
          <<exp(-kappa*z/mus)<<"\n";
        if(k%max(1,Nz/10)==0) cout<<z<<"\t"<<0.5*(m+dir)<<"\t"<<0.5*e<<"\t"<<0.5*dir<<endl;
    }
    fJ.close();

    // ------------------------------------- the surface I(z,mu): file + gnuplot script
    // column 3 is the diffuse (collided) intensity, column 4 adds the direct beam
    // delta(mu-mu_s)exp(-kappa z/mu_s) averaged over the mu bin which contains mu_s
    int jbeam = int((mus+1.)/dmu); if(jbeam>=Nmu) jbeam=Nmu-1;
    ofstream fI("milneMC_I.txt");
    fI<<"# z   mu   I_diffuse(z,mu)   I_diffuse+direct beam smeared on its mu bin\n";
    for(int k=0;k<Nz;k++){
        double z=(k+0.5)*dz, beam=exp(-kappa*z/mus)/dmu;
        for(int j=0;j<Nmu;j++){
            double v=Spsi[k*Nmu+j]/Nbatch;
            fI<<z<<" "<<-1.+(j+0.5)*dmu<<" "<<v<<" "<<v+(j==jbeam)*beam<<"\n";
        }
        fI<<"\n";                       // gnuplot needs a blank line between the scans
    }
    fI.close();

    ofstream fG("milneMC_I.gp");         // surface plot, column 3 (diffuse) by default
    fG<<"# gnuplot milneMC_I.gp\n"
        "# col=3: diffuse intensity, the direct beam being drawn apart as a Dirac;\n"
        "# col=4: the beam is instead smeared over the mu bin which contains mu_s.\n"
        "col = 3\n"
        "mus = "<<mus<<"\n"
        "set terminal pngcairo size 1500,760 enhanced font 'Helvetica,13'\n"
        "set output 'milneMC_I.png'\n"
        "set palette defined (0 '#2b1b6e', 0.35 '#1f8ac0', 0.6 '#7fd06a', 0.8 '#f6d746', 1 '#b3202a')\n"
        "set xlabel 'z'; set ylabel '{/Symbol m}'\n"
        "set xrange [0:"<<Z<<"]; set yrange [-1:1]\n"
        "set multiplot layout 1,2 title \"Milne problem by Monte Carlo: diffuse intensity "
        "I(z,{/Symbol m}),   Z="<<Z<<", {/Symbol k}="<<kappa<<", {/Symbol s}_s="<<sigmas
      <<", {/Symbol m}_s="<<mus<<",  "<<Nhist<<" photons\"\n"
        "  set pm3d depthorder\n"
        "  set hidden3d; set style fill transparent solid 0.9\n"
        "  set contour base; set cntrparam levels 12; unset clabel\n"
        "  set zlabel 'I' rotate by 90\n"
        "  set view 55,340,1,1.15\n"
        "  set ticslevel 0.05\n"
        "  unset colorbox; unset key\n"
        "  set zrange [0:1.05]\n"
        "  set label 11 'Dirac {/Symbol d}({/Symbol m}-{/Symbol m}_s), of weight "
        "e^{-{/Symbol k}z/{/Symbol m}_s}' at 0.05*"<<Z<<",mus,1.0 tc 'red' front\n"
        "  set label 12 'crest = {/Symbol s}_sJ_0(z)/{/Symbol k} : the grazing limit "
        "I(z,0^{/*0.8 +}_{/*0.8 -})' at 0.45*"<<Z<<",0,0.80 tc 'black' front\n"
        "  splot 'milneMC_I.txt' using 1:2:col with pm3d notitle, \\\n"
        "        'milneMC_I.txt' using 1:2:col with lines lc 'black' lw 0.2 nogrid notitle, \\\n"
        "        'milneMC_J0.txt' using 1:(0):5 with lines lw 3 lc 'black' nogrid nocontours notitle, \\\n"
        "        'milneMC_J0.txt' every 2 using 1:(mus):6 with impulses lw 1.5 lc 'red' nogrid nocontours notitle\n"
        "  set colorbox; set view map; set pm3d interpolate 2,2\n"
        "  unset label 11; unset label 12; unset zrange\n"
        "  set arrow 1 from 0,mus,9e9 to "<<Z<<",mus,9e9 nohead lc 'red' lw 2 dt 2 front\n"
        "  set label 1 '  {/Symbol d}({/Symbol m}-{/Symbol m}_s), the incoming beam' at 0,mus,9e9 "
        "front tc 'red' offset 0,0.7\n"
        "  set arrow 2 from 0,0,9e9 to "<<Z<<",0,9e9 nohead lc 'black' lw 1 dt 3 front\n"
        "  set label 2 '  {/Symbol m}=0: I = {/Symbol s}_sJ_0/{/Symbol k}' at 0.55*"<<Z
      <<",0,9e9 front offset 0,0.7\n"
        "  set title 'seen from above, with the level lines'\n"
        "  splot 'milneMC_I.txt' using 1:2:col with pm3d notitle, \\\n"
        "        'milneMC_I.txt' using 1:2:col with lines lc 'black' lw 0.6 nosurface notitle\n"
        "unset multiplot\n";
    fG.close();

    if(doplot){
        if(system("gnuplot milneMC_I.gp"))
            cout<<" gnuplot not found: run it later on milneMC_I.gp"<<endl;
        else cout<<" surface plot written on milneMC_I.png"<<endl;
    }

    // ---------------------------------------------------- emergent intensities
    ofstream fE("milneMC_exit.txt");
    fE<<"# mu   I(Z,mu) mu>0   error     I(0,mu) mu<0   error\n";
    cout<<"\n mu\t\tI(Z,mu),mu>0\terror\t\tI(0,mu),mu<0\terror\n";
    double curT=0, curB=0;
    for(int j=0;j<Nmu;j++){
        double m=-1.+(j+0.5)*dmu, a,ea,b,eb;
        stat(Sup[j],Qup[j],a,ea);  stat(Sdn[j],Qdn[j],b,eb);
        if(m>0){ fE<<m<<" "<<a<<" "<<ea<<" 0 0\n"; curT += m*a*dmu; }
        else   { fE<<m<<" 0 0 "<<b<<" "<<eb<<"\n"; curB += (-m)*b*dmu; }
        if(j%max(1,Nmu/10)==0) cout<<m<<"\t"<<a<<"\t"<<ea<<"\t"<<b<<"\t"<<eb<<endl;
    }
    fE.close();

    // ---------------------------------------------------- balance
    double eT,eB,mT,mB; stat(Sesc[0],Qesc[0],mT,eT); stat(Sesc[1],Qesc[1],mB,eB);
    double dirT = mus*exp(-kappa*Z/mus);   // current of the transmitted direct beam
    cout<<"\n incoming current             mu_s          = "<<mus<<endl;
    cout<<" escaping current at z=Z (analog)           = "<<mT<<" +/- "<<eT<<endl;
    cout<<"                         (next event+beam)  = "<<curT+dirT<<endl;
    cout<<" escaping current at z=0 (analog)           = "<<mB<<" +/- "<<eB<<endl;
    cout<<"                         (next event)       = "<<curB<<endl;
    cout<<" absorbed                                   = "<<Sabs/Nbatch<<endl;
    cout<<" balance (out+absorbed)/in                  = "<<(mT+mB+Sabs/Nbatch)/mus<<endl;
    cout<<"\n files milneMC_J0.txt, milneMC_I.txt, milneMC_I.gp, milneMC_exit.txt written"<<endl;
    return 0;
}
