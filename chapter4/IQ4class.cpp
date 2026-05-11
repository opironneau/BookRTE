//
//  Created by Olivier Pironneau on 02/05/2026.
//  used in chapter 4
//  This program is the same as IQchap4.cpp but C++ classes to improve readibility
//  A factor of 2 in speed is lost because of the arrays are X<vector<vector<double>> instead of X[][]

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <array>
#include "IQ4class.hpp"
//const double pi = std::atan(1.)*4, stefan=pi*pi*pi*pi/15;

namespace physics{
    const double Ce= 2.4, Cs=5.0e-6 ; // intensity of earth and sun radiation
    const double TOA=1.2;// Top of Atmosphere
    const double rho0=1.,drho=-0.7; // fake density and density gradient
    const double B0=1.4744e6, T0=4798; // scaling
    const double Te=(273+18)/T0, Ts=5798/T0;// Earth and Sun temperature in °K
    const double q0=-0.3; // albedo intensity
    const double mus=0.5; // minus direction of collimated light
    const double beta=0.5;  // % of isotropic to Rayleigh scattering
    const double nuRayleighm=0.5,nuRayleighM=1; // frequency range of Rayleigh scattering
    const double nuWaterS = 3./20, nuWaterm = 3./10, nuWaterM = 3./1; // frequency range of water vapor absorption
    const double cloud = 1.5; // strength of cloud in (zcloudm,zcloudM)
    double zcloudm = 0.4, zcloudM = 0.7, zRayleigh = 0.8; // Rayleigh scatt is applied above zRayleigh
    const double CO2m=3./18, CO2M=3./14, dark=0.2; // ghg C02 opaque for nu in (CO2m,CO2M)
}

namespace algo{
    const double Tstart = physics::Te/2; // initial T for ISIF
    const std::string basedir("/Users/pironneau/Dropbox/aranger/TeX2026/BookVrte/prog4/IQ4chap4/IQ4chap4/");
    const std::string kappafile(basedir+"_kappa.txt"); // contains nu and kappa
     const int Nz=80;  // size of altitude grid
    const int kmax=9;  // nb of ISIF iterations
    const bool verbose=false; // more output for debugging
//   Parameters belows need not be changed
    const int newton = 50; // max iter # to compute T from the temperature equation
    const double epsdycho=1e-4, epsnewton=1.e-10;  // precision for dychotomy & Newton
    const double kappamin=0.001;  // if kappa - read - is too small, max it with kappamin
}
using namespace std;

// -------------------------------------------------------
//   Main class to define global arrays
// -------------------------------------------------------
class rte {
public:
    static std::vector<double> T;
    static std::vector<std::vector<double>> J0,J2,K0,K2;
    rte(){ ; }
};

// -------------------------------------------------------
//   Global arrays
// -------------------------------------------------------
std::vector<double> rte::T;
std::vector<std::vector<double>> rte::J0,rte::J2,rte::K0,rte::K2;

// -------------------------------------------------------
//   Defines things related to the altitude grid and optical thickness
// -------------------------------------------------------
class AtmosphereGrid {
public:
    // From physical Altitude h  →  optical coordonnate z:  z = h * rho0 * (1 + drho * h / 2)
    double zprime(double h) const {
        return h * physics::rho0 * (1.0 + physics::drho * h / 2.0);
    }
    // From optical coordinate  z  →  to physical altitude h : solves h*rho0*(1 + drho*h/2) = z
    double backtoz(double z) const {
        if (physics::drho == 0.0)  return z / physics::rho0;
        double a = physics::drho * physics::rho0 / 2.0, b = physics::rho0, c = -z;
        return (-b + std::sqrt(b*b - 4.0*a*c)) / (2.0*a);
    }
    const double Z = zprime(physics::TOA);
    const double dz = Z/(algo::Nz-1);
    const double Nz = algo::Nz;
    const double z1 = zprime(physics::zcloudm);
    const double z2 = zprime(physics::zcloudM);
    const double z3 = zprime(physics::zRayleigh);
};

// -------------------------------------------------------
//   Read kappa(nu) and build kappa(nu,z) and its z-integral. Cloud and GHG affect kappa
// -------------------------------------------------------
class OpacityTable {
public:
    OpacityTable(const AtmosphereGrid& grid,
        const std::string&    filepath): grid_(grid)    { readAndBuild(filepath); }
    // -------------------------------------------------------
    // Access Functions
    // -------------------------------------------------------
    int    jmax()             const { return static_cast<int>(nu_.size()); }
    double nu    (int j)      const { return nu_[j];      }
    double kappanu(int j)     const { return kappanu_[j]; }
    double kappaz (int j, int i) const { return kappaz_[j][i];  }
    double kappasum(int j, int i) const { return kappasum_[j][i]; }

    double kapparhoz(double nu, double z) const {
        double cloudMask = static_cast<double>(z  > grid_.z1)
                        * static_cast<double>(z  < grid_.z2);
        double freqMask  = static_cast<double>(nu > physics::nuWaterm)
                        * static_cast<double>(nu < physics::nuWaterM)
                        + static_cast<double>(nu < physics::nuWaterS);
        return 1.0 + physics::cloud * cloudMask * freqMask;// total kappa is kappa(n) times this
    }

    // -------------------------------------------------------
    // Read j → nu(j) and j  →  kappa(nu(j))=kappa(j), j=1..jmax
    // -------------------------------------------------------
    void readAndBuild(const std::string& filepath) {
        std::ifstream file(filepath);
        if (!file)
            throw std::runtime_error("OpacityTable: cannot open " + filepath);

        std::cout << "Read kappa " << filepath << "\n";

        std::vector<double> nuRaw, kappanuRaw;
        double waveno, kappaux;
        while (file >> waveno >> kappaux) {
            nuRaw    .push_back(3.0 / waveno);
            double ghgMask = static_cast<double>(3.0 / waveno > physics::CO2m) *
            static_cast<double>(3.0 / waveno < physics::CO2M);
            kappanuRaw.push_back(std::max(kappaux*(1+ghgMask*physics::dark), algo::kappamin));
        }
        file.close();

        const int jmax = static_cast<int>(nuRaw.size());
        if (jmax == 0)
            throw std::runtime_error("OpacityTable: kappa file contains nothing");

        std::cout << "Number of frequencies read : " << jmax << "\n";

        // Transfert dans les membres
        nu_      = std::move(nuRaw);
        kappanu_ = std::move(kappanuRaw);

        // Allocation of tables 2D [jmax][Nz]
        kappaz_.assign(jmax, std::vector<double>(algo::Nz, 0.0));
        kappasum_.assign(jmax, std::vector<double>(algo::Nz, 0.0));

        // Precompute kappaz and its z-intégral
        const double dz = grid_.dz;
        for (int j = 0; j < jmax; ++j) {
            kappasum_[j][0] = 0.0;
            for (int i = 0; i < algo::Nz; i++) {
                double z      = i * dz;
                kappaz_[j][i] = kapparhoz(nu_[j], z);
                if (i > 0)
                    kappasum_[j][i] = kappasum_[j][i-1]
                        + (kappaz_[j][i] + kappaz_[j][i-1]) * dz / 2.0;
            }
        }
    }
    // Dépendencies
    const AtmosphereGrid& grid_;

    // Spectral tables  [jmax]
    std::vector<double> nu_;
    std::vector<double> kappanu_;

    // Spatiospectral tables  [jmax][Nz]
    std::vector<std::vector<double>> kappaz_;
    std::vector<std::vector<double>> kappasum_;
};

// -------------------------------------------------------
// Defines J0,K0... and contains the function to update these. Initialize T
// -------------------------------------------------------
class RadiativeField:rte {
public:
    double z1,z2,z3;
    RadiativeField(const AtmosphereGrid& grid, const OpacityTable& opacity)
    : rte(),grid_(grid), opacity_(opacity), jmax_(opacity.jmax())//,z1(grid_.z1),z2(grid_.z2),z3(grid_.z3)
    {
        T.assign(algo::Nz, algo::Tstart);
        J0.assign(jmax_, std::vector<double>(algo::Nz, 0.0));
        J2.assign(jmax_, std::vector<double>(algo::Nz, 0.0));
        K0.assign(jmax_, std::vector<double>(algo::Nz, 0.0));
        K2.assign(jmax_, std::vector<double>(algo::Nz, 0.0));
    }
    
    static inline double sqr(const double x) { return x*x; }
    
    double scattering(double z, double nu) const {
        const double nu1= physics::nuRayleighm,
                     nu2= physics::nuRayleighM ;
        const double cloudMask = static_cast<double>(z  > z1) * static_cast<double>(z  < z2);
        double RayleighMask = static_cast<double>(z  > z3);
        double freqMask  = static_cast<double>(nu > nu1) * static_cast<double>(nu < nu2);
        return 0.3 + 0.3 *physics::cloud * cloudMask*4*(z2-z)*(z-z1)/(z1-z2)/(z1-z2)
                   + 0.3 *freqMask*RayleighMask*16*sqr((nu-nu1)*(nu-nu2)/sqr(nu1-nu2));
     }
    // -------------------------------------------------------
    // Implements ISIF
    // -------------------------------------------------------
    void updateJK()//(const std::vector<double>& T)
    {
         const double mus=physics::mus,
                     beta = physics::beta,
                     Cs = physics::Cs,
                     Ce = physics::Ce,
                     Ts = physics::Ts,
                     Te = physics::Te,
                     q0 = physics::q0;
        const double dz = grid_.dz;

        std::vector<double> H(Nz_), S(Nz_);

        // Main frequency loop
        for (int jnu = 0; jnu < jmax_; ++jnu) {
                 const double kappanuj = opacity_.kappanu(jnu);
                 const double nuj      = opacity_.nu(jnu);

                 // albedo
                 double c0 = -q0 * mus * Cs * Bsun(nuj, Ts)
                            * std::exp(-kappanuj * opacity_.kappasum(jnu, Nz_-1) / mus);

                 // Remplissage H et S
                 for (int i = 0; i < Nz_; ++i) {
                     double z      = i * dz;
                     double kappaz = opacity_.kappaz(jnu, i) * kappanuj;
                     double asz    = scattering(z, nuj);
                     double sigs   = kappaz * asz;
                     double siga   = kappaz * (1.0 - asz);

                     S[i] = sigs * J0[jnu][i] + siga * BB(nuj, T[i]);
                     H[i] = 9.0 * beta * sigs / 8.0
                          * (J2[jnu][i] - J0[jnu][i]/3.0
                             - K0[jnu][i] + K2[jnu][i]);

                     c0 -= q0 * ExpIntegral::E2(opacity_.kappasum(jnu, i) * kappanuj)
                              * J0[jnu][i] * dz;
                 }

                 // Main altitude loop
                 for (int i = 0; i < Nz_; ++i) {
                     double ksz = opacity_.kappasum(jnu, i) * kappanuj;
                     double aux = Cs * Bsun(nuj, Ts)
                                * std::exp(-(opacity_.kappasum(jnu, Nz_-1)
                                            - opacity_.kappasum(jnu, i))
                                            * kappanuj / mus) / 2.0;

                     double J0z = Ce * Bearth(nuj, Te) * ExpIntegral::E3(ksz) / 2.0;
                             J0z+=   + aux
                                + c0 * ExpIntegral::E2(ksz) / 2.0;
                     double J2z = Ce * Bearth(nuj, Te) * ExpIntegral::E5(ksz) / 2.0
                                + sqr(mus) * aux
                                + c0 * ExpIntegral::E4(ksz) / 2.0;
                     double K0z = 0.0, K2z = 0.0;

                     // Convolution en z'
                     for (int j = 1; j < Nz_; ++j) {
                         double Hj  = (H[j] + H[j-1]) / 2.0;
                         double Sj  = (S[j] + S[j-1]) / 2.0 - Hj / 3.0;
                         int ip = (i < j) ? j : i;
                         int jp = (i < j) ? i : j;
                         double arg = kappanuj
                                    * (opacity_.kappasum(jnu, ip)
                                      - opacity_.kappasum(jnu, jp) - dz / 2.0);

                         J0z += (ExpIntegral::E1(arg) * Sj + ExpIntegral::E3(arg) * Hj) * dz / 2.0;
                         J2z += (ExpIntegral::E3(arg) * Sj + ExpIntegral::E5(arg) * Hj) * dz / 2.0;
                         K0z -= (ExpIntegral::E1(arg) - ExpIntegral::E3(arg)) * Hj * dz / 2.0;
                         K2z -= (ExpIntegral::E3(arg) - ExpIntegral::E5(arg)) * Hj * dz / 2.0;
                     }

                     J0[jnu][i] = J0z;
                     J2[jnu][i] = J2z;
                     K0[jnu][i] = K0z;
                     K2[jnu][i] = K2z;
                 }
             }
        }

    //   Planck functions and derivative
         static double Bsun(double nu, double Ts)  { return sqr(nu)*nu / (std::exp(nu/Ts)  - 1.0); }
         static double Bearth(double nu, double Te){ return sqr(nu)*nu / (std::exp(nu/Te)  - 1.0); }
         static double BB(double nu, double T)     { return sqr(nu)*nu / (std::exp(std::fmin(nu/T, 100.0)) - 1.0); }
        static double dBB(double nu, double T)     {
            double a = exp(nu/T);
            if(a>1e100) return sqr(nu*nu/T)*1e-30; // needed but could endanger the Newton iterations
            return a*sqr(nu*nu/fmax(1.e-30,a -1)/T);
    }
         // Dépendencies
         const AtmosphereGrid& grid_;
         const OpacityTable&   opacity_;
         int jmax_, Nz_=algo::Nz;
};

// -------------------------------------------------------
// Manage the Temperature equation
// -------------------------------------------------------
class Tsolver:rte{
   private:
        const AtmosphereGrid& grid_;
        const OpacityTable& opacity_;
        const RadiativeField& field_;
    int jmax_;
        double dz_;
    
   public:
    Tsolver(const AtmosphereGrid& grid, const OpacityTable& opacity, const RadiativeField& field)
    : rte(), grid_(grid), opacity_(opacity), field_(field)
    { jmax_ = opacity_.jmax(); dz_=grid_.dz;  }

   private:
    double root(const double rhs, const double T0,int i) {
        double myeq=-rhs;
        for(int j=1; j<jmax_;j++){
            double nuj = opacity_.nu(j), nuj1 = opacity_.nu(j-1);
            myeq += (1-field_.scattering(i*dz_,nuj))*opacity_.kappanu(j)*field_.BB((nuj+nuj1)/2,T0)*(nuj-nuj1);
        }
        return myeq;
    }
    double rhsTeq(const int i){
        double rhs=0;
        for(int j=1; j<jmax_;j++){
            double nuj = opacity_.nu(j), nuj1 = opacity_.nu(j-1);
            rhs+=(1-field_.scattering(i*dz_,nuj))*opacity_.kappanu(j)*(J0[j][i]+J0[j-1][i])/2*(nuj-nuj1);
        }
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
            if(algo::verbose) cout<<T0<<" down "<<myeq0<<endl;
        }
        while (myeq1<0 && counter<100){
            T1=2*T1; myeq1=root(rhs,T1,i); counter++;
            if(algo::verbose) cout<<T1<<" up "<<myeq1<<endl;
        }
        if(algo::verbose){
            myeq0=root(rhs,T0,i);
            myeq1=root(rhs,T1,i);
            if(myeq0>myeq1) cout<<"BUG in dychotomy\n";
        }
        while (T1-T0 > algo::epsdycho && counter<100){
            Taux=(T1+T0)/2; counter++;
            myeq0=root(rhs,Taux,i);
            if(myeq0>0) T1=Taux; else T0=Taux;
            if(algo::verbose){
                myeq0=root(rhs,T0,i);
                myeq1=root(rhs,T1,i);
                cout<<T0<<" middle "<<T1<<" "<<myeq0 << " "<< myeq1<<endl;
            }
        }
        if(counter>99) cout<<"Divergence in dichotomy\n";
        return (T1+T0)/2;
    }

    public:
    
    int updateT(){ // Solves the temperature equation by dychotomy + Newton
        for(int i=0;i<grid_.Nz;i++){
            double rhs=rhsTeq(i);
            T[i] = getTbydycho(rhs,i,T[i]);
            double presfunc = 1;
            int inewton =0;
            while(inewton++<algo::newton && fabs(presfunc) > algo::epsnewton){
                double T0 = T[i];
                double left=0,  deriv=0;
                for(int j=1; j<jmax_;j++){
                    double dnu=opacity_.nu(j)-opacity_.nu(j-1), nu11=(opacity_.nu(j)+opacity_.nu(j-1))/2,
                     kappaux=(1-field_.scattering(i*dz_,opacity_.nu(j)-opacity_.nu(j-1)))*opacity_.kappanu(j);
                     left += kappaux*field_.BB(nu11,T0)*dnu;
                    deriv += kappaux*field_.dBB(nu11,T0)*dnu;
                }
                presfunc = rhs-left;
                if(fabs(deriv)>1e-10) T[i] = T0+presfunc/deriv;
                if(algo::verbose)
                    cout<<i<<" "<<" T= "<<T[i]<<" residue= "<<presfunc<<" deriv= "<<deriv<<" Tstefan= "<< sqrt(sqrt(rhs*15*2))/3.1416<<endl; // the 2 is kappanuj
            }
            if(inewton>=algo::newton)cout << "Newton precision doubtful "<<endl;
        }
        return 0;
    }
};
    
int main(){
    rte problem;
    AtmosphereGrid zgrid;
    OpacityTable nugrid(zgrid,algo::kappafile);
    RadiativeField rfield(zgrid, nugrid);
    Tsolver solve(zgrid, nugrid, rfield);
    
    cout<<"Iteration  T[2]"<<endl;
    for(int k=0;k<algo::kmax;k++){
        rfield.updateJK();
        solve.updateT();
        cout<<k<<" "<<problem.T[2]*physics::T0-273<<endl;
    }
// Print results
    cout <<endl<<"Altitude z T(z) and int_R^+ J_0(nu,z)dnu"<<endl;
    double meanJ0=0;
    for(int i=0;i<zgrid.Nz;i++){
        for(int j=1;j<nugrid.jmax();++j)
            meanJ0 += problem.J0[nugrid.nu(j)][i*zgrid.dz]*(nugrid.nu(j)-nugrid.nu(j-1));
        cout <<std::fixed << std::setprecision(4)
             <<  zgrid.backtoz(i*zgrid.dz)
             <<"\t"<<problem.T[i]*physics::T0-273
             <<"\t\t"<<meanJ0<<endl;
    }
    return 0;
}

