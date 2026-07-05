//
//  IQ4class.hpp
//  IQ4++
//
//  Created by Olivier Pironneau on 02/05/2026.
//

#ifndef IQ4class_hpp
#define IQ4class_hpp

#include <stdio.h>

#ifndef IQ4_HAVE_GSL
#  if __has_include(<gsl/gsl_sf_expint.h>)
#    define IQ4_HAVE_GSL 1
#  else
#    define IQ4_HAVE_GSL 0
#  endif
#endif

#if IQ4_HAVE_GSL
#include <gsl/gsl_sf_expint.h>
#endif

namespace ExpIntegral {
#if IQ4_HAVE_GSL
    inline double E1(double x) { return gsl_sf_expint_E1(x); }
    inline double En(int n, double x) { return gsl_sf_expint_En(n, x); }
    inline double E2(double x) { return En(2, x); }
    inline double E3(double x) { return En(3, x); }
    inline double E4(double x) { return En(4, x); }
    inline double E5(double x) { return En(5, x); }
#else
    inline double E1(const double t=1){
         double t1=fabs(t);
        const int Kexpint = 9+(t1-1)*4;
        const double  gaNtaua =0.577215664901533;
        if(t1<1e-5) return 0;
        if(t1>4) {
            double tx=1./t1;
            return exp(-t1)*tx * (1 +(-1+(2+(-6+(24+(-120+720*tx)*tx)*tx)*tx)*tx ));
        }
        double ak=t1, soNtaue=-gaNtaua - log(t1)+ak;
        for(int k=2;k<Kexpint;k++){
            ak *= -t1*(k-1)/k/k;
            soNtaue += ak;
        }
        return soNtaue;
    }
    inline double E2(const double t=1){
        double t1=fabs(t);
        return exp(-t1) - t1*E1(t1);
    }
    inline double E3(const double t=1){
        double t1=fabs(t);
        return (exp(-t1) - t1*E2(t1))/2;
    }
    inline double E4(const double t=1){
        double t1=fabs(t);
        return (exp(-t1) - t*E3(t1))/3;
    }
    inline double E5(const double t=1){
        double t1=fabs(t);
        return (exp(-t1) - t*E4(t1))/4;
    }
#endif
}

#endif /* IQ4class_hpp */
