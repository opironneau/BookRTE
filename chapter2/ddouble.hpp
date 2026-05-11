// file ddouble.h, for automatic differentiation (single variable)
// adapted from M. Grundmann's MIPAD

#ifndef _DDOUBLE__H_
#define _DDOUBLE__H_

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <time.h>
using namespace std;
#
static const double sqrtpi = sqrt(4*atan(1.0));

class ddouble 
{
 public:  double val[2]; // val[0] is value and val[1] is the derivative
 ddouble() { val[0] = 0; val[1] = 0;}
 ddouble(const ddouble& a){ val[0] = a.val[0]; val[1] = a.val[1]; }
 ddouble(double a, double b=0){ val[0]=a; val[1]=b;} // specifies var for diff
 
 ddouble& operator=(double a)
 	{ val[0] = a; val[1] = 0.0; return *this;} 		
 ddouble& operator=(const ddouble& a)
 	{ val[0] = a.val[0]; val[1] = a.val[1]; return *this;}

		
 double& operator[](const int ii){ return this->val[ii];}
 double operator[](const int ii) const { return this->val[ii];}

  ddouble& operator  + (){return *this;};
  ddouble& operator += (double);
  ddouble& operator += (const ddouble&);
  ddouble& operator -= (double );
  ddouble& operator -= (const ddouble&);
  ddouble& operator *= (double);
  ddouble& operator *= (const ddouble&);
  ddouble& operator /= (double) ;
  ddouble& operator /= (const ddouble&) ;

  ddouble operator++(int);
  ddouble operator--(int);
  ddouble&  operator++();
  ddouble&  operator--();
  
  friend ostream& operator << (ostream&, const ddouble&);
  friend ddouble& operator << (ddouble&,double);
  
  friend ddouble parameter(double);
  
  friend int operator != (const ddouble&,const ddouble&);
  friend int operator != (double,const ddouble&);
  friend int operator != (const ddouble&,double);
  friend int operator == (const ddouble&,const ddouble&);
  friend int operator == (double,const ddouble&);
  friend int operator == (const ddouble&,double);
  friend int operator >= (const ddouble&,const ddouble&);
  friend int operator >= (double,const ddouble&);
  friend int operator >= (const ddouble&,double);
  friend int operator <= (const ddouble&,const ddouble&);
  friend int operator <= (double,const ddouble&);
  friend int operator <= (const ddouble&,double);
  friend int operator > (const ddouble&,const ddouble&);
  friend int operator > (double,const ddouble&);
  friend int operator > (const ddouble&,double);
  friend int operator < (const ddouble&,const ddouble&);
  friend int operator < (double,const ddouble&);
  friend int operator < (const ddouble&,double);


  friend ddouble operator + (const ddouble& x); //{return x + 0.0 ;} ; 
  friend ddouble operator + (const ddouble&,const ddouble&); 
  friend ddouble operator + (double, const ddouble&); 
  friend ddouble operator + (const ddouble&, double); 
  friend ddouble operator - (const ddouble& x ,double y); //{return (-y)+x;}; 
  friend ddouble operator - (const ddouble&,const ddouble&); 
  friend ddouble operator - (double, const ddouble&); 
  friend ddouble operator - ( const ddouble& );
  friend ddouble operator * (const ddouble&,const ddouble&); 
  friend ddouble operator * (double, const ddouble& ); 
  friend ddouble operator * (const ddouble& x, double y);  
  friend ddouble operator / (const ddouble& x, double y); 
  friend ddouble operator / (const ddouble&,const ddouble&); 
  friend ddouble operator / (double,const ddouble&) ; 
  friend ddouble exp (const ddouble&); 
  friend ddouble log (const ddouble&) ; 
  friend ddouble sqrt (const ddouble&) ;
  friend ddouble sin (const ddouble&); 
  friend ddouble cos (const ddouble&);
  friend ddouble tan (const ddouble&);
  friend ddouble pow (const ddouble&,double);
//  friend ddouble pow (const ddouble&,const ddouble&);
//  friend ddouble pow (const ddouble&, const int);
  friend ddouble abs (const ddouble&) ;
  friend ddouble erfc (const ddouble&) ;
 
};

inline double sign(const ddouble& x)
{ return ( x < 0.0 ? -1.0 : 1.0); }

inline double sign(const ddouble& x, double y)
{ return ( x < 0.0 ? -fabs(y) : fabs(y)); }

// f2c

inline ddouble d_abs(ddouble * x){ return abs(*x);	}
inline ddouble d_cos(ddouble * x){ return cos(*x);	}
inline ddouble d_sin(ddouble * x){ return sin(*x);	}
inline ddouble d_tan(ddouble * x){ return tan(*x);	}
inline ddouble d_exp(ddouble * x){ return exp(*x);	}
inline ddouble d_log(ddouble * x){ return log(*x);	}

inline ddouble d_sign(ddouble * x){ return sign(*x);	}
inline ddouble d_sign(ddouble * x,double*y){ return sign(*x,*y);	}

inline ddouble d_sqrt(ddouble * x){ return sqrt(*x);	}

// inline ddouble pow_dd(ddouble * x,ddouble*y)	{ return pow(*x,*y);	}
// inline ddouble pow_dd(double * x,ddouble*y) 	{ return pow(*x,*y);	}
// inline ddouble pow_dd(ddouble * x,double*y) 	{ return pow(*x,*y);	}
// inline ddouble pow_di(ddouble * x,int*y)	{ return pow(*x,*y);	}
		
/////////////////// Implementation Part /////////////////////////
// file ddouble.cpp, for automatic differentiation (single variable)

const double eps = 1.0e-50; // to avoid NaN when y=0 in sqrt(y)

ostream& operator<<(ostream& f, const ddouble& a)
{ f << "[" << a[0] << ','<< a[1] << "]"; return f;}	

ddouble ddouble::operator++(int)  
{	ddouble r=(*this);  r[0]++; return r;}

ddouble ddouble::operator--(int)  
{	ddouble r=(*this);  r[0]--; return r;}

ddouble& ddouble::operator++() 
{  (*this)[0]++; return *this;}

ddouble& ddouble::operator--() 
{  (*this)[0]--; return *this;}


ddouble& ddouble::operator += (double y) 
{  (*this)[0] += y; return *this; } 


ddouble operator - (const ddouble& a)
{ ddouble r; r[0] = -a[0]; r[1] = -a[1];  return r;} 	

ddouble& ddouble::operator -= (double y) 
{  (*this)[0]-=y; return *this;}

ddouble& ddouble::operator += (const ddouble& y) 
{ 	(*this)[0]+=y[0];(*this)[1]+=y[1];  return *this; }

ddouble& ddouble::operator -= (const ddouble& y) 
{ 	(*this)[0]-=y[0];(*this)[1]-=y[1];  return *this; }


ddouble& ddouble::operator *= (double y) 
{ (*this)[0] *=y; (*this)[1] *=y;  return *this;}

ddouble& ddouble::operator *= (const ddouble& y) 
{ return *this = *this * y;}

ddouble& ddouble::operator /= (const ddouble& y)  
{  return *this = *this / y;}

ddouble& ddouble::operator /= (double y) 
{ const double inv = 1.0 / y; 	
	(*this)[1]	*= inv; (*this)[1]	*= inv;
return *this;
}

int operator != (const ddouble& u,const ddouble& v)
{  return u[0] != v[0];}

int operator != (double u,const ddouble& v)
{  return u != v[0];}

int operator != (const ddouble& v,double u)
{  return v[0] != u;}

int operator == (const ddouble& u,const ddouble& v)
{  return u[0] == v[0];}

int operator == (double u,const ddouble& v)
{  return u == v[0];}

int operator == (const ddouble& v,double u)
{  return v[0] == u;}


int operator <= (const ddouble& u,const ddouble& v)
{  return u[0] <= v[0];}

int operator <= (double u,const ddouble& v)
{  return u <= v[0];}

int operator <= (const ddouble& v,double u)
{  return v[0] <= u;}

int operator >= (const ddouble& u,const ddouble& v)
{  return u[0] >= v[0];}

int operator >= (double u,const ddouble& v)
{  return u >= v[0];}

int operator >= (const ddouble& v,double u)
{  return v[0] >= u;}


int operator > (const ddouble& u,const ddouble& v)
{  return u[0] > v[0];}

int operator > (double u,const ddouble& v)
{  return u > v[0];}

int operator > (const ddouble& v,double u)
{  return v[0] > u;}

int operator < (const ddouble& u,const ddouble& v)
{  return u[0] < v[0];}

int operator < (double u,const ddouble& v)
{  return u < v[0];}

int operator < (const ddouble& v,double u)
{  return v[0] < u;}

ddouble operator + (const ddouble& x, const ddouble& y) 
{
  ddouble r;
  	r[0]	= x[0] + y[0];r[1]	= x[1] + y[1];
  return r;
}

ddouble operator + (double x, const ddouble& y) 
{   ddouble	r(y);  r[0]	+= x;  return r;}

ddouble operator + (const ddouble& y, double x) 
{   ddouble	r(y);  r[0]	+= x;  return r;}

ddouble operator - (const ddouble& x, const ddouble& y)
{   ddouble r;  r[0] = x[0] - y[0];r[1]	= x[1] - y[1]; return r;}

ddouble operator - (double x, const ddouble& y)
{  ddouble r; r[1]	= - y[1];  r[0] = x - y[0];	  return r;}

ddouble operator - (const ddouble& x, double y)
{  ddouble	r(x);  r[0]	-= y;  return r;		}

ddouble operator * (const ddouble& x, const ddouble& y)
{ ddouble r; r[0] = x[0]*y[0]; r[1]=x[0]*y[1]+x[1]*y[0];return r;}


ddouble operator * (double x, const ddouble& y)
{ddouble r;	r[0]	= x * y[0];	r[1]	= x * y[1];	return r;}

ddouble operator * ( const ddouble& y, double x)
{return x * y;}


ddouble operator / (const ddouble& x, const ddouble& y) 
{  ddouble r; r[0] = x[0]/y[0]; r[1]=(x[1]-x[0]*y[1]/y[0])/y[0];return r;}


ddouble operator / (double x, const ddouble& y) 
{ ddouble r; r[0] = x/y[0]; r[1]=-x*y[1]/y[0]/y[0];return r;}

ddouble operator/(const ddouble& x, double y)
{ddouble r; r[0] = x[0]/y; r[1]=x[1]/y;return r;}
	

ddouble exp (const ddouble& x)
{  ddouble r;r[0] = exp(x[0]); r[1] = x[1]*r[0];return r;}

ddouble log (const ddouble& x) 
{ ddouble r;r[0] = log(x[0]);r[1]	  = x[1]/x[0];return r;}

ddouble sqrt (const ddouble& x) 
{ ddouble r;r[0] = sqrt(x[0]); r[1] = 0.5*x[1]/(eps+r[0]);return r;}


ddouble sin (const ddouble& x) 
{   ddouble r; r[0]=sin(x[0]); r[1]=x[1]*cos(x[0]);  return r;}


ddouble cos (const ddouble& x) 
{   ddouble r; r[0]=cos(x[0]); r[1]=-x[1]*sin(x[0]);  return r;}


ddouble tan (const ddouble& x) { return (sin(x) / cos(x));}

ddouble pow (const ddouble& x,double y) {return exp(log(x) * y);}
// ddouble pow (const ddouble& x,const int y) {return exp(log(x) * (double)y);}
// ddouble pow (const ddouble& x,const ddouble& y)  {return exp(log(x) * y);}

ddouble abs (const ddouble& x) 
{ ddouble y;if(x[0] >= 0) y=x; else y = -x; return y;}

ddouble erfc (const ddouble& x) 
{ ddouble y;y[1] = 4*x[0]*x[1]*exp(-x[0]*x[0])/sqrtpi; y[0] = erfc(x[0]); return y;}

#endif

