The Milne problem of interest to us is

$$
	\mu\partial_x I+ I = \frac12\int_{-1}^1 I d\mu~~\tau\in(0,1),\qquad 
	I(0,\mu)|_{\mu>0}=0,~I(1,\mu)|_{\mu<0} =c_s\delta(\mu-\mu_s).
$$
It is an equation for the radiance $I$ versus altitude $x$ and ray direction cosine $\mu$. The radiation source is a collimated sunlight in direction $\mu_s$ of intensity $c_s$. Denote the directional averaged radiance
$$
J_0(x) = \frac12\int_{-1}^1 I d\mu.
$$
By the method of characteristics one shows that it is the solution of the integral equation
$$
J_0(x) = \frac12\int_0^1 E(|x-x'|)J_0(x')d x' + \frac{c_s}2 e^{-\frac{x}{\mu_s}}.
$$
~~~freefem
load "gsl" // gnu scientific library for E_1
meshL Th = segment(80,[x]); // 20 segments for the altitude grid x in (0,1)
fespace Vh(Th,P1);      // continuous pieccewise linear functions on (0,1)
Vh J0=0;                // Declares J0  in the FEM space
real cs=1, mus=1/sqrt(3.);
func real E1B(real xx){ return int1d(Th)(gslsfexpintE1(abs(x-xx))*J0(x));}
                        // computes int_0^1 E1(|x-xx|) S(x) dx. at xx
cout << "Convergence of J0(0)" <<endl;
for(int i=0; i<15;i++){ // ISIF
    J0= E1B(x)/2 + cs/2 * exp(-x/mus); // computes expression at all x
    cout<<J0(0)<<endl;  // to check convergence
}
cout<<" tau  J0(tau)"<<endl;
for(real tau=0;tau<1; tau+=0.03)
	cout<<tau<<" "<<J0(tau) <<endl;

~~~~

To compute $I(x,\mu)$ at some given $x,\mu$ one must approximate the Dirac mass at $\mu_s$ by an exponential and use the following when $\mu>0$

$$
I(x,\mu)\approx \int_0^x \frac1\mu e^{-\frac {|x-x'|}\mu}J_0(x')d x' +  c_s  e^{-\frac x{\mu_s}} e^{-\lambda(\mu-\mu_s)}
/\int_0^1 e^{-\lambda(\mu-\mu_s)}d\mu.
$$

~~~freefem
// compute I from J0 when mu>0
real mu;

func real convolp(real xx){ return  int1d(Th)((x<xx)*exp((x-xx)/mu)*J0(x)/mu);}

Vh I ;
cout<<endl;
real C=0, dmu=0.01, lambda=500;
for(mu=0.; mu<1;mu+=dmu)
	C += exp(-lambda*(mu-mus)^2)*dmu;
C = cs/C;
~~~
 When $\mu<0$,

$$
I(x,\mu)= -\int_x^1 \frac1\mu e^{\frac {|x-x'|}\mu}J_0(x')d x'.
$$

~~~freefem
// compute I from J0 when mu<0
func real convoln(real xx){ return - int1d(Th)((x>xx)*exp((x-xx)/mu)*J0(x)/mu);}

bool withSurfaceI = false;
if(withSurfaceI){
	for(mu=0.03; mu<1;mu+=0.03){
	 	for(real xx=0.;xx<1;xx+=0.03)
		  	cout<<mu<<" "<<xx<<" "<<convolp(xx) + C*exp(-lambda*(mu-mus)^2)*exp(-xx/mus)<<endl;
		cout<<endl;
	}
  for(mu=-1; mu<-0.03;mu+=0.03){
	for(real xx=0.;xx<1;xx+=0.03)
		cout<<mu<<" "<<xx<<" "<<convoln(xx)<<endl;
	cout<<endl;
  }
}
~~~
The last chunk of printed data can be displayed as a surface $x,\mu \to I(x,\mu)$ by gnuplot. Copy paste the data into a text file called "surfdata.txt"; launch gnuplot from the terminal window and type

splot "surfdata.txt" w l

If gnuplot is not installed you can uses "brew install gnuplot" on the mac.
As most of the radiance is into  the input 

C*exp(-lambda*(mu-mus)\^2)*exp(-xx/mus),

 you may set $C=0$ to see the rest of the signal.