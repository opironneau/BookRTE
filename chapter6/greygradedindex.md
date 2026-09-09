## Solution of the VRTE  in the grey case when the refactive index varies (without discontinuity) by solving the PDEs with SUPG regularization and iterations on the source.

For compatibility with freefem++ $\mu\in(-1,1)$ is $ x $ and $z\in(0,Z)$ is $y$. The $\tilde~$ stands for the integral over $(0,\infty)$ with respect to $\nu$ divided by $n(z)^2$.
$$
\left\{
\begin{aligned}&
 x  \partial_y\tilde I + \partial_y\log n (1- x ^2)\partial_ x \tilde I+ \kappa\tilde I 
\cr&
\hskip 1cm=\kappa_a \frac\sigma{n^2} T^4+ \frac{\kappa_s}2\int_{-1}^1 \tilde I d x '
+\frac{\beta\kappa_s}4 P_2( x )\int_{-1}^1 [P_2\tilde I-(1-P_2 )\tilde Q] d x ',
\cr&
 x  \partial_y \tilde Q + \partial_y\log n (1- x ^2)\partial_ x \tilde Q + \kappa\tilde Q 
\cr&
\hskip 1cm= -\frac{\beta\kappa_s}4 (1-P_2( x ))\int_{-1}^1 [P_2\tilde I-(1-P_2 )\tilde Q] d x ',
\end{aligned}
\right.
$$

where  $P_2(x)=\tfrac12(3x^2-1)$, $\kappa=\kappa_a+\kappa_s$ is the extinction, $\kappa_s=a_s\kappa$ is scattering and $\kappa_a=(1-a_s)\kappa$.

Temperature is computed from the equation
$$
\sigma T^4 = n^2 J_0(y), \quad J_0(y) := \tfrac12\int_{-1}^1 \tilde I(x,y)d x, ~~\sigma =\frac{\pi^4}{15}
$$
is the Stefan constant. Note that $\kappa_a \frac{\sigma T^4}{n^2}+ \frac{\kappa_s}2\int_{-1}^1 \tilde I d x '= \kappa \frac{\sigma T^4}{n^2}$.  Hence when $\beta=0$, $a_s$ has no effect on the results.

We consider the following boundary conditions:

$$
I(x,0)|_{x>0} = x C_e \sigma T_e^4
\quad I(x,Z)|_{x<0} = 0
$$

~~~freefem
// I-Q VRRTE Handles continuous refractive index and constant kappa
// No albedo
// See same file with .md suffix for more explanations
// See in Chapter 7 for the same using DG and .md for more expanations
real Z=1.2;
real stefan=pi^4/15;
real Te = (18+273.)/4798, Ts = 5800./4798;
real Cs = 0.1e-5, Ce=2.;
real q0=-0.3, xs=0.5; // cosine incident angle of sunlight
real beta=0.5;
real z1=0.5, z2=0.8;

border a1(t=0,1) {x=t; y=0;}
border b(t=0,Z)  {x=1; y=t;}
border c1(t=1,0) {x=t;y=Z;}
border c2(t=0,-1){x=t;y=Z;}
border d(t=Z,0)  {x=-1; y=t;}
border a2(t=-1,0){x=t;y=0;}
int n=2;
mesh Th=buildmesh(a1(10*n)+b(10*n)+c1(10*n)+c2(10*n)+d(10*n)+a2(10*n));
fespace Vh(Th,P1);

real eps=0, kappanu=0.5;
real a=0.4;               // refractive index param (nn is n^2); a_n of (6.15)
string outfile="Iout4.txt"; // rename to match when a is changed (0 -> Iout0, 0.2 -> Iout2)

Vh  nn=1+a*4*max(y-z1,0.)*max(z2-y,0.)/sqr(z2-z1), 
	dlogn = a*2*(y>z1)*(y<z2)*(z1+z2-2*y)/sqr(z2-z1)/nn, // = (1/2)(n^2)'/n^2
	kappa=kappanu*(1-0.7*y), 
	as=0.3 + 0.3*(4*max(y-z1,0.)*max(z2-y,0.)/sqr(z2-z1)) ,
	sigmas=kappa * as, bk2= beta*sigmas/2, 
	T=Te/2, I=0, Q=0, Ih, Qh;
	  
real  I0=stefan*Ce*Te^4, Is=stefan*Cs*Ts^4;

macro PP(x) (1.5*x*x-0.5) //
macro kaint(xp,xpp) (kappanu*(xpp-xp)*(1-0.35*(xp+xpp))) // int_xp^xpp kappa
macro BB(T) (stefan*sqr(sqr(T))) //
func real J0(real yaux){return 0.5*int1d(Th,a1,a2)(I(x,yaux));}
func real QP(real yaux){return 0.5*int1d(Th,a1,a2)(PP(x)*I(x,yaux)-(1-PP(x))*Q(x,yaux));}

for(int i=0; i<10;i++){
	
	// int_{-1}^0 mu.[c_s.sigma.Ts^4.e.delta(mu+mu_s)]dmu = -mu_s.c_s.sigma.Ts^4.e
	real JmT=-xs*Is*exp(-kaint(0,Z)/xs) + int1d(Th,a2)(I(x,0)*x);
	
	solve BBI(I,Ih) = int2d(Th)( (x*dy(I) + kappa*I +  dlogn * (1-x*x)*dx(I))
						* (Ih + eps*(x*dy(Ih) + kappa*Ih +  dlogn * (1-x*x)*dx(Ih)))  )
	- int2d(Th)(  (kappa*BB(T)/nn  + bk2*PP(x)*(QP(y) +0.5*PP(xs)*Is*exp(-kaint(y,Z)/xs))) 
					*(Ih + eps*(x*dy(Ih) + kappa*Ih +  dlogn * (1-x*x)*dx(Ih))))
    + on(a1,I=x*I0 + q0*JmT) + on(c2, I=0) ;

	solve BBQ(Q,Qh) = int2d(Th)((x*dy(Q) + kappa*Q +  dlogn * (1-x*x)*dx(Q))
					* (Qh + eps*(x*dy(Qh) + kappa*Qh +  dlogn * (1-x*x)*dx(Qh)))  )
	+ int2d(Th) (bk2*(1-PP(x))*(QP(y) +0.5*PP(xs)*Is*exp(-kaint(y,Z)/xs))
				*(Qh+ eps*(x*dy(Qh) + kappa*Qh +  dlogn * (1-x*x)*dx(Qh)) ))
	+ on(a1,c2,Q=0);
	
	T=sqrt(sqrt(nn*(J0(y)+0.5*Is*exp(-kaint(y,Z)/xs))/stefan ));
	cout<<T(0.,0.)*4798-273<<endl;
	
    if(i==7) {
        Th = adaptmesh(Th,I,Q,verbosity=1, abserror=1, nbjacoby=2,
            err=0.0015, nbvx=5000, omega=1.8, ratio=1.8, nbsmooth=3,
            splitpbedge=1, maxsubdiv=5, rescaling=1);
        plot(Th);
		I=I; Q=Q;
		eps=0.05;
    }
}
I = I*nn; Q=Q*nn; T=T*4798-273;
plot(I, dim=3, fill=true, value=true);
plot(Q, dim=3, fill=true, value=true);
plot(T, dim=3, fill=true, value=true);
cout<<endl;
ofstream Iout(outfile);
for(int i=0;i<60;i++){
    real z=Z*i/60.;
    Iout<<z<<" "<<T(0.1,z)<<" "<<I(0.1,z)<<" "<<Q(0.1,z)<<endl;
    cout << T(0.1,z)<<endl;
}

~~~~
Case 1:  $\beta=0, \kappa_s=0.2, a=0$

|------------------------|
|![](as2.png)         |
|------------------------|

Case 2:  $\beta=0.5, \kappa_s=0.2$, a=0$

|------------------------|
|![](betaas.png)         |
|------------------------|

Case 3:  $\beta=0.2, \kappa_s=0.2, a=0.2$


|------------------------|
|![](betaasa.png)         |
|------------------------|

Case 4:  $\beta=0.2, \kappa_s=0.2, a=-0.2$


|------------------------|
|![](betaasam.png)         |
|------------------------|
