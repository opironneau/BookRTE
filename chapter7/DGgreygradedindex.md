## Solution of the VVRTE by solving the PDE with DG-FEM

$$
{\mathcal L} u := \mu\partial_z u + \partial_z\log n (1-\mu^2)\partial_\mu u + \kappa(z) u  = f, \hbox{ in } \Omega,
u|_\Sigma=u_\Sigma,
$$

applied to

$$
{\mathcal L}([\tilde I,\tilde Q]^T) = S \hbox{ in }\Omega, \quad \tilde I|_\Sigma=I_\Sigma,\quad \tilde Q|_\Sigma=0,
$$
 
with 
 
$$
 \Omega=(-1,1)\times(0,Z),
$$
 
$$
 \Sigma=\left([0,1]\times\{0\}\right)\cup\left([-1,0]\times(\{Z\}\right).
$$

$$
 I(0,\mu>0)=C_e\mu -q_0\int_{-1}^0 I(0,\mu d\mu),~Q(0,\mu>0)=0,
$$
 
$$
 I(Z,\mu<0)= C_s e^{-\lambda(\mu+\mu_s)^2}, ,~Q(Z,\mu<0)=0
$$

The problem is defined by the following data

~~~freefem

real Z=1.2;  // TAO 
real Te = (18+273.)/4798, Ts = 5798./4798; // Earth and sun temperatures
real Cs = 0.2e-5, Ce=2.; // radiation intensities of the Sun and Earth
real q0=-0.3, xs=0.5; // cosine incident angle of sunlight
real beta=0.5;  // Proportion of Rayleigh scattering
real z1=0.5, z2=0.8;  // Cloud is between z1 and z2

border a1(t=0,1) {x=t; y=0;} // The computational domain is a rectangle
border b(t=0,Z)  {x=1; y=t;} // right vertical side
border c1(t=1,0) {x=t;y=Z;} // first half of top horizontal boundary
border c2(t=0,-1){x=t;y=Z;} // second half of horizontal boundary
border d(t=Z,0)  {x=-1; y=t;} // left vertical side
border a2(t=-1,0){x=t;y=0;} // left part of horizontalground boundary
int n=3; // controls the fineness of the triangulation
mesh Th=buildmesh(a1(10*n)+b(10*n)+c1(10*n)+c2(10*n)+d(10*n)+a2(10*n));
fespace Vh(Th,P1dc);  // P1 discontinuous functions on Th

real kappanu=0.5;  // absorption is this times density
real a=0.4; // Diffraction index parameter: nn is n^2

Vh  nn=1+a*4*max(y-z1,0.)*max(z2-y,0.)/sqr(z1+z2), 
	dlogn = a*2*(y>z1)*(y<z2)*(z1+z2-2*y)/sqr(z1+z2),
	kappa=kappanu*(1-0.7*y), 
	as=0.3 + 0.3*(4*max(y-z1,0.)*max(z2-y,0.)/sqr(z1+z2)) ,
	sigmas=kappa * as, bk2= beta*sigmas/2, // as is scattering coeff.
	T=Te, I=0, Q=0, Ih, Qh; // Temperature, radiation and polarization intensities
	
~~~ 

Notice that the scattering coefficient $a_s$ and the refraction index nn are made bigger in the cloud $(z_1,z_2)$.
The boundary functions and the right hand sides of the PDEs are defined below.

~~~freefem  

real  stefan = pi^4/15, I0=stefan*Ce*Te^4, Is=stefan*Cs*Ts^4;

macro PP(x) (1.5*x*x-0.5) // EOM
macro kaint(xp,xpp) (kappanu*(xpp-xp)*(1-0.35*(xp+xpp))) // int_xp^xpp kappa
macro BB(T) (stefan*sqr(sqr(T))) // nu integral of the Planck function

~~~

The functions $J_0$ and $Q_p$ return the $\mu$-integrals of $I$ and of
the terms that appear in the right hand side of the PDE.
~~~freefem
func real J0(real yaux){return 0.5*int1d(Th,a1,a2)(I(x,yaux));}
func real QP(real yaux){return 0.5*int1d(Th,a1,a2)(PP(x)*I(x,yaux)-(1-PP(x))*Q(x,yaux));}

~~~

The 2 PDEs are now defined and solved with the addition of the integrals of the jump accross edges as needed by DG. Theses are embedded in the source iteration loop

~~~freefem

macro bx(x) (dlogn*(1-x*x)) //
macro by(x) (x) //

for(int i=0; i<12;i++){
	
	real JmT=Is*exp(-kaint(0,Z)/xs) + int1d(Th,a2)(I(x,0)*x);
	
	real al=1; // recommanded penalty coefficient
	solve BBI(I,Ih) = int2d(Th)(( kappa*I +  bx(x)*dx(I) + by(x)*dy(I) )*Ih)
	//				+ intalledges(Th)((1-nTonEdge)*Ih*(al*abs(N.x*dlogn * (1-x*x) + N.y*x)-(N.x*dlogn * (1-x*x) + N.y*x)/2)*jump(I))
					- intalledges(Th)(al*(1-nTonEdge)*Ih*(min(0.,N.x*bx(x) + N.y*by(x)))*jump(I))
					- int2d(Th)(  (kappa*BB(T)  + bk2*PP(x)*(QP(y) +PP(xs)*Is*exp(-kaint(y,Z)/xs))) *Ih)
	 				- int1d(Th,a1,c2)(min(0.,N.x*bx(x)+N.y*by(x))*I*Ih) 
					+ int1d(Th,a1)(min(0.,N.x*bx(x)+N.y*by(x))*(x*I0 + q0*JmT)*Ih) 
;

	solve BBQ(Q,Qh) = int2d(Th)((x*dy(Q) + kappa*Q +  dlogn * (1-x*x)*dx(Q) )*Qh)
					- intalledges(Th)(al*(1-nTonEdge)*Qh*(min(0.,N.x*bx(x) + N.y*by(x)))*jump(Q))
					- int1d(Th,a1,c2)(min(0.,N.x*bx(x)+N.y*by(x))*Q*Qh) 
					+ int2d(Th) (bk2*(1-PP(x))*QP(y)*Qh )
;

~~~
The temperature equation has an explicit solution 
~~~freefem

	T=sqrt(sqrt((abs(J0(y))+0.5*Is*exp(-kaint(y,Z)/xs))/stefan ));
	cout<<T(0.,0.)*4798-273<<endl;
~~~
At iteration 7 the mesh is adapted to the solution
~~~freefem	

    if(i==7) {
        Th = adaptmesh(Th,I,verbosity=1, abserror=1, nbjacoby=2,
            err=0.0015, nbvx=5000, omega=1.8, ratio=1.8, nbsmooth=3,
            splitpbedge=1, maxsubdiv=5, rescaling=1);
        plot(Th);
		I=I; Q=Q;
    }
}
~~~
Results are ploted and written on a file
~~~freefem

I = I*nn; Q=Q*nn; T=T*4798-273;
plot(I, dim=3, fill=true, value=true, cmm="I");
plot(Q, dim=3, fill=true, value=true, cmm="Q");
plot(T, dim=3, fill=true, value=true, cmm="T");
cout<<endl;
ofstream Iout("Iout4.txt");
for(int i=0;i<60;i++){
    real z=Z*i/60.;
    Iout<<z<<" "<<T(0.1,z)<<" "<<I(0.1,z)<<" "<<Q(0.1,z)<<endl;
    cout << T(0.1,z)<<endl;
}

~~~
Use P and <CR> to navigate the figue. End it with <escape>