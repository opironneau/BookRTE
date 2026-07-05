## Solution of the VRTE  in the grey case when the refactive index varies (without discontinuity) by solving the PDEs with SUPG regularization and iterations on the source.

For compatibility with freefem++ $\mu\in(-1,1)$ is $ x $ and $z\in(0,Z)$ is $y$. The $\tilde~$ stands for the integral over $(0,\infty)$ with respect to $\nu$ divited by $n(z)^2$.
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

where  $P_2(x)=\tfrac12(3x^2-1)$, $\kappa$ is absorption, $\kappa_s=a_s\kappa$ is scattering and $\kappa_a=(1-a_s)\kappa$.

Temperature is computed from the equation
$$
\sigma T^4 = J_0(y) := \tfrac12\int_{-1}^1 \tilde I(x,y)d x, ~~\sigma =\frac{\pi^4}{15}
$$
is the Stefan constant (divided by $n^2$). Note that $\kappa_a \sigma T^4+ \frac{\kappa_s}2\int_{-1}^1 \tilde I d x '= \kappa \sigma T^4$.  Hence when $\beta=0$, $a_s has no effect on the results.

We consider the following boundary conditions:

$$
I(x,0)|_{x>0} = x C_e \sigma T_e^4
\quad I(x,Z)|_{x<0} = 0
$$

~~~freefem
verbosity=0;
real Z=1.2,
    kappanu = 0.5, 
    as=0.3,
    kappas=kappanu * as;
real stefan=pi^4/15;
real Te = (18+273.)/4798;
real Ce=2.5;
real I0=stefan*Ce*Te^4;
real beta=0.5;

border a1(t=0,1) {x=t; y=0;}
border b(t=0,Z)  {x=1; y=t;}
border c1(t=1,0) {x=t;y=Z;}
border c2(t=0,-1){x=t;y=Z;}
border d(t=Z,0)  {x=-1; y=t;}
border a2(t=-1,0){x=t;y=0;}
int n=2;
mesh Th=buildmesh(a1(10*n)+b(10*n)+c1(10*n)+c2(10*n)+d(10*n)+a2(10*n));
fespace Vh(Th,P1);
fespace V2h(Th,[P1,P1]);
Vh T=Te/2;
V2h [I,Q], [Ih,Qh];

real a=0.2, z1=0.5, z2=0.8; // Diffraction index 
Vh nn=1+a*4*max(y-z1,0.)*max(z2-y,0.)/sqr(z1+z2), dlogn= a*2*(y>z1)*(y<z2)*(z1+z2-2*y)/sqr(z1+z2);

real eps=0.0;
macro kappay()(1-0.7*y) // altitude effect on kappa		
macro PP(x) (1.5*x*x-0.5) //
macro BB(T) (stefan*sqr(sqr(T))) //
func real J0(real yaux){return 0.5*int1d(Th,a1,a2)(I(x,yaux));}
func real QP(real yaux){return 0.5*int1d(Th,a1,a2)(PP(x)*I(x,yaux)-(1-PP(x))*Q(x,yaux));}


varf AA([I,Q],[Ih,Qh]) = int2d(Th)(x*dy(I)*Ih +x*dy(Q)*Qh+ kappay*kappanu*(I*Ih +Q*Qh)
                                +  dlogn * (1-x*x)*(dx(I)*Ih + dx(Q)*Qh)
    + eps*(x*x+0.001)*(dy(I)*dy(Ih)+dx(I)*dx(Ih) + dy(Q)*dy(Qh)+dx(Q)*dx(Qh))
    + eps*x*(kappay*kappanu*I+  dlogn * (1-x*x)*dx(I))*dy(Ih)
    + eps*x*(kappay*kappanu*Q + dlogn * (1-x*x)*dx(Q))*dy(Qh))
    +on(a1,I=x*I0,Q=0) + on(c2, I=0,Q=0) ;

 matrix A  =  AA(V2h,V2h) ;

 varf aa([I,Q],[Ih,Qh])  = 
	 	 int2d(Th)(kappay*(kappanu*BB(T)*(Ih+eps*x*dy(Ih)) 
                    + beta*kappas/2* ( PP(x)*QP(y)*(Ih+eps*x*dy(Ih)) 
                                       -(1-PP(x))*QP(y)*(Qh+eps*x*dy(Qh)))))
	+on(a1,I=x*I0,Q=0) + on(c2, I=0,Q=0);


for(int i=0; i<10;i++){
	real[int] bb = aa(0,V2h) ;
    I[]=A^-1*bb; 
	T=sqrt(sqrt( J0(y)/stefan ));
	cout<<T(0.,0.)*4798-273<<endl;
        if(i==7) {
            Th = adaptmesh(Th,I,verbosity=1, abserror=1, nbjacoby=2,
                err=0.003, nbvx=5000, omega=1.8, ratio=1.8, nbsmooth=3,
                splitpbedge=1, maxsubdiv=5, rescaling=1);
 //           plot(Th);
            A  =  AA(V2h,V2h) ;
            [I,Q]=[I,Q];
        eps=0.01;
        }
}
[I,Q]=[I*nn,Q*nn];
plot(I, dim=3, fill=true, value=true);
plot(Q, dim=3, fill=true, value=true);
T=T*4798-273;
plot(T, dim=3, fill=true, value=true);

ofstream Iout("Iout2.txt");
for(int i=0;i<=30;i++){
    real z=Z*i/30.;
    Iout<<z<<" "<<T(0.,z)<<" "<<I(0.5,z)<<" "<<Q(0.5,z)<<endl;
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
