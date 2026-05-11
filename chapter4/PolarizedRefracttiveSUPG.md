Solution of the VRTE  in the grey case by SUPG and iterations on the source.
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
\frac\sigma{n^2} T^4 = J_0(y) := \tfrac12\int_{-1}^1 \tilde I(x,y)d x, ~~\sigma =\frac{\pi^4}{15}
$$
is the Stefan constant. Note that $\kappa_a \frac\sigma{n^2} T^4+ \frac{\kappa_s}2\int_{-1}^1 \tilde I d x '= \kappa \frac\sigma{n^2} T^4$.  Hence when \beta=0$, $a_s has no effect on the results.

We consider the following boundary conditions:

$$
I(x,0)|_{x>0} = x C_e \sigma T_e^4
\quad I(x,Z)|_{x<0} = 0
$$

~~~freefem
verbosity=0;
real kappa = 0.5, 
    as=0.5,
    kappas=kappa * as, 
    kappaa=kappa*(1-as);
real stefan=pi^4/15;
real Te = (18+273.)/4798;
real Ce=2.4;
real I0=stefan*Ce*Te^4;
real beta=0.5;

border a1(t=0,1) {x=t; y=0;}
border b(t=0,1)  {x=1; y=t;}
border c1(t=1,0) {x=t;y=1;}
border c2(t=0,-1){x=t;y=1;}
border d(t=1,0)  {x=-1; y=t;}
border a2(t=-1,0){x=t;y=0;}
int n=2;
mesh Th=buildmesh(a1(10*n)+b(10*n)+c1(10*n)+c2(10*n)+d(10*n)+a2(10*n));
fespace Vh(Th,P1);
fespace V2h(Th,[P1,P1]);
Vh T=Te*(1-0.7*y);
V2h [I,Q], [Ih,Qh];
real a=0.; // Diffraction index is 1+a*z
//Vh nn=sqrt(1+a*sqr(y-0.5)), dlognn=a*(y-0.5)/(1+a*sqr(y-0.5));
Vh nn=1+a*y, dlognn=a;
real eps=0.02;
//Vh phi=sqr(1+0.01*(y>0.5)*(y<0.7))*(1-x*x); plot(phi,wait=1);
//exit(0);

macro kappay()(1-0.7*y) // altitude effect on kappa		
macro PP(x) (1.5*x*x-0.5) //
macro BB(T) (stefan*sqr(sqr(T))) //
func real J0(real yaux){return 0.5*int1d(Th,a1,a2)(I(x,yaux));}
func real QP(real yaux){return 0.5*int1d(Th,a1,a2)(PP(x)*I(x,yaux)-(1-PP(x))*Q(x,yaux));}


varf AA([I,Q],[Ih,Qh]) = int2d(Th)(x*dy(I)*Ih +x*dy(Q)*Qh+ kappay*kappa*(I*Ih +Q*Qh)
                                +  dlognn * (1-x*x)*(dx(I)*Ih + dx(Q)*Qh)
    + eps*(x*x+0.001)*(dy(I)*dy(Ih)+dx(I)*dx(Ih) + dy(Q)*dy(Qh)+dx(Q)*dx(Qh))
    + eps*x*(kappa*I+  dlognn * (1-x*x)*dx(I))*dy(Ih)
    + eps*x*(kappa*Q + dlognn * (1-x*x)*dx(Q))*dy(Qh))
    +on(a1,I=x*I0/sqr(nn),Q=0) + on(c2, I=-x*I0/sqr(nn),Q=0) ;

 matrix A  =  AA(V2h,V2h) ;

 varf aa([I,Q],[Ih,Qh])  = 
	 	 int2d(Th)(kappay*(kappa*BB(T)/sqr(nn)*(Ih+eps*x*dy(Ih)) 
                    + beta*kappas/2* ( PP(x)*QP(y)*(Ih+eps*x*dy(Ih)) 
                                       -(1-PP(x))*QP(y)*(Qh+eps*x*dy(Qh)))))
	+on(a1,I=x*I0/sqr(nn),Q=0) + on(c2, I=0,Q=0);


for(int i=0; i<10;i++){
	real[int] bb = aa(0,V2h) ;
    I[]=A^-1*bb; 
	T=sqrt(nn*sqrt( J0(y)/stefan ));
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
plot(I, dim=3, fill=true, value=true);
plot(Q, dim=3, fill=true, value=true);
T=T*4798-273;
plot(T, dim=3, fill=true, value=true);

ofstream Iout("Iout2.txt");
for(int i=0;i<=30;i++){
    real z=i/30.;
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
