# Navier-Stokes Equations

In chapter 4 we will compute the temperature in a water pond subject to wind and sunlight. We begiin witht he computation of the velocity field by solving
$$
\frac{\partial u}{\partial t}
\partial_t{\bf u}+{\bf u}\cdot\nabla{\bf u}
-\nu\Delta {\bf u}  -\nabla p =0,~~~\nabla\cdot{\bf u}=0~~~in ~ \Omega, ~~t>0.
$$
The initial velocity is given and the velocity on the surface is equal to the wind velocity. A no-slip condition is imposed on the boundary in contact with the ground: $u=v=0$.

The domain is defined by 2 curves $t\to x(t),y(t)$
~~~freefem
int n=4;
border a1(t=-1,1){x=t; y= -sqrt(sqrt((1-x*x)));}
border a2(t=1,-1){x=t; y=0;}
mesh Th=buildmesh(a1(40*n)+a2(10*n));
Th = movemesh(Th,[3*x,(1-0.2*x)*y]);
plot(Th,wait=true,cmm="Hit the Return key");
~~~~
The pressure is approximated by a FEM $P^1$ function on the triangulation $T_h$. The velocity components are approximated by a FEM $P^1$ function on the triangulation $T_h$.
~~~freefem
fespace Vh(Th,P1);
fespace Uh(Th,P2);
Vh p,q;
Uh u=1,v=0,uh,vh,ux,uy;

real visc=0.05 ,eps=1.e-5, dt=0.02;// viscosity, pregularization and time step
~~~~
A Characteristic Galerkin discretization is applied:
$$
\partial_t{\bf u}+{\bf u}\cdot\nabla{\bf u}|_{x,t}\approx \frac1{\delta t}[{\bf u}(x,t)-{\bf u}(x-{\bf u}(x,t-\delta t)\delta t,t-\delta t)]
$$
An implicit scheme used, leading to
~~~freefem
for(int tstep=1;tstep<70;tstep++)
{
	ux=u; uy=v;
	solve NS(u,v,p,uh,vh,q) 
	   = int2d(Th)((u*uh+v*vh)/dt 
	   + (dx(q)*u+dy(q)*v)-(dx(p)*uh+dy(p)*vh)  + eps*p*q
	   + visc*(dx(u)*dx(uh)+dy(u)*dy(uh)
	   		+dx(v)*dx(vh)+dy(v)*dy(vh)))
	   -int2d(Th)((convect([ux,uy],-dt,ux)*uh
	            +  convect([ux,uy],-dt,uy)*vh)/dt)
	+on(1,u=0,v=0)+on(2,u=10,v=0);// boundary conditions

        if(tstep==20){ // build a better mesh at time step 20
	        Th = adaptmesh(Th, [u, v], q, err=0.005); 
	       plot(Th,cmm="Hit the Return key",wait=true);
	    }
     plot([u,v],cmm="Time step "+tstep, wait=0);
}
plot([u,v],wait=0,cmm="Vector [u,v] is shown. Hit the Escape key");// , ps="lakevelocity.ps"); 
// remove comment to obtain a .ps copy
~~~






