# CHAPTER 3

The RTE in 3 dimensions is
The RTE in 3 dimensions is

$$
\frac1c\partial_t {I_\nu} + \omega\cdot\nabla {I_\nu}+\rho\kappa a_s{I_\nu} = {\frac{\rho\kappa a_s }{4\pi}\int_{S_2}}p_\nu(\omega,\omega'){I_\nu}(\omega'){d}\omega'
+ \rho\kappa(1-a_s) B_\nu(T),
$$

$c$ is the speed of light, $\omega$ is the ray direction, $I_\nu$ is the radiative intensity, $\rho$ is the air density, $\kappa$ is the absorption, $a_s$ is the scattering coefficent, $p_\nu$ is the probability that a ray in direction $\omega'$scatters in direction $\omega$. Many coefficeints depend on frequency $\nu$.  After scaling the Planck function is

$$
B_\nu(T) = \frac{\nu^3}{ e^{\frac\nu T} -1}
$$

It is coupled with the temperature equation

$$
\rho (\partial_t T+u\cdot\nabla T) -\nabla\cdot( \kappa_T\nabla T )
+ b\nabla\cdot\int_0^\infty\frac1{4\pi}\int_{S_2}
 {I_\nu}(\omega)\omega d\omega d\nu =0,
$$
where $b$ depends on the thermodynamic properties of air.  The thermal diffusion is $\kappa_T$ and the wind velocity $u$ is given by the Navier-Stokes equations with a Boussinesq source,
$$
\partial_t u+ u\cdot\nabla u-\frac{\mu_F}\rho\Delta u+ \nabla p={\bf g}(T),\quad \nabla\cdot u=0.
$$

In the book it is shown how to derive an integral equation for the angular mean of $I_\nu$, $J_0(x,y,z) = \frac12\int_{S_2}I_\nu d\omega$.

This integral equation is solved iteratively (ISIF), in the freefem script chamonix.xxx, xxx=edp or md.  The script reads a tetraedization of the Chamonix valley and use a, H-matrix library, htool and an interface to it for FreeFem written in C++ .

The absorption $\kappa$ as a function of $\nu$ is read from the file geminitransmittance.txt.
