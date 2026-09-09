# CHAPTER 4. VRTE

When polarization is important the radiative transfer equations have 2 unknowns, the radiative intensity I and the polarization intensity Q. Other polarization parameters are zero when Q is zero on the input boundaries.  VRTE is

$$
\mu \partial_z I + \kappa I =\sigma_a B_\nu(T) + \frac{\sigma_s}2\int_{-1}^1 I d\mu'
+ \frac{\beta\sigma_s}4 P_2(\mu)\int_{-1}^1 [P_2 I-(1-P_2 )Q] d\mu',
$$

$$
\mu \partial_z Q + \kappa Q = -\frac{\beta\sigma_s}4 (1-P_2(\mu))\int_{-1}^1 [P_2 I-(1-P_2 )Q] d\mu',
$$

$$
  \int_{\R_+}\sigma_a\big(B_\nu(T)
 - \tfrac12\int_{-1}^1 I d\mu\big) d\nu=0.
$$

where  $P_2(\mu)=\tfrac12(3\mu^2-1)$. 

The boundary conditions for a collimated sunlight at cosine angle $\mu_s$ and intensity $c_s$ is

$$
 I(Z,\mu)=c_s B_\nu(T_s) \delta(\mu+\mu_s), \quad Q(Z,\mu)=0,~~\forall\mu<0.
$$

An infrared Lanbertian source from the ground of intensity $c'_e$ combine with albedo of intensity $q_0$ is modelled by

$$
I(0,\mu) = c'_e B_\nu(T_e)\mu + q_0\int_{-1}^0 \mu I(0,\mu) d\mu,\quad Q(0,\mu)=0, ~~\forall\mu>0,
$$

## ISIF
 
Denote

$$
	J_q(z) = \tfrac12\int_{-1}^1 \mu^qI d\mu
\quad
	K_q(z) = \tfrac12\int_{-1}^1 \mu^qQ d\mu .
$$

Then, 

$$
J_q(z) = R_q +\tfrac{1}2\int_0^Z \left( E_{q+1}(\kappa_\nu|z-z'|)S_0(z')+ E_{q+3}(\kappa_\nu|z-z'|)S_2(z')\right) dz',
$$

$$
K_q(z) =\tfrac{1}2\int_{-1}^1\mu^q Q(z,\mu) d\mu =\tfrac12\int_0^Z \left(E_{q+1}(\kappa_\nu|z-z'|)S'_0(z') + E_{q+3}(\kappa_\nu|z-z'|)S'_2(z')\right) dz',
$$

with $c_0$ a function of $(1-a_s)B_\nu(T)+a_s J_0$, and with 

$$
R_q(z,\nu) = \tfrac{c'_e}2  B_\nu(T_e) E_{q+3}(\kappa_\nu z)  + \tfrac{c_s}2  B_\nu(T_s)\mu_s^q\kappa_\nu a_s e^{-\kappa_\nu\frac{Z-z}{\mu_s}}+
\tfrac{c_0}2 E_{q+2}(\kappa_\nu z),
$$

$$
H(z,\nu) = \frac{9\beta\sigma_s}8(J_2-\tfrac13 J_0 - K_0 + K_2),
$$

and where

$$
S_0=\sigma_a B_\nu + \sigma_s J_0 - \frac13 H,
\qquad
S_2 = H,
\qquad
S'_0=-H,
\qquad
S'_2 = H.
$$

So at each iteration we only need to compute, for $q=0,2$,

$$
J_q(z) = R_q +\tfrac12 \int_0^Z \left(E_{q+1}(\kappa_\nu|z-z'|)(\sigma_a B_\nu + \sigma_s J_0 - \frac13 H) +  E_{q+3}(\kappa_\nu|z-z'|)H \right) dz',
$$

$$
 K_q(z) =  -\tfrac12 \int_0^Z \left(E_{q+1}(\kappa_\nu|z-z'|) - E_{q+3}(\kappa_\nu|z-z'|)\right)H dz',
 $$
 and update  $T$ by  solving the temperature equation.  It gives the following iterative scheme
 
- Initialize $J_q,K_q,q=0,2$.
-  Compute $c_0$  and $R_q,~ q=0,2$ .
-  Compute $T$.
-  Update $J_q,K_q,~q=0,2$ .
  
## List of files

In this repository all files  with IQ in their names solve the above problem. Some are in C++ , one is in Matlab another in Fortran and others are in Python (trying to improve the speed).  Finally one is by the Monte Carlo method for comparison (written by Anthropic/claude.AI).

The file "twoRTlacJK.edp" is a combination of VRTE with the Navier-Stokes equations to simulate the heating of water of a pool by the sun. It is FreeFem++ file.
 


