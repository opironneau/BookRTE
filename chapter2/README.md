# CHAPTER 2

The temperature $T$ in the atmosphere is the main objective of the chapter.  The atmosphere is wroughly a black body emitting EM waves according to Planck's law. A new equation is added to the RTE: the conservation of energy.  In the stratified case, the system is now for $[I_\nu(z,\mu),T(z)]$ in $(0,Z)\times(-1,1)$,

$$
\mu\partial_z I_\nu + \kappa I_\nu =\sigma_a B_\nu(T)+  {\sigma_s}\tfrac12\int_{-1}^1 I_\nu(z,\mu')d\mu ,
$$

$$
I_\nu(0,\mu)|_{\mu>0} = \mu c'_e B_{\nu}(T_e) + q_0\int_{-1}^0 \mu'I_\nu(0,\mu')d\mu',
$$

$$
  I_\nu(Z,\mu)|_{\mu<0}=\delta(\mu+\mu_s)c_s B_\nu(T_s), \quad B_\nu(T) = \frac{\nu^3}{e^{\frac \nu T}-1}.
$$

$$
\int_0^\infty \sigma_a B_\nu(T)d\nu = \int_0^\infty\sigma_a J_0(z)d\nu.
$$ 

$\mu$ is the cosine of the angle of the ray with the vertical, $z$ is the altitude, $\kappa$ the absorption, $\sigma_s=a\kappa,~\sigma_a=\kappa(1-a)$ where $a$ is the scattering parameter. $T_e,T_s$ are the temperature of Earth and Sun.

The boundary conditions correspond to a Lambert radiation from Earth with albedo $q_0$ and collimated sunlight of cosine direction $\mu_s$.

The grey case is when $\kappa$ and $a$ are constant. In practice they usually are functions of the frequency $\nu$. This archive provides values for $\nu\to\kappa$ from the Gemini website, in the file $\_kappa.txt$.

The integral formulation for $J_0(z)=\frac12\int_{-1}^1 I_\nu d\mu$ is 

$$
	 J_0(z) =  \tfrac{c'_e}2 B_\nu(T_e)E_3(\kappa_\nu z)
	+\tfrac{c_0}2 E_2( \kappa_\nu z)
	+ \tfrac{c_s}2 B_\nu(T_s)\sigma_s e^{-\frac{\kappa_\nu(Z-z)}{\mu_s}}
$$

$$
	+  \tfrac{1}2\int_0^Z  E_1( \kappa_\nu|z'-z|)\left[ \sigma_a B_\nu(T(z')) + \sigma_s J_0(z')\right] d z',
$$
$$
	c_0  =  -q_0 \left(\mu_s c_s B_\nu(T_s) e^{-\frac{\kappa_\nu Z}{\mu_s}}
+ \int_0^Z  E_2(\kappa_\nu z)\left[ \sigma_a B_\nu(T(z)) +  \sigma_s J_0(z)\right]d z\right).
$$

The grey case is the easiest: see

- earthsungrey.cpp

The same with a other quadrature formulae 

- earthsunquadra.cpp

Other programs have variable $\kappa$ read from _kappa.txt.  To obtain the sensitivity of the solution with respect to some parameter use

- earthsungreyAD.cpp  together with the library ddouble.hpp

Finally a cloud which makes $a$ and $\kappa$ altitude dependent has been added. See

- earthsunnucloud.cpp

The easiest way to try the program is to use a terminal and type

g++ earthsungrey.cpp
./a.out

and similar for others. But to work on the program we recommand using an IDE as discussed in the helpxxx.md in the main section of this reporitory.


