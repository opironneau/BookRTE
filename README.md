# Mathematical and Numerical Methods for Radiative Transfer in the Atmosphere — Code Companion
## Introduction
Atmospheric temperature is part of a very complex system, modelled by the equations of meteorology for short-term predictions and by climatology for long-term statistical estimations. The 3 fundamenta[...]

## Chapter 1

The first chapter introduces the problem and some basic numerical methods for the crudest model, Milne's equation.

$$
\mu\partial_\tau I^1 +  I^1 = \frac{1}{2}\int_{-1}^1 I^1 d\mu, 
	\quad I^1(\tau_1,\mu)|_{\mu<0} =I_s\delta(\mu+\mu_s), 
	\quad  I^1(0,\mu)_{\mu>0}  =0
$$

To avoid the Dirac mass, Chandrasekhar proposed to translate the boundary condition and solve 

$$
\mu\partial_\tau I' +  I' = \frac{1}{2}\int_{-1}^1 I' d\mu 
	 + \frac{1}{2} I_s e^{-\frac{\tau_1-\tau}{\mu_s}}, 
	\quad I'(\tau_1,\mu)|_{\mu<0} = 0, 
	\quad  I'(0,\mu)_{\mu>0}  =0.
$$

Both are solved using ISIF (iterations on the source in the integral formulation):
	- Choose $b=0$ or $1$ and $c\ge0$. Set 
  
  $$ J_0^0(\tau)= b\left(\tfrac{1}{2} e^{-\frac{\tau_1-\tau}{\mu_s}}I_s+c\right)$$
	 
Loop in m with
  
 $$ S^m=J_0^m$$
    
 $$ J_0^{m+1} =\frac{1}{2} e^{-\frac{\tau_1-\tau}{\mu_s}}I_s +  \tfrac{1}{2} \int_0^{\tau_1} E_1(|\tau'-\tau|)S^m(\tau')d\tau'$$

where $E_1$ is the first exponential integral function.
