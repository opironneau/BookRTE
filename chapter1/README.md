# Programs from CHAPTER 1 

Except files lake4book.xxx, all files are about solving Milne's problem.

## The Milne problem

Of interest to us is

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
To compute $I(x,\mu)$ at some given $x,\mu$ one must approximate the Dirac mass at $\mu_s$ by an exponential and use the following: if $\mu> > 0$,

$$
I(x,\mu)\approx \int_0^x \frac1\mu e^{-\frac {|x-x'|}\mu}J_0(x')d x' +  c_s  e^{-\frac x{\mu_s}} e^{-\lambda(\mu-\mu_s)}
/\int_0^1 e^{-\lambda(\mu-\mu_s)}d\mu.
$$

When $\mu<0$,

$$
I(x,\mu)= -\int_x^1 \frac1\mu e^{\frac {|x-x'|}\mu}J_0(x')d x'.
$$

