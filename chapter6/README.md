# CHAPTER 6 
## Graded diffractive index

When the refractive index is a function of altitude $z$, VRTE for the stratified case is

$$
{\mathcal L}(\tilde I) =\sigma_a \tilde B_\nu + \frac{\sigma_s}2\int_{-1}^1 \tilde I d\mu'
+ \frac{\beta\sigma_s}4 P_2(\mu)\int_{-1}^1 [P_2\tilde I-(1-P_2 )\tilde Q] d\mu',
$$

$$
{\mathcal L}\tilde Q) = -\frac{\beta\sigma_s}4 (1-P_2(\mu))\int_{-1}^1 [P_2\tilde I-(1-P_2 )\tilde Q] d\mu',
$$
where $P_2(\mu)=\tfrac12(3\mu^2-1)$ and

$$
{\mathcal L}(I):=\mu\partial_z I + \partial_z\log n (1-\mu^2)\partial_\mu I
+\kappa I.
$$

If $n$ is discontinuous at $z=Y$ then the Fresnel transmission conditions must be applied.

All programs in this folder are FreeFem++ scripts.

The grey case without discontinuity are treated by 

- greygradedindex.edp, where $\kappa_\nu$ is constant.
  
- GreyPolarizedSmooRefracttiveSUPG.md is the same but it has embedded math comments for clarity.
  
- nongreygeneral.edp  solve the case of a $\nu$ depended $\kappa_\nu$.
  
- VRTEFresnel allows for a discontinuity at $z=Y$.

