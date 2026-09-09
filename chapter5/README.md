# CHAPTER 5
## Discontinuous Refractive Index

The VRTE of chapter 4 are valid for a piecewise constant refractive index $n$ if $I,Q$ are replaced by $\tilde I=I/n^2, \tilde Q=Q/n^2$. The programs deal with the case $n=n_-$ in $(0,Y)$ and and $n=n^+$ in $(Y,Z)$.  The difficulty is that the Fresnel transmission conditions at $z=Y$ are quite complex. They are implemented in snell2.edp (FreeFem++)

ISIF with the Fresnel conditions is implemented in IQn.cpp.  It is possible to simulate system with discontinuous density at $z=Y$.

The programs of Chapter 6 can also apply to the piecewise constant case.
