MilneMC.cpp was written by Claude AI.  To use it do the following

1. in Terminal do a cd to this folder

2. Check in MilneMC.cpp that the data are suitable to you

3. Compile MineMC.cpp as said in the comments of MilneMC.cpp

4. Run MilneMC

5. in the terminal type gnuplot milneMC_I.gp

6. Look at milneMC.png

**Method.** Photons are born at $z=0$ with $\mu=\mu_s$, carrying weight $w_0=\mu_s/N$ (the incoming current $\int_0^1\mu I(0,\mu)d\mu=\mu_s$). Free paths are sampled as $s=-\log U/\kappa$ with $z\!\to\!z+\mu s$; scattering is isotropic ($\mu'=2U-1$), as dictated by the $\frac12\int I\,d\mu'$ kernel. Absorption uses implicit capture, $w\to a w$ with albedo $a=\sigma_s/\kappa$ (so $a=1$ for the equation exactly as you wrote it), plus Russian roulette. Estimators:

- $\int_{-1}^{1} I\,d\mu = 2J_0$ and $I(z,\mu)$ by track length: $\sum w\ell/\Delta z$ (and $/\Delta\mu$ for the angular bin);
- emergent $I(Z,\mu),\ \mu>0$ and $I(0,\mu),\ \mu<0$ by a next-event estimator at each collision, $\frac{aw}{2\mu}e^{-\kappa(Z-z)/\mu}$ — much smoother than counting escapes;
- the direct beam is *not* sampled: the first flight is excluded from the tallies and $I_{\rm dir}=\delta(\mu-\mu_s)e^{-\kappa z/\mu_s}$ is added analytically, which removes the $\delta$ spike from the angular bins.

Error bars come from batch means (default 20 batches). Parameters are `key=value` on the command line: `Z kappa sigmas mus N Nz Nmu batches seed`. Output goes to `milneMC_J0.txt`, `milneMC_I.txt`, `milneMC_exit.txt`.

**Verification.** Three checks, all passing:

1. $\sigma_s=0$: $J_0$ reproduces $\tfrac12 e^{-\kappa z/\mu_s}$ exactly, and the transmitted current matches $\mu_s e^{-\kappa Z/\mu_s}=0.10215$ (analog tally $0.10227\pm0.00044$).
2. Conservative case ($\kappa=\sigma_s=1$, $Z=1$, $\mu_s=1/\sqrt3$): I solved the equivalent integral equation $J_0(z)=\tfrac12 e^{-\kappa z/\mu_s}+\tfrac{a}{2}\int_0^Z \kappa E_1(\kappa|z-z'|)J_0\,dz'$ to high accuracy (cell-exact integration of $E_1$ via $E_2$). MC agrees within $\pm0.1\%$, i.e. inside the statistical error, for both $J_0(z)$ and $I(Z,\mu)$. Energy balance in/out $=1.0000$.
3. Absorbing case ($\kappa=2,\ \sigma_s=0.8,\ Z=1.5,\ \mu_s=0.6$): same agreement on $J_0$ and on both emergent intensities.

For reference, [chapter1/MilneProblem.cpp](chapter1/MilneProblem.cpp) sits about 0.4–1.5% below the converged answer — that is its midpoint quadrature of the log-singular $E_1$ kernel (and `kmax=15` is short of convergence for $a=1$), not a defect in the Monte Carlo.

[MilneMC.cpp](chapter4/MilneMC.cpp) now writes a gnuplot script and calls it, so a run produces the figure directly:

```bash
g++ -O3 -std=c++17 -o MilneMC MilneMC.cpp && ./MilneMC N=4000000 Nmu=48
```

- `milneMC_I.txt` — columns `z  mu  I_diffuse  I_diffuse+beam`, with a blank line between the $z$ scans so gnuplot reads it as a grid.
- `milneMC_I.gp` — the script (editable), `milneMC_I.png` — the figure: left, the pm3d surface with contours projected on the base (view `55,340`, chosen so both vacuum edges are visible); right, the same seen from above with level lines.
- `plot=0` on the command line skips the gnuplot call and only writes the script.

The plotted quantity is the diffuse part, since the direct beam is a $\delta(\mu-\mu_s)$. Setting `col = 4` at the top of the script instead of `col = 3` adds it back, spread over the $\mu$ bin containing $\mu_s$ (i.e. $e^{-\kappa z/\mu_s}/\Delta\mu$ in that bin), which is what a raw track-length tally would have shown.

The surface reads as expected: $I\to0$ on the two inflow edges where the diffuse field has no source, $z=0$ with $\mu>0$ and $z=Z$ with $\mu<0$; the ridge near $\mu\simeq0^-$ at small $z$ is the grazing radiation backscattered out of the illuminated face, and the edge values agree with the emergent intensities computed by the next-event estimator ($I(0,-0.975)=0.4217$, $I(Z,0.975)=0.3618$).


Let me check the numbers before changing anything, because I think that ridge is physical.

The numbers say the plot is right — but let me show you why, and then make the figure say it out loud.

The boundary condition is fine — two different things were confusing the picture, and only the first was my doing.

**The Dirac at $\mu_s$ was deliberately absent.** The plotted column is the *diffuse* field: I remove the direct beam from the tallies and carry it analytically as $\delta(\mu-\mu_s)e^{-\kappa z/\mu_s}$, precisely because a $\delta$ cannot be drawn on a surface (smeared over one bin its height is $1/\Delta\mu\approx24$, and everything else becomes flat). That was a silent choice on my part. The figure now draws it explicitly: red impulses standing at $\mu=\mu_s=1/\sqrt3$, of height $e^{-\kappa z/\mu_s}$ — the weight of the Dirac. `col = 4` in the script still puts it back inside the surface if you want it that way.

**The ridge at $\mu\to0^-$ is not a Dirac, it is the grazing limit.** For the diffuse field
$$I(z,\mu)=\frac1{|\mu|}\int S(z')e^{-\kappa|z-z'|/|\mu|}dz'\ \xrightarrow[\mu\to0]{}\ \frac{S(z)}{\kappa}=\frac{\sigma_s J_0(z)}{\kappa},$$
because the exponential concentrates on the local source. So along $\mu=0$ the surface must ride at $\sigma_s J_0/\kappa$ — that is the black crest line I added, and the computed crest sits on it. At the corner $z=0$ this collides with $I(0,\mu)=0$ for all $\mu>0$ (no diffuse inflow): the exact solution genuinely jumps from $0$ at $\mu=0^+$ to $\approx0.81$ at $\mu=0^-$. That cliff is what looked like a spike.

It is a smooth monotone rise across the 48 bins, not a peak, and it matches the integral equation:

| $\mu$ | $-0.979$ | $-0.500$ | $-0.104$ | $-0.062$ | $-0.021$ |
|---|---|---|---|---|---|
| MC $I(0,\mu)$ | 0.4205 | 0.6020 | 0.8067 | 0.8175 | 0.8199 |
| exact | 0.4208 | 0.6020 | 0.8067 | 0.8179 | 0.8194 |

with $\sigma_sJ_0(0)/\kappa = 0.807$ as the $\mu\to0^-$ asymptote. Physically it is the radiation that scattered once or twice just under the illuminated face and creeps out nearly parallel to it — the same effect that makes limb darkening, and the reason $I(0,\mu)$ doubles between normal and grazing exit.



\section*{Appendix}
\subsection{Continuity of Solutions of the Characteristic Equation}

\begin{proposition}
Let $\kappa>0$, $I_0,I_Z>0$, and $S\in C^0([0,Z])$, $S>0$. Let $I$ defined by
\begin{align}
	I (z,\mu) &= \One_{\mu>0}\left[ \e^{-\frac1\mu \kappa  z }I_0\mu
	+ \int_{0}^z \frac1\mu\e^{-\frac1\mu \kappa (z-z')}S(z')\d z'\right]
\cr&
I(z,0) =\frac{S(z)}\kappa, \hbox{consistently with the transport equation
$\mu\partial_z I+\kappa I=S$ at $\mu=0$.}
\cr&
	+ \One_{\mu<0}\left[ \e^{\frac1\mu \kappa (Z-z) }I_Z
	- \int_{z}^Z \frac1\mu\e^{\frac1\mu \kappa (z'-z)}S(z')\d z'\right].
\end{align}
$z,\mu \mapsto I$ is continuous in $(0,Z)\times[-1,1]$.
\end{proposition}
%
\begin{proof}

Since $S$ is continuous on the compact $[0,Z]$, it is bounded, $M:=\max_{[0,Z]}S$,
and uniformly continuous. Continuity for $\mu\neq0$ is clear; the point is the
behaviour as $\mu\to0$.

\medskip
\noindent\textbf{Step 1 (rewriting, $\mu>0$).}
The change of variable $u=(z-z')/\mu$ gives
\[
\int_{0}^z \frac1\mu\,\e^{-\frac{\kappa}{\mu}(z-z')}S(z')\d z'
=\int_0^{z/\mu}\e^{-\kappa u}\,S(z-\mu u)\d u ,
\]
and since $S(z)/\kappa=S(z)\int_0^\infty\e^{-\kappa u}\d u$,
\begin{equation}\label{eq:split}
I(z,\mu)-\frac{S(z)}{\kappa}
=\underbrace{I_0\,\mu\,\e^{-\kappa z/\mu}}_{T_1}
+\underbrace{\int_0^{z/\mu}\e^{-\kappa u}\big[S(z-\mu u)-S(z)\big]\d u}_{T_2}
-\underbrace{S(z)\int_{z/\mu}^{\infty}\e^{-\kappa u}\d u}_{T_3}.
\end{equation}

\medskip
\noindent\textbf{Step 2 (estimates).}
Clearly
\[
|T_1|\le I_0\,\mu,
\qquad
|T_3|\le \frac{M}{\kappa}\,\e^{-\kappa z/\mu}.
\]
For $T_2$, let $\varepsilon>0$ and choose $\delta>0$ such that
$|S(a)-S(b)|<\varepsilon$ whenever $|a-b|<\delta$ (uniform continuity).
Splitting the integral at $u=\delta/\mu$,
\[
|T_2|\le \varepsilon\int_0^{\delta/\mu}\e^{-\kappa u}\d u
+2M\int_{\delta/\mu}^{\infty}\e^{-\kappa u}\d u
\le \frac{\varepsilon}{\kappa}+\frac{2M}{\kappa}\,\e^{-\kappa\delta/\mu}.
\]
Hence, for every compact set $[a,Z-a]\subset(0,Z)$,
\begin{equation}\label{eq:unif}
\sup_{z\in[a,Z-a]}\Big|I(z,\mu)-\frac{S(z)}{\kappa}\Big|
\le I_0\,\mu+\frac{M}{\kappa}\,\e^{-\kappa a/\mu}
+\frac{\varepsilon}{\kappa}+\frac{2M}{\kappa}\,\e^{-\kappa\delta/\mu},
\end{equation}
so $\limsup_{\mu\to0^+}\sup_{[a,Z-a]}|I(\cdot,\mu)-S/\kappa|\le\varepsilon/\kappa$
for every $\varepsilon>0$, i.e.\ the supremum tends to $0$.
(Alternatively, $T_2\to0$ follows from dominated convergence: the integrand tends
to $0$ pointwise and is dominated by $2M\e^{-\kappa u}\in L^1(0,\infty)$.)

The case $\mu<0$ is symmetric: writing $\mu=-|\mu|$ and $u=(z'-z)/|\mu|$,
\[
I(z,\mu)=I_Z\,\e^{-\kappa(Z-z)/|\mu|}
+\int_0^{(Z-z)/|\mu|}\e^{-\kappa u}\,S(z+|\mu| u)\d u ,
\]
and the same three estimates hold with $z$ replaced by $Z-z$ and $I_0$ by $I_Z$.
Therefore \eqref{eq:unif} holds with $\mu$ replaced by $|\mu|$, uniformly for
$z\in[a,Z-a]$.

\medskip
\noindent\textbf{Step 3 (joint continuity).}
Let $(z_n,\mu_n)\to(z_0,0)$ with $z_0\in(0,Z)$; choose $a>0$ with
$z_n,z_0\in[a,Z-a]$ for $n$ large. Then
\[
\Big|I(z_n,\mu_n)-\frac{S(z_0)}{\kappa}\Big|
\le \Big|I(z_n,\mu_n)-\frac{S(z_n)}{\kappa}\Big|
+\frac{|S(z_n)-S(z_0)|}{\kappa}\longrightarrow 0,
\]
the first term by the uniform estimate \eqref{eq:unif}, the second by continuity
of $S$. Together with continuity for $\mu\neq0$, this proves the claim.
\end{proof}

\begin{proof}[Remark]
\renewcommand{\qedsymbol}{}
Continuity fails at the corners $(0,0)$ and $(Z,0)$: for instance
$I(0,\mu)=I_0\mu\to0$ as $\mu\to0^+$ while $I(0,\mu)\to S(0)/\kappa$ as
$\mu\to0^-$. This is reflected in \eqref{eq:unif}, where the bound degrades as
$a\to0$: the terms $\e^{-\kappa a/\mu}$ are not small when $a\lesssim\mu/\kappa$,
i.e.\ there are boundary layers of width $O(|\mu|/\kappa)$ near $z=0$ and $z=Z$.
\end{proof}
