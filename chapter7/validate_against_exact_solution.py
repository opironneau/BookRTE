import numpy as np
import dg_transport_solver as dg

# Override problem data for the manufactured-solution test:
# n(z) = const  => dlogn/dz = 0 ;  kappa = kappa0 = const ; f = 1
kappa0 = 0.5

def dlogn_dz_zero(z):
    z = np.asarray(z, dtype=float)
    return 0.0*z

def kappa_const(z):
    z = np.asarray(z, dtype=float)
    return kappa0 + 0.0*z

dg.dlogn_dz = dlogn_dz_zero
dg.kappa_func = kappa_const

Zdom = 1.2
mesh, space, U = dg.solve(p_deg=3, n_mu=20, n_z=20, Z=Zdom, verbose=True)

def exact(mu, z, Z):
    mu = np.asarray(mu,dtype=float)
    z = np.asarray(z,dtype=float)
    out = np.zeros_like(mu)
    pos = mu>0
    neg = mu<0
    out[pos] = 1.0/kappa0 + (mu[pos]-1.0/kappa0)*np.exp(-kappa0*z[pos]/mu[pos])
    out[neg] = 1.0/kappa0 + (mu[neg]-1.0/kappa0)*np.exp(-kappa0*(z[neg]-Z)/mu[neg])
    out[mu==0] = 1.0/kappa0
    return out

# sample many interior points away from mu=0 (weak singularity direction) and away from very edges
rng = np.random.default_rng(0)
mus = rng.uniform(-1,1,4000)
zs = rng.uniform(0,Zdom,4000)
mask = np.abs(mus) > 0.05   # avoid near-singular mu=0 region for the ODE exponential
mus, zs = mus[mask], zs[mask]

u_num = dg.eval_solution(mesh, space, U, mus, zs)
u_ex = exact(mus, zs, Zdom)

err = np.abs(u_num-u_ex)
print("max abs error         :", err.max())
print("mean abs error         :", err.mean())
print("max |u_exact|          :", np.max(np.abs(u_ex)))
# worst points
worst = np.argsort(-err)[:5]
for i in worst:
    print(f" mu={mus[i]:+.4f} z={zs[i]:+.4f}  u_num={u_num[i]:+.5f} u_ex={u_ex[i]:+.5f} err={err[i]:.2e}")
