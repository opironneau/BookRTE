# Solve   mu du/dz + dlogn(z)(1-mu^2) du/dmu + kappa(z) u = 1   in (-1,1)x(0,Z)
# with inflow condition u = mu on Sigma = {z=0, mu>=0} U {z=Z, mu<=0}.
# First order upwind finite differences on a tensor grid, sparse direct solve.
# No condition is needed at mu = +-1 because the mu-advection vanishes there.

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

Z = 1.0
aref = 0.2                                   # refraction index n(z) = 1 + aref*z
kappa0 = 0.5

def dlogn(z): return aref / (1 + aref * z)   # d/dz log n
def kappa(z): return kappa0 * (1 - 0.7 * z)  # altitude effect on kappa


def solve(Nmu, Nz, dlogn, kappa, Z):
    mu = np.linspace(-1, 1, Nmu + 1)
    z = np.linspace(0, Z, Nz + 1)
    dmu, dz = 2.0 / Nmu, Z / Nz
    MU, ZG = np.meshgrid(mu, z, indexing='ij')   # shape (Nmu+1, Nz+1), j fastest
    K = np.arange(MU.size).reshape(MU.shape)

    dirich = ((ZG == 0) & (MU >= 0)) | ((ZG == Z) & (MU <= 0))
    free = ~dirich
    rows, cols, vals = [], [], []

    def add(mask, col, val):
        rows.append(K[mask]); cols.append(col[mask]); vals.append(val[mask])

    add(free, K, kappa(ZG) + 0 * MU)             # absorption

    m = free & (MU > 0)                          # z-advection, upwind = below
    add(m, K, MU / dz); add(m, K - 1, -MU / dz)
    m = free & (MU < 0)                          # z-advection, upwind = above
    add(m, K, -MU / dz); add(m, K + 1, MU / dz)

    c = dlogn(ZG) * (1 - MU**2)                  # mu-advection (zero at mu=+-1)
    m = free & (c > 0)
    add(m, K, c / dmu); add(m, K - (Nz + 1), -c / dmu)
    m = free & (c < 0)
    add(m, K, -c / dmu); add(m, K + (Nz + 1), c / dmu)

    add(dirich, K, 1 + 0 * MU)                   # Dirichlet rows
    b = np.where(dirich, MU, 1.0).ravel()

    A = sp.coo_matrix((np.concatenate(vals),
                       (np.concatenate(rows), np.concatenate(cols))),
                      shape=(MU.size, MU.size)).tocsr()
    U = spla.spsolve(A, b).reshape(MU.shape)
    return mu, z, U


def exact_flat(mu, z, kap, Z):
    # exact solution when dlogn = 0 and kappa constant
    MU, ZG = np.meshgrid(mu, z, indexing='ij')
    with np.errstate(divide='ignore', invalid='ignore'):
        U = np.where(MU > 0, 1/kap + (MU - 1/kap) * np.exp(-kap * ZG / MU),
            np.where(MU < 0, 1/kap + (MU - 1/kap) * np.exp(-kap * (ZG - Z) / MU),
                     1/kap))
    return U


if __name__ == "__main__":
    # validation on the flat case dlogn=0, kappa constant, against the exact solution
    for N in (50, 100, 200):
        mu, z, U = solve(N, N, lambda z: 0 * z, lambda z: kappa0 + 0 * z, Z)
        err = np.max(np.abs(U - exact_flat(mu, z, kappa0, Z)))
        print(f"flat case  N={N:4d}   max error = {err:.4f}")

    # refractive case
    Nmu = Nz = 200
    mu, z, U = solve(Nmu, Nz, dlogn, kappa, Z)

    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    pc = ax[0].pcolormesh(mu, z, U.T, shading='gouraud', cmap='turbo')
    fig.colorbar(pc, ax=ax[0])
    ax[0].set_xlabel(r'$\mu$'); ax[0].set_ylabel('z')
    ax[0].set_title(rf'$u(\mu,z)$,  $n=1+{aref}z$, $\kappa={kappa0}(1-0.7z)$')
    for zc in (0.0, 0.25, 0.5, 0.75, 1.0):
        j = int(round(zc / Z * Nz))
        ax[1].plot(mu, U[:, j], label=f'z={z[j]:.2f}')
    ax[1].set_xlabel(r'$\mu$'); ax[1].set_ylabel('u'); ax[1].legend()
    ax[1].set_title(r'sections $u(\cdot,z)$')
    plt.tight_layout()
    plt.savefig('refractive_upwind.png', dpi=150)
    plt.show()
