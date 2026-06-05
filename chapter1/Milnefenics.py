# Solves Milne's Problem by ISIF
# --------------------------
import numpy as np
from mpi4py import MPI
from dolfinx import mesh
import scipy.special as sc
# -------------------------
# Parameters
# -------------------------
Nz = 200
kmax = 15
tau1 = 1.0
dtau = tau1 / (Nz - 1)
mus = 1 / np.sqrt(3.0)
Is = 1.0

# 1D mesh (only for consistency with FEniCSx)
domain = mesh.create_interval(MPI.COMM_WORLD, Nz - 1, [0.0, tau1])
# -------------------------
# Arrays
# -------------------------
J00 = np.zeros(Nz)
J0  = np.zeros(Nz)
S   = np.zeros(Nz)
# -------------------------
# Exponential integral E1 from scipy.special
# -------------------------
def expint_E1(t):
    if t>0 :
        return sc.exp1(t)
    return 0
# -------------------------
# Initial condition
# -------------------------
for i in range(Nz):
    J00[i] = np.exp(-i * dtau / mus) * Is / 2
    J0[i] = J00[i]
# -------------------------
# Fixed-point iteration
# -------------------------
for k in range(kmax):

    S[:] = J0[:]

    for i in range(Nz):
        J0[i] = J00[i]
        for j in range(Nz):
            J0[i] += (
                dtau / 2
                * expint_E1((j - i + 0.5) * dtau)
                * S[j]
            )

    if MPI.COMM_WORLD.rank == 0:
        print(k, J0[Nz // 2])
# -------------------------
# Output result
# -------------------------
if MPI.COMM_WORLD.rank == 0:
    print()
    for i in range(Nz):
        print(i * dtau, J0[i])