# filepath: milnejax/src/milnejax.py
import jax
import jax.numpy as jnp
import numpy as np
import scipy.special as sc
import matplotlib.pyplot as plt

# Parameters
Nz, Nmu = 30, 100   # nb points in z and tau in (0,1)
kmax = 15           # nb fixed point iterations
tau1 = 2
dtau = tau1 / (Nz - 1)
M = 0.0             # start iterations from above or below
mus = 1.0 / jnp.sqrt(3.0)  # inclination 60°
Is = 1.0            # Sunlight intensity

# Arrays
J00 = jnp.zeros(Nz)
J0 = jnp.zeros(Nz)
S = jnp.zeros(Nz)

# Exponential integral approximation (E1)
def expint_E1(t):
    if t>0 :
        return sc.exp1(t)
    return 0

# Main iteration loop
def run_iterations():
    global J0, J00, S

    J00 = jnp.array([jnp.exp(-i * dtau / mus) * Is / 2 for i in range(Nz)])
    J0 =  J00 + M

    def body_fun(J0):
        S = J0 

        def update_J0(i):
            integral = jnp.sum(jnp.array([
                dtau / 2 * expint_E1(jnp.abs(j - i) * dtau) * S[j]
                for j in range(Nz)
            ]))
            return  J00[i] + integral

        J0_new = jnp.array([update_J0(i) for i in range(Nz)])
        return J0_new

    for k in range(kmax):
        J0 = body_fun(J0) 
    return J0

# Compute I(mu)
def getI(S):
    dmu = 0.051
    muspace = jnp.arange(-1.0, 1.0, dmu)

    def compute(j, mu):
        i_arr = jnp.arange(Nz)
        aux = ((i_arr < j) * (mu < 0) + (i_arr > j) * (mu > 0))
        aux = aux*dtau*jnp.exp(-jnp.abs((i_arr-j)*dtau/mu))*S[i_arr]/jnp.abs(mu)
        return jnp.sum(aux)

    results = [(j * dtau, mu, compute(j, mu)) for mu in muspace for j in range(Nz)]
    return results

if __name__ == "__main__":
    J0 = run_iterations()
    for i in range(Nz):
        print(i * dtau, J0[i])
    
    Ivals = getI(S)
    # optional: print or analyze Ivals

    z = np.linspace(0, tau1, Nz)
    plt.plot(z, np.array(J0))
    plt.xlabel("tau")
    plt.ylabel("J0(tau)")
    plt.title("Mean Radiance J0")
    plt.grid(True)
    plt.show()