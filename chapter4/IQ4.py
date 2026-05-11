import os # use Xcode python on Mac+VS Code
import time
import jax
import jax.numpy as jnp
from jax.scipy.special import expn
import numpy as np
from functools import partial


# Configuration pour utiliser des flottants 64 bits (similaire au 'double' en C++)
jax.config.update("jax_enable_x64", True)

# --- Constantes ---
pi = jnp.pi
stefan = (pi**4) / 15

Ce = 2.0
Cs = 2.0e-6
rho0 = 1.0
drho = -0.7
TOA = 1.2
Z = TOA * rho0 * (1 + drho * TOA / 2)

B0 = 1.4744e-8
T0 = 4798.0

Te = (273 + 18) / T0
Ts = 5778 / T0
q0 = -0.3
mus = 0.5
beta = 0.5

verbose = True
Nz = 10
dz = Z / (Nz - 1)
kmax = 14
newton = 50
epsdycho = 1e-4
epsnewton = 1.e-10
kappamin = 0.001
z1 = Z * 0.5
z2 = Z * 0.7
z3 = 0.8 * Z
nu1 = 0.5
nu2 = 1.0
cloud = 1.0

# --- Chemins ---
basedir = "/Users/pironneau/Dropbox/aranger/TeX2026/BookVRTE/prog4/"
mykappafile = os.path.join(basedir, "_kappa.txt")
myresulttemperature = os.path.join(basedir, "temperatureyyx")

# --- Fonctions Physiques ---
#@jax.jit
def kapparhoz(nu, z):
    return 1.0 + cloud * ((z > z1) & (z < z2)) * ((nu > 3/10.) & (nu < 3/1.) + (nu < 3/20.))

#@jax.jit
def as_scatter(z, nu):
    return 0.3 + \
           0.3 * 4 * (z2 - z) * (z - z1) * (z > z1) * (z < z2) / ((z1 - z2)**2) + \
           0.3 * 16 * (nu < nu2) * (nu > nu1) * (((nu - nu1) * (nu - nu2))**2 / ((nu1 - nu2)**2)) * (z > z3)

#@jax.jit
def Bsun(nu):
    return (nu**3) / (jnp.exp(nu / Ts) - 1.0)

#@jax.jit
def Bearth(nu):
    return (nu**3) / (jnp.exp(nu / Te) - 1.0)

#@jax.jit
def BB(nu, T):
    return (nu**3) / (jnp.exp(jnp.fmin(nu / T, 100.0)) - 1.0)

#@jax.jit
def dBB(nu, T):
    a = jnp.exp(nu / T)
    # Condition pour éviter les dépassements (overflow)
    return jax.lax.cond(
        a > 1e100,
        lambda _: ((nu**2) / T) * 1e-30,
        lambda _: a * ((nu**2) / jnp.fmax(1.e-30, a - 1) / T)**2,
        operand=None
    )


# --- Initialisation des données ---
def readkappa(filepath):
    print(f"reading kappa(nu) from file {filepath}")
    if not os.path.exists(filepath):
        raise RuntimeError(f"Cannot open file: {filepath}")
    
    # Lecture avec numpy standard
    data = np.loadtxt(filepath)
    waveno = data[:, 0]
    kappaux = data[:, 1]
    jmax = len(waveno)
    kappanu_np = np.maximum(kappaux, kappamin)
    nu_np = 3.0 / waveno
    
    kappasum_np = np.zeros((jmax, Nz))
    kappaz_np = np.zeros((jmax, Nz))
    
    z_vals = np.arange(Nz) * dz
    
    for j in range(jmax):
        kappaz_np[j, :] = kapparhoz(nu_np[j], z_vals)
        # Intégration trapézoïdale cumulée
        kappasum_np[j, 0] = 0
        for it in range(1, Nz):
            kappasum_np[j, it] = kappasum_np[j, it-1] + (kappaz_np[j, it] + kappaz_np[j, it-1]) * dz / 2.0
            
    print(f"Number of frequencies {jmax}")
    return jmax, jnp.array(kappanu_np), jnp.array(nu_np), jnp.array(kappaz_np), jnp.array(kappasum_np)

# --- Cœur de l'algorithme ---

#@partial(jax.jit, static_argnums=0)
def updateJK_step(jmax, kappanu, nu, kappaz, kappasum, J0old, J2old, K0old, K2old, T):
    J0 = jnp.zeros_like(J0old)
    J2 = jnp.zeros_like(J2old)
    K0 = jnp.zeros_like(K0old)
    K2 = jnp.zeros_like(K2old)
    
    z_arr = jnp.arange(Nz) * dz

    for jnu in range(jmax):
        kappanuj = kappanu[jnu]
        nuj = nu[jnu]
        
        c0 = -q0 * mus * Cs * Bsun(nuj) * jnp.exp(-kappanuj * kappasum[jnu, Nz-1] / mus)
        
        kappazj = kappaz[jnu, :] * kappanuj
        asz = jax.vmap(as_scatter, in_axes=(0, None))(z_arr, nuj)
        sigs = kappazj * asz
        siga = kappazj * (1 - asz)
        
        S = sigs * J0old[jnu, :] + siga * jax.vmap(BB, in_axes=(None, 0))(nuj, T)
        H = 9 * beta * sigs / 8 * (J2old[jnu, :] - J0old[jnu, :] / 3 - K0old[jnu, :] + K2old[jnu, :])
        
        # Intégrale pour ajuster c0
        E2_vals = expn(2, kappasum[jnu, :] * kappanuj)
        c0 -= q0 * jnp.sum(E2_vals * S * dz)

        # Calcul principal sur l'altitude
        for i in range(Nz):
            ksz = kappasum[jnu, i] * kappanuj
            
            J0z = Ce * Bearth(nuj) * expn(3, ksz) / 2 + \
                  Cs * Bsun(nuj) * jnp.exp(-(kappasum[jnu, Nz-1] - kappasum[jnu, i]) * kappanuj / mus) / 2 + \
                  c0 * expn(2, ksz) / 2
                  
            J2z = Ce * Bearth(nuj) * expn(5, ksz) / 2 + \
                  Cs * (mus**2) * Bsun(nuj) * jnp.exp(-(kappasum[jnu, Nz-1] - kappasum[jnu, i]) * kappanuj / mus) / 2 + \
                  c0 * expn(4, ksz) / 2
                  
            K0z = 0.0
            K2z = 0.0
            
            # Convolution
            for j in range(1, Nz):
                Hj = (H[j] + H[j-1]) / 2
                Sj = (S[j] + S[j-1]) / 2 - Hj / 3
                aux = kappanuj * (kappasum[jnu, i] - (kappasum[jnu, j] + kappasum[jnu, j-1]) / 2)
                
                expE1ij = expn(1, jnp.abs(aux))
                expE3ij = expn(3, jnp.abs(aux))
                expE5ij = expn(5, jnp.abs(aux))
                
                # JAX expn retourne 0 pour x très grand, C++ gérait un seuil.
                J0z += (expE1ij * Sj + expE3ij * Hj) * dz / 2
                J2z += (expE3ij * Sj + expE5ij * Hj) * dz / 2
                K0z -= (expE1ij - expE3ij) * Hj * dz / 2
                K2z -= (expE3ij - expE5ij) * Hj * dz / 2
                
            J0 = J0.at[jnu, i].set(J0z)
            J2 = J2.at[jnu, i].set(J2z)
            K0 = K0.at[jnu, i].set(K0z)
            K2 = K2.at[jnu, i].set(K2z)
            
    return J0, J2, K0, K2

def root(rhs, T0_val, i, jmax, nu, kappanu):
    myeq = -rhs
    for j in range(1, jmax):
        nu_mid = (nu[j] + nu[j-1]) / 2
        dnu = nu[j] - nu[j-1]
        myeq += (1 - as_scatter(i * dz, nu[j])) * kappanu[j] * BB(nu_mid, T0_val) * dnu
    return myeq

def rhsTeq(i, jmax, nu, kappanu, J0):
    rhs = 0.0
    for j in range(1, jmax):
        rhs += (1 - as_scatter(i * dz, nu[j])) * kappanu[j] * J0[j, i] * (nu[j] - nu[j-1])
    return rhs

def getTbydycho(rhs, i, Tstart, jmax, nu, kappanu):
    if rhs == 0: return 0.0
    T0_val = Tstart
    T1_val = 2 * T0_val
    counter = 0
    
    myeq0 = root(rhs, T0_val, i, jmax, nu, kappanu)
    myeq1 = root(rhs, T1_val, i, jmax, nu, kappanu)
    
    while myeq0 > 0 and counter < 100:
        T0_val /= 2
        counter += 1
        myeq0 = root(rhs, T0_val, i, jmax, nu, kappanu)
        
    while myeq1 < 0 and counter < 100:
        T1_val *= 2
        myeq1 = root(rhs, T1_val, i, jmax, nu, kappanu)
        counter += 1
        
    while (T1_val - T0_val) > epsdycho and counter < 100:
        Taux = (T1_val + T0_val) / 2
        counter += 1
        myeq0 = root(rhs, Taux, i, jmax, nu, kappanu)
        if myeq0 > 0:
            T1_val = Taux
        else:
            T0_val = Taux
            
    return (T1_val + T0_val) / 2

def genT(jmax, nu, kappanu, J0, T_old):
    T_new = np.array(T_old)
    for i in range(Nz):
        rhs = rhsTeq(i, jmax, nu, kappanu, J0)
        T_new[i] = getTbydycho(rhs, i, T_new[i], jmax, nu, kappanu)
        
        presfunc = 1.0
        inewton = 0
        while inewton < newton and abs(presfunc) > epsnewton:
            inewton += 1
            T0_val = T_new[i]
            left = 0.0
            deriv = 0.0
            
            for j in range(1, jmax):
                dnu = nu[j] - nu[j-1]
                nu11 = (nu[j] + nu[j-1]) / 2
                kappaux = (1 - as_scatter(i * dz, nu[j])) * (kappanu[j] + kappanu[j-1]) / 2
                left += kappaux * BB(nu11, T0_val) * dnu
                deriv += kappaux * dBB(nu11, T0_val) * dnu
                
            presfunc = rhs - left
            if abs(deriv) > 1e-10:
                T_new[i] = T0_val + presfunc / deriv
                
    return jnp.array(T_new)

def backtoz(z_val):
    if drho == 0:
        return z_val
    return (jnp.sqrt(1 + 2 * drho * z_val / rho0) - 1) / drho

def multiBlock(initT, jmax, kappanu, nu, kappaz, kappasum):
    T = jnp.full(Nz, initT)
    J0old = jnp.zeros((jmax, Nz))
    J2old = jnp.zeros((jmax, Nz))
    K0old = jnp.zeros((jmax, Nz))
    K2old = jnp.zeros((jmax, Nz))

    for k in range(kmax):
        J0, J2, K0, K2 = updateJK_step(jmax, kappanu, nu, kappaz, kappasum, J0old, J2old, K0old, K2old, T)
        T = genT(jmax, nu, kappanu, J0, T)
        
        J0old, J2old, K0old, K2old = J0, J2, K0, K2
        
        normG = 0.0
        for j in range(jmax):
            normG += jnp.sum(jnp.abs(J0[j, :]) * (nu[j] - nu[j-1]) * dz)
            
        print(f"k= {k}  {(T[2] * 4798) - 273:.4f}  {normG:.6e}")
        
    return T

# --- Main ---
def main():
    try:
        jmax, kappanu_base, nu_base, kappaz, kappasum = readkappa(mykappafile)
    except Exception as e:
        print(f"Erreur lors de la lecture: {e}")
        return

    for K in range(2):
        kappanu = np.array(kappanu_base)
        
        for j in range(jmax):
            if K == 1:
                if (nu_base[j] > 3./18) and (nu_base[j] < 3./14):
                    kappanu[j] = max(0.5, 1.2 * kappanu_base[j])
            elif K == 2:
                kappanu[j] = 0.5
                
        kappanu = jnp.array(kappanu)

        output_file = f"{myresulttemperature}{K}.txt"
        print(f"\n iterations \t [T] near earth [T] far ||G|| and ||S||")
        
        t0 = time.time()
        T_final = multiBlock(Te / 2, jmax, kappanu, nu_base, kappaz, kappasum)
        cpu_time = time.time() - t0
        
        print(f" Time CPU = {cpu_time:10.6f}\n")
        print(" tau\t [T]:")
        
        with open(output_file, 'w') as f:
            for i in range(Nz - 1):
                z_val = backtoz(i * dz)
                t_val = T_final[i] * T0 - 273
                print(f"{z_val:.6f} {t_val:.6f}")
                f.write(f"{z_val:.6f} {t_val:.6f}\n")

if __name__ == "__main__":
    main()