import os
import time
import jax
import jax.numpy as jnp
from jax.scipy.special import expn
import numpy as np

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

verbose = 1
Nz = 50
dz = Z / (Nz - 1)
kmax = 9
newton = 50
epsdycho = 1e-4
epsnewton = 1.e-10
kappamin = 0.001

def zprime(z):
    return z*rho0*(1+drho*z/2)

z1 = zprime(0.4)
z2 = zprime(0.7)
z3 = zprime(0.8)
nu1 = 0.5
nu2 = 1.0
cloud = 1.5

# --- Chemins ---
basedir = ""
mykappafile = os.path.join(basedir, "_kappa.txt")
myresulttemperature = os.path.join(basedir, "temperatureyyx")

# --- Fonctions Physiques ---
# These work with scalar or array inputs via broadcasting
def kapparhoz(nu, z):
    return 1.0 + cloud * ((z > z1) & (z < z2)) * ((nu > 3/10.) * (nu < 3/1.) + (nu < 3/20.))

def as_scatter(z, nu):
    return 0.3 + \
           cloud * 0.3 * 4 * (z2 - z) * (z - z1) * (z > z1) * (z < z2) / ((z1 - z2)**2) + \
           0.3 * 16 * (nu < nu2) * (nu > nu1) * (((nu - nu1) * (nu - nu2))**2 / ((nu1 - nu2)**4)) * (z > z3)

def expint_E1(t):
    # Matches C++ exactly: small→0, large→asymptotic, mid→series with Kexpint=9+(t-1)*4 terms
    t1 = np.abs(np.asarray(t, dtype=float))
    scalar = t1.ndim == 0
    t1 = np.atleast_1d(t1)
    small = t1 < 1e-5
    large = t1 > 4.0
    series_t = np.where(small | large, 1.0, t1)
    # C++ uses Kexpint = 9+(t1-1)*4; max over array elements to cover all
    Kexpint = int(9 + (min(float(np.max(np.where(small | large, 0.0, t1))), 4.0) - 1) * 4)
    Kexpint = max(Kexpint, 3)
    ak = series_t.copy()
    soNtaue = -0.577215664901533 - np.log(series_t) + ak
    for k in range(2, Kexpint):
        ak = ak * (-series_t * (k - 1) / (k * k))
        soNtaue = soNtaue + ak
    tx = 1.0 / np.where(large, t1, 1.0)
    asym = np.exp(-t1) * tx * (1 + (-1 + (2 + (-6 + (24 + (-120 + 720*tx)*tx)*tx)*tx)*tx)*tx)
    result = np.where(small, 0.0, np.where(large, asym, soNtaue))
    return float(result[0]) if scalar else result

def expint_E2(t):
    t1 = np.abs(t)
    return np.exp(-t1) - t1 * expint_E1(t1)

def expint_E3(t):
    t1 = np.abs(t)
    return (np.exp(-t1) - t1 * expint_E2(t1)) / 2.0

def expint_E4(t):
    t1 = np.abs(t)
    return (np.exp(-t1) - t1 * expint_E3(t1)) / 3.0

def expint_E5(t):
    t1 = np.abs(t)
    return (np.exp(-t1) - t1 * expint_E4(t1)) / 4.0

def Bsun(nu):
    return (nu**3) / (jnp.exp(nu / Ts) - 1.0)

def Bearth(nu):
    return (nu**3) / (jnp.exp(nu / Te) - 1.0)

def BB(nu, T):
    return (nu**3) / (jnp.exp(jnp.fmin(nu / T, 100.0)) - 1.0)

def dBB(nu, T):
    # Vectorizable version: derivative of BB w.r.t. T
    safe_x = jnp.minimum(nu / T, 100.0)
    a = jnp.exp(safe_x)
    return a * (nu**2 / jnp.fmax(1.e-30, a - 1.0) / T)**2


# --- Initialisation des données ---
def readkappa(filepath):
    print(f"reading kappa(nu) from file {filepath}")
    if not os.path.exists(filepath):
        raise RuntimeError(f"Cannot open file: {filepath}")

    data = np.loadtxt(filepath)
    waveno = data[:, 0]
    kappaux = data[:, 1]
    jmax = len(waveno)
    kappanu_np = np.maximum(kappaux, kappamin)
    nu_np = 3.0 / waveno

    z_vals = np.arange(Nz) * dz

    # Vectorized: kappaz_np[j, i] = kapparhoz(nu_np[j], z_vals[i])
    kappaz_np = kapparhoz(nu_np[:, None], z_vals[None, :])  # (jmax, Nz)

    # Vectorized trapezoidal cumulative sum
    increments = (kappaz_np[:, 1:] + kappaz_np[:, :-1]) * dz / 2.0  # (jmax, Nz-1)
    kappasum_np = np.concatenate(
        [np.zeros((jmax, 1)), np.cumsum(increments, axis=1)], axis=1
    )  # (jmax, Nz)

    print(f"Number of frequencies {jmax}")
    return jmax, jnp.array(kappanu_np), jnp.array(nu_np), jnp.array(kappaz_np), jnp.array(kappasum_np)


# --- Cœur de l'algorithme ---
# Fully vectorized: no Python loops over frequencies or altitudes.
def updateJK_step(kappanu, nu, kappaz, kappasum, J0old, J2old, K0old, K2old, T):
    """
    kappanu: (jmax,)
    nu:      (jmax,)
    kappaz:  (jmax, Nz)
    kappasum:(jmax, Nz)
    J0old .. K2old: (jmax, Nz)
    T:       (Nz,)
    """
    z_arr = jnp.arange(Nz) * dz  # (Nz,)

    # --- Scattering and source terms: (jmax, Nz) ---
    kappazj = kappaz * kappanu[:, None]                          # (jmax, Nz)
    asz     = as_scatter(z_arr[None, :], nu[:, None])            # (jmax, Nz)
    sigs    = kappazj * asz
    siga    = kappazj * (1.0 - asz)

    BB_vals = BB(nu[:, None], T[None, :])                        # (jmax, Nz)
    S = sigs * J0old + siga * BB_vals                           # (jmax, Nz)
    H = 9.0 * beta * sigs / 8.0 * (J2old - J0old / 3.0 - K0old + K2old)

    # --- c0: (jmax,) ---
    ks_top = kappanu * kappasum[:, Nz - 1]                      # (jmax,)
    c0 = -q0 * mus * Cs * Bsun(nu) * jnp.exp(-ks_top / mus)    # (jmax,)
    E2_top = expint_E2(kappasum * kappanu[:, None])                # (jmax, Nz)
    c0 -= q0 * jnp.sum(E2_top * S * dz, axis=1)                 # (jmax,)

    # --- Boundary terms: (jmax, Nz) ---
    ksz = kappasum * kappanu[:, None]                            # (jmax, Nz)
    sun_attn = jnp.exp(-(kappasum[:, -1:] - kappasum) * kappanu[:, None] / mus)

    J0 = (Ce * Bearth(nu[:, None]) * expint_E3( ksz) / 2.0
          + Cs * Bsun(nu[:, None]) * sun_attn / 2.0
          + c0[:, None] * expint_E2( ksz) / 2.0)

    J2 = (Ce * Bearth(nu[:, None]) * expint_E5( ksz) / 2.0
          + Cs * (mus**2) * Bsun(nu[:, None]) * sun_attn / 2.0
          + c0[:, None] * expint_E4(ksz) / 2.0)

    # --- Convolution over j: (jmax, Nz, Nz-1) ---
    # C++ loop: for j=1..Nz-1, ip=max(i,j), jp=min(i,j)
    # arg = kappanu*(kappasum[ip]-kappasum[jp]-dz/2)
    Hj  = (H[:, 1:] + H[:, :-1]) / 2.0                          # (jmax, Nz-1)
    Sj  = (S[:, 1:] + S[:, :-1]) / 2.0 - Hj / 3.0               # (jmax, Nz-1)

    i_idx = jnp.arange(Nz)        # (Nz,)    observation points i=0..Nz-1
    j_idx = jnp.arange(1, Nz)     # (Nz-1,)  interval indices  j=1..Nz-1
    ip = jnp.maximum(i_idx[:, None], j_idx[None, :])  # (Nz, Nz-1)
    jp = jnp.minimum(i_idx[:, None], j_idx[None, :])  # (Nz, Nz-1)
    ks_ip = kappasum[:, ip]        # (jmax, Nz, Nz-1)
    ks_jp = kappasum[:, jp]        # (jmax, Nz, Nz-1)
    aux = kappanu[:, None, None] * (ks_ip - ks_jp - dz / 2.0)   # (jmax, Nz, Nz-1)
    absaux = jnp.abs(aux)

    E1 = expint_E1(absaux)                                        # (jmax, Nz, Nz-1)
    E3 = expint_E3(absaux)
    E5 = expint_E5(absaux)

    # Sum over j (axis=-1)
    J0 += jnp.sum((E1 * Sj[:, None, :] + E3 * Hj[:, None, :]) * (dz / 2.0), axis=-1)
    J2 += jnp.sum((E3 * Sj[:, None, :] + E5 * Hj[:, None, :]) * (dz / 2.0), axis=-1)
    K0  = -jnp.sum((E1 - E3) * Hj[:, None, :] * (dz / 2.0), axis=-1)
    K2  = -jnp.sum((E3 - E5) * Hj[:, None, :] * (dz / 2.0), axis=-1)

    return J0, J2, K0, K2


def genT(jmax, nu, kappanu, J0, T_old):
    T_new = np.array(T_old)
    z_pts = np.arange(Nz) * float(dz)

    # Precompute frequency-grid quantities shared across altitudes
    nu_np      = np.array(nu)
    kappanu_np = np.array(kappanu)
    J0_np      = np.array(J0)

    dnu         = nu_np[1:] - nu_np[:-1]                              # (jmax-1,)
    nu_mid      = (nu_np[1:] + nu_np[:-1]) / 2.0                      # (jmax-1,) midpoints
    kappa_mid   = (kappanu_np[1:] + kappanu_np[:-1]) / 2.0            # (jmax-1,)

    for i in range(Nz):
        zi = z_pts[i]

        # --- rhs (uses nu[j] not midpoints, kappanu[j] not averaged) ---
        asz_rhs = as_scatter(zi, nu_np[1:])                            # (jmax-1,)
        kappaux_rhs = (1.0 - asz_rhs) * kappanu_np[1:]
        rhs = float(np.sum(kappaux_rhs * (J0_np[1:, i] + J0_np[:-1, i]) / 2 * dnu))
        if rhs == 0:
            continue

        # --- kappaux for Newton: C++ uses kappanu[j] (not averaged) ---
        asz_mid = as_scatter(zi, nu_np[1:])                            # scattering at nu[j]
        kappaux_nwt = (1.0 - asz_mid) * kappanu_np[1:]                # (jmax-1,)

        # --- Bisection for initial bracket ---
        T0_val = float(T_new[i])
        T1_val = 2.0 * T0_val
        counter = 0

        def _root(Tv):
            BB_v = np.array(BB(jnp.array(nu_mid), Tv))
            return float(np.sum(kappaux_nwt * BB_v * dnu)) - rhs

        myeq0 = _root(T0_val)
        while myeq0 > 0 and counter < 100:
            T0_val /= 2.0
            counter += 1
            myeq0 = _root(T0_val)

        myeq1 = _root(T1_val)
        while myeq1 < 0 and counter < 100:
            T1_val *= 2.0
            myeq1 = _root(T1_val)
            counter += 1

        while (T1_val - T0_val) > epsdycho and counter < 100:
            Taux = (T1_val + T0_val) / 2.0
            counter += 1
            if _root(Taux) > 0:
                T1_val = Taux
            else:
                T0_val = Taux

        T_new[i] = (T1_val + T0_val) / 2.0

        # --- Newton refinement (vectorized over frequencies) ---
        presfunc = 1.0
        inewton = 0
        while inewton < newton and abs(presfunc) > epsnewton:
            inewton += 1
            Tv = T_new[i]
            nu_jax = jnp.array(nu_mid)
            BB_v  = np.array(BB(nu_jax, Tv))
            dBB_v = np.array(dBB(nu_jax, Tv))
            left  = float(np.sum(kappaux_nwt * BB_v  * dnu))
            deriv = float(np.sum(kappaux_nwt * dBB_v * dnu))
            presfunc = rhs - left
            if abs(deriv) > 1e-10:
                T_new[i] = Tv + presfunc / deriv

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
        J0, J2, K0, K2 = updateJK_step(kappanu, nu, kappaz, kappasum, J0old, J2old, K0old, K2old, T)
        T = genT(jmax, nu, kappanu, J0, T)

        J0old, J2old, K0old, K2old = J0, J2, K0, K2

        normG = float(jnp.sum(jnp.abs(J0) * jnp.diff(jnp.concatenate([nu[:1], nu])) [:, None] * dz))
        print(f"k= {k}  {(T[2] * 4798) - 273:.4f}  {normG:.6e}")

    return T


# --- Main ---
def main():
    try:
        jmax, kappanu_base, nu_base, kappaz, kappasum = readkappa(mykappafile)
    except Exception as e:
        print(f"Erreur lors de la lecture: {e}")
        return

    for K in range(0,2):
        kappanu = np.array(kappanu_base)

        for j in range(jmax):
            if K == 1:
                if (nu_base[j] > 3./18) and (nu_base[j] < 3./14):
                    kappanu[j] = max(0.5, 1.2 * float(kappanu_base[j]))
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
