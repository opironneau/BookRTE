#ISIF applied to the RTE in the grey case
 # $$
 #  \mu\partial_z I + \kappa I = \kappa(1-a) \sigma T^4 + a\kappa\int_{-1}^1I d\mu
 #    $$


import numpy as np
import time

# --- Constants & Configuration ---
Ce, Cs = 2.0, 0.2e-5    # Intensities of EM radiation from Earth and Sun
drho = -0.7  # density gradient
Zatmo = 1.2 # TAO (top of atmosphere in km)
Z = Zatmo * (1 + drho * Zatmo / 2)  # optical thickness

Te = (273 + 18) / T0 # Earth surface temperature in scaled units
Ts = 5800 / T0 # Sun apparent temperature in scaled units
kappanu = 0.5 # absorption coefficient (divided by density)
q0 = -0.3   # intensity of albeedo
mus = 0.5  # cosine of direction of collimated solar radiation

Nz = 300  # nb points in z
quadra = 1  # 0,1 quadrature precision order
kmax = 15  # nb fixed point iterations
dz = Z / (Nz - 1)

B0, T0 = 1.4744e-8, 4798.0 # scaling constants 
pi, stefan = np.pi,(pi**4) / 15

# --- Helper Functions ---

def expint_E1(t):
    """ Précision controled by Kexpint """
    t1 = abs(t)+1e-15 # strictly >0 argumeent only
    Kexpint = 10
    gaNtaua = 0.577215664901533  
    ak = t1
    soNtaue = -gaNtaua - np.log(t1) + ak
    for k in range(2, Kexpint):
        ak *= -t1 * (k - 1) / (k*k)
        soNtaue += ak
        
    return soNtaue

def expint_E2(t):
    t1=abs(t);
    return np.exp(-t1) - t1*expint_E1(t1);

def expint_E3(t):
    t1=abs(t);
    return (np.exp(-t1) - t1*expint_E2(t1))/2;

def expint_E1b(t):
    t1 = abs(t)+1e-15 
    Kexpint = 10
    gaNtaua = 0.577215664901533
    ak = t1
    soNtaue = -gaNtaua + ak # log part handled by hand
    for k in range(2, Kexpint):
        ak *= -t1 * (k - 1) / (k*k)
        soNtaue += ak
        
    return soNtaue


# --- Initialization ---
z_mesh = np.linspace(0, Z, Nz)
T = np.zeros(Nz)    # Temperature array
S = np.zeros(Nz)    # Energy source array
J0 = np.zeros(Nz)   # Radiative flux array
J0log = np.zeros(Nz) # Logarithmic part of the quadrature for E1 if quadra=1
I0 = Ce*stefan*Te**4
IZ = Cs*stefan*Ts**4
# Precompute J00 and J0Z arrays vectorially
J00Z = 0.5 * I0 * expint_E3(kappanu * z_mesh) + 0.5 * IZ * np.exp(-kappanu * (Z - z_mesh)/mus) 
# Vectorized  for integrations
# Creates 2D grids of z (rows) and zp (columns) to eliminate internal loops
z_grid, zp_grid = np.meshgrid(z_mesh, z_mesh, indexing="ij")


# --- Main Solver ---
def solve():
    global T, S, J0
    # Initialize T and S
    T[:] = Te/2
    S[:] = stefan * T**4
    J0 = J00Z
    
    start_time = time.time()
        # ISIF loop
    for k in range(kmax):
        Q0 = -q0*kappanu*mus*IZ*np.exp(-kappanu*Z/mus) 
        Q0 -= q0*kappanu*dz*np.dot(expint_E2(kappanu * z_mesh),J0)
        
        diff_z = zp_grid - z_grid
        J0z = J00Z + 0.5*Q0*expint_E2(kappanu*z_mesh)
        # quadrature control
        if quadra == 0:
            J0 = J0z  +  0.5*dz*np.dot(expint_E1(kappanu *(-dz/2+diff_z)), S)

        elif quadra == 1:
            # the log part of E_1 is treated separately 
            J0log = (diff_z) * (np.log(kappanu * ( 1e-10+ np.abs(diff_z))) - 1) - (diff_z - dz) * (np.log(kappanu * ( 1e-10+ np.abs(diff_z-dz))) - 1)
            J0 = J0z +  0.5*dz*np.dot(expint_E1b(kappanu *(-dz/2+diff_z)), S) - 0.5 * np.dot(J0log, S)

        # Update Temperature via Fixed Point
        T = (J0/stefan)**0.25

        #  updates to S
        S[0]=0
        for i in range(1,Nz):
            S[i] = kappanu * (J0[i]+J0[i-1])/2

        print(f"T(0)= {T[0]*T0-273:.4f} at iteration {k}")

    print(f"Calculation finished in {time.time() - start_time:.4f} seconds.")
    
    print("Altitude  T  J0")
    for i in range(0,Nz,Nz-1):
        print(f"{(np.sqrt(1+2*drho*z_mesh[i])-1)/drho :.4f}   {T[i]*T0-273  :.4f}   {J0[i]/B0:.4f}")


if __name__ == "__main__":
    solve()