# Vectorized version of earthsunnu.py. Much faster

import numpy as np
import os
import sys
import time

verbose =0 # Set to =1 for detailed output during execution

# Constants from original C++ code
Ce = 2.0 # Infrared intensity
Cs = 2e-6 # Solar intensity of collimated light
drho = -0.7 # gradient of density with altitude
Zatmo = 1.2 #Top of Atmosphere
Z = Zatmo * (1.0 + drho * Zatmo / 2.0) # computational TAO due to optical thickness
Bscale = 1.4744e-8 # scaling factor for I and J0...
Tscale = 4798.0 #scaling factor for the temperature
pi = np.pi
stefan = pi**4 / 15.0 # Stefan constant

Te = (273.0 + 18.0) / Tscale
Ts = 5798.0 / Tscale
q0 = -0.3  # Albedo intensity
mus = 0.5  # Angle of collimated light
z1 = Z * 0.5 # Cloud base
z2 = Z * 0.7 # Cloud top
z3 = Z * 0.8 # Rayleigh scttering base
nu1 = 0.5 # cloud is opaque between nu1 and nu2
nu2 = 1.0
lambda_val = 1.5  # Cloud opacity parameter 
cscat=0.3 # Rayleigh and enhanced scattering in clouds

# Grid parameters
Nz = 60 # Number of discretization points in the vertical direction
dz = Z / (Nz - 1)
kmax = 8 # Number of ISIF iterations
newton = 50 # Max number of Newton iterations
epsdycho = 1e-4  # Precision for dichotomy method
epsnewton = 1e-10 # Precision for Newton method
kappamin = 0.001 # Avoid zero opacity

# Initialize arrays
nu = []
kappanu = []
T = np.zeros(Nz)
z_mesh = np.linspace(0, Z, Nz)
zp_mesh = np.linspace(0, Z, Nz)
z_grid, zp_grid = np.meshgrid(z_mesh, z_mesh, indexing="ij")


# File paths
basedir = ""
mykappafile = os.path.join(basedir, "_kappa.txt")
myresulttemperature = os.path.join(basedir, "temperaturec")
myresultmeanintensity = os.path.join(basedir, "imean0")

def sqr(x):
    return x * x

# Scattering function
def ascat(z,nu):  # background scattering + cloud scattering
    aux= 0.3 + cscat*(z2-z)*(z-z1)*(z>z1)*(z<z2)*4/sqr(z1-z2) \
        + cscat*(nu<nu2)*(nu>nu1)*sqr(4*(nu-nu1)*(nu-nu2)/sqr(nu2-nu1))*(z>z3)  #+ Rayleigh
    return aux

# Cloud functions
def cloud(z):
    return 1.0 + lambda_val *((z>z1)^(z>z2))#( (z > z1)-(z > z2))

def scloud(z, zp):  # integral of cloud() over z,zp (zp>z)
    aux = zp - z + lambda_val * ((zp>z1)* (zp - z1) - (zp > z2) * (zp - z2) - (z > z1) * (z - z1) + (z > z2) * (z - z2))
    return np.abs(aux)

# Exponential integrals (optimized for speed)
def expint_E1(t):
    t1 = abs(t)+1e-10
    ak = t1 * 1.0 # to force ak to be an independent copy
    soNtaue =- 0.577215664901533 - np.log(t1) + ak
    for k in range(2, 10):
        ak *= -t1 * (k - 1) / (k * k)
        soNtaue += ak
    return soNtaue

def expint_E2(t):
    t1 = abs(t)
    return np.exp(-t1) - t1 * expint_E1(t1)

def expint_E3(t):
    t1 = abs(t)
    return (np.exp(-t1) - t1 * expint_E2(t1)) / 2.0

# Blackbody functions
def Bsun(nu):
    return nu**3 / (np.exp(nu / Ts) - 1.0)

def Bearth(nu):
    return nu**3 / (np.exp(np.fmin(nu / Te, 100.0)) - 1.0)

def BB(nu, T):
    return nu**3 / (np.exp(np.fmin(nu / T, 100.0)) - 1.0)

def dBB(nu, T):
    a =np.exp(np.fmin(nu / T, 100.0))
    return a * sqr(nu**2 / np.fmax(1e-30, a - 1.0) / T)

# Read kappa from file
def readkappa(mykappafile):
    print(f"Reading kappa(nu) from file {mykappafile}")
    with open(mykappafile, 'r') as kappafile:
        j = -1
        for line in kappafile:
            parts = line.split('\t')
            wavel, kappaux = map(float, line.strip().split())
            kappanu.append(max(kappaux, kappamin))
            nu.append(3.0/wavel)
            j += 1
    print(f"Number of frequencies {j}")
    dnu = np.diff(nu)
    auxe = 0.0
    auxs = 0.0
    for k in range(1, j):
        auxe += BB(nu[k], Te) * (nu[k] - nu[k-1])
        auxs += BB(nu[k], Ts) * (nu[k] - nu[k-1])
    
    print(f"Check Stefan id for the sun {auxs} {stefan * (Ts**4)}")
    print(f"Check Stefan id for the earth {auxe} {stefan * (Te**4)}")
    return j

# Update J0 (the mean angulat radiation intensity)
def updateJ():
    for jnu in range(jmax):
        kappanuj = kappanu[jnu]
        nuj = nu[jnu]
        bearthnu = Ce * Bearth(nuj)
        bsunnu = Cs * Bsun(nuj)
        c0 = -q0 * bsunnu * mus * np.exp(-scloud(0.0, Z) * kappanuj / mus)
        
        kappanuz = kappanuj * cloud(z_mesh)
        sigs = kappanuz * ascat(z_mesh, nuj)
        siga = kappanuz - sigs
        S = sigs * J0[jnu] + siga * BB(nuj, T)
        c0 -= q0 * np.dot(expint_E2(kappanuj * scloud(0.0, z_mesh)), S) * dz 
       
        J0z1 = bearthnu * expint_E3(kappanuj * scloud(0.0, z_mesh)) / 2.0
        J0z2 = bsunnu * np.exp(-kappanuj * scloud(z_mesh, Z) / mus) / 2.0
        J0z3 = c0 * expint_E2(kappanuj * scloud(0.0, z_mesh)) / 2.0
        J0z4 = np.dot(expint_E1(kappanuj * scloud(zp_grid, z_grid) + 0.5 * dz) ,S)*dz/2 
        J0[jnu] = J0z1+J0z2+J0z3 +J0z4

# Root function for dychotomy
def root(rhs, Tin, i):
    myeq = -rhs
    for j in range(1, jmax):
        myeq += (1.0 - ascat(i * dz, nu[j])) * kappanu[j] \
            * BB((nu[j] + nu[j-1]) / 2.0, Tin) * (nu[j] - nu[j-1])
    return myeq

# RHS for temperature equation
def rhsTeq(i):
    rhs = 0.0
    for j in range(1, jmax):
        rhs += (1.0 - ascat(i * dz, nu[j])) * kappanu[j] * J0[j][i] * (nu[j] - nu[j-1])
    return rhs

# Get temperature via dichotomy
def getTbydycho(rhs, i, Tstart):
    Taux, T0, T1 = Tstart, Tstart, 2.0 * Tstart
    if rhs == 0.0:
        return 0.0
    
    counter = 0
    myeq0 = root(rhs, T0, i)
    myeq1 = root(rhs, T1, i)
    
    while myeq0 > 0.0 and counter < 100:
        T0 /= 2.0
        counter += 1
        myeq0 = root(rhs, T0, i)
        if verbose:
            print(f"{T0} down {myeq0}")
    
    while myeq1 < 0.0 and counter < 100:
        T1 *= 2.0
        myeq1 = root(rhs, T1, i)
        counter += 1
        if verbose:
            print(f"{T1} up {myeq1}")
    
    if verbose:
        myeq0 = root(rhs, T0, i)
        myeq1 = root(rhs, T1, i)
        if myeq0 > myeq1:
            print("BUG in dychotomy")
    
    while abs(T1 - T0) > epsdycho and counter < 100:
        Taux = (T1 + T0) / 2.0
        counter += 1
        myeq0 = root(rhs, Taux, i)
        if myeq0 > 0.0:
            T1 = Taux
        else:
            T0 = Taux
        if verbose:
            myeq0 = root(rhs, T0, i)
            myeq1 = root(rhs, T1, i)
            print(f"{T0} middle {T1} {myeq0} {myeq1}")
    
    if counter > 99:
        print("Divergence in dichotomy")
    
    return (T1 + T0) / 2.0

# Generate temperature field
def genT():
    for i in range(Nz):
        rhs = rhsTeq(i)
        T[i] = getTbydycho(rhs, i, T[i])
        
        presfunc = 1.0
        inewton = 0
        while inewton < newton and abs(presfunc) > epsnewton:
            T0 = T[i]
            left = 0.0
            deriv = 0.0
            for j in range(1, jmax):
                dnu = nu[j] - nu[j-1]
                nu11 = (nu[j] + nu[j-1]) / 2.0
                kappaux = (1.0 - ascat(i * dz, nu[j])) * (kappanu[j] + kappanu[j-1]) / 2.0
                left += kappaux * BB(nu11, T0) * dnu
                deriv += kappaux * dBB(nu11, T0) * dnu
            
            presfunc = rhs - left
            if abs(deriv) > 1e-10:
                T[i] = T0 + presfunc / deriv
            
            if verbose:
                print(f"{i} T={T[i]} residue={presfunc} deriv={deriv} \
                      Tstefan={np.sqrt(np.sqrt(rhs * 15 * 2)) / 3.1416}")
        
        if inewton >= newton:
            print("Newton precision doubtful")

# Multi-block iteration
def multiBlock(initT):
    for i in range(Nz):
        T[i] = initT

    for k in range(kmax):
        updateJ()
        genT() 

        normG = 0.0
        for j in range(jmax):
            normG += np.dot(J0[j],J0[j])*dz
        print(f"{k}     {T[2]*4798-273}     {normG/Bscale}")

# Back to z function
def backtoz(z):
    return (np.sqrt(1.0 + 2.0 * drho * z) - 1.0) / drho

# Main function
def main():
    global jmax,J0, nu,kappanu
    jmax = readkappa(mykappafile)
    J0 =np.zeros((jmax,Nz))

    for K in range(1): # change to range(2) or range(2,3,1) to test different kappa
        for j in range(jmax):
            if K == 1:
                if nu[j] > 3.0 / 18.0 and nu[j] < 3.0 / 14.0:
                    kappanu[j] = 1.2 * kappanu[j]  # enhance kappa where CO2 absorbs nu
            elif K == 2:
                kappanu[j] = 0.5  # constant kappa
        
        resstream = myresulttemperature + "5" + str(K) + ".txt"
        if lambda_val == 0:
            resstream = myresulttemperature + str(K) + ".txt"
        
        with open(resstream, 'w') as resultfile:
            print(f"\n iterations   T[2]    Norm of J0")
            t0 = time.time()
            multiBlock(Te / 2.0)
            elapsed = (time.time() - t0) / 60.0
            print(f" Time CPU = {elapsed:.6f}")
            
            print(f"\n tau [T]:")
            for i in range(0, Nz, 1):
                z = i * dz
                print(f"{backtoz(z):.4f} {T[i] * Tscale - 273:.4f}")
            
            for i in range(Nz):
                resultfile.write(f"{backtoz(i * dz):.4f} {T[i] * Tscale - 273:.4f}\n")
        
        rb = Cs * stefan * (Ts**4)
        for j in range(jmax):
            rb -= J0[j][Nz-1] * (nu[j] - nu[j-1])
        print(f" Radiation Budget {rb}")

if __name__ == "__main__":
    main()
