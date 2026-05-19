import math
import numpy as np
import matplotlib.pyplot as plt

# Constantes globales 
Nz = 60
Nmu = 100
kmax = 15
tau1 = 1.0
dtau = tau1 / (Nz - 1)
mus = 1.0 / math.sqrt(3.0) # Une inclinaison de 60°
Is = 1.0 # Intensité de la lumière du soleil

# Initialisation des tableaux
J00 = [0.0] * Nz
J0 = [0.0] * Nz
S = [0.0] * Nz

def sqr(x):
    """Fonction utilitaire équivalente à la macro #define sqr(x)"""
    return x * x

def expint_E1(t=1.0):
    """Précision dans la fonction d'intégrale exponentielle E1"""
    t1 = abs(t) # permet les appels avec un argument négatif
    Kexpint = int(9 + (t1 - 1) * 4) 
    gaNtaua = 0.577215664901533 # intégration spéciale pour log(t)
    
    if t1 < 1e-5:
        return 0.0 # E1(t), t petit est toujours appelé avec t*E1(t)
    
    if t1 > 4.0:
        tx = 1.0 / t1
        return math.exp(-t1) * tx * (1 + (-1 + (2 + (-6 + (24 + (-120 + 720 * tx) * tx) * tx) * tx) * tx) * tx)
    
    ak = t1
    soNtaue = -gaNtaua - math.log(t1) + ak
    for k in range(2, Kexpint):
        ak *= -t1 * (k - 1) / sqr(k)
        soNtaue += ak
        
    return soNtaue

def getI():
    pi = 4.0 * math.atan(1.0)
    dmu = 2.0 / Nmu
    b = 1 # 1 pour Chandrasekhar, 0 pour Dirac
    lam = 12.0
    w = math.sqrt(pi) / 2.0 / lam * (math.erf(mus * lam) + math.erf((1.0 - mus) * lam))
    
    print(" tau,mu --> I(tau,mu) ")
    
    mu = -1.0
    while mu < 1.0:
        for j in range(Nz):
            # En Python, int(True) vaut 1 et int(False) vaut 0, ce qui remplace la syntaxe (mu<0) du C++
            Iss = int(mu < 0) * Is * math.exp(-sqr((mu + mus) * lam)) * math.exp(-j * dtau / mus) / w
            I = (1 - b) * Iss
            
            for i in range(Nz):
                condition = int(i < j) * int(mu < 0) + int(i > j) * int(mu > 0)
                if abs(mu) > 0: # Évite la division par zéro
                    I += condition * dtau * math.exp(-abs((i - j) * dtau / abs(mu))) \
                         * (J0[i] + b * Is * math.exp(-i * dtau / mus) / 2.0) / abs(mu)
                         
            print(f"{j * dtau:.5f}\t{mu:.5f}\t{I + b * Iss:.5f}")
        print()
        mu += dmu
        
    return 0
def getI2():
    pi = 4.0 * math.atan(1.0)
    dmu = 2.0 / Nmu
    b = 1  # 1 pour Chandrasekhar, 0 pour Dirac
    lam = 12.0
    w = math.sqrt(pi) / 2.0 / lam * (math.erf(mus * lam) + math.erf((1.0 - mus) * lam))
    
    # Listes pour stocker les données du graphique
    mu_vals = []
    tau_vals = [j * dtau for j in range(Nz)]
    I_surface = []
    
    print("Calcul de la surface I(tau, mu) en cours...")
    
    mu = -1.0
    while mu <= 1.0:
        mu_vals.append(mu)
        I_row = []
        
        for j in range(Nz):
            Iss = int(mu < 0) * Is * math.exp(-sqr((mu + mus) * lam)) * math.exp(-j * dtau / mus) / w
            I = (1 - b) * Iss
            
            for i in range(Nz):
                condition = int(i < j) * int(mu < 0) + int(i > j) * int(mu > 0)
                if abs(mu) > 1e-9:  # Sécurité pour éviter la division par zéro
                    I += condition * dtau * math.exp(-abs((i - j) * dtau / abs(mu))) \
                         * (J0[i] + b * Is * math.exp(-i * dtau / mus) / 2.0) / abs(mu)
                         
            I_row.append(I + b * Iss)
            
        I_surface.append(I_row)
        mu += dmu

    # Conversion en tableaux NumPy pour la reconstruction de la surface
    Mu, Tau = np.meshgrid(mu_vals, tau_vals, indexing='ij')
    I_surface = np.array(I_surface)
    
    # --- Création du graphique 3D ---
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111, projection='3d')
    
    # Tracé de la surface
    surf = ax.plot_surface(Tau, Mu, I_surface, cmap='viridis', edgecolor='none', alpha=0.9)
    
    # Personnalisation des axes
    ax.set_xlabel(r'$\tau$ (Profondeur optique)', fontsize=12)
    ax.set_ylabel(r'$\mu$ ($\cos(\theta)$)', fontsize=12)
    ax.set_zlabel(r'$I(\tau, \mu)$', fontsize=12)
    ax.set_title("Surface de l'intensité spécifique $I(\tau, \mu)$ (Problème de Milne)", fontsize=14, pad=20)
    
    # Ajout d'une barre de couleur
    fig.colorbar(surf, ax=ax, shrink=0.5, aspect=10, label="Intensité $I$")
    
    # Orientation initiale de la vue 3D (ajustable à la souris)
    ax.view_init(elev=30, azim=-60)
    
    plt.tight_layout()
    plt.show()
    return 0

def schwarzschild(tau):
    J00_val = Is / 2.0
    B = J00_val / tau1 / 2.0
    A = J00_val - B
    exp2tau = math.exp(-2.0 * tau)
    return 2.0 * (J00_val * exp2tau + A * (1.0 - exp2tau) - B * (2.0 * tau - 1.0 + exp2tau))

def main():
    global J0, S
    
    # Initialisation
    for i in range(Nz):
        J0[i] = math.exp(-i * dtau / mus) * Is / 2.0
        
    print("k     J0[1]")
    
    # Itérations sur la source (point fixe)
    for k in range(kmax):
        for i in range(Nz):
            S[i] = J0[i]
            
        for i in range(Nz):
            J0[i] = math.exp(-i * dtau / mus) * Is / 2.0
            for j in range(1, Nz):
                J0[i] += dtau / 4.0 * expint_E1((j - i - 0.5) * dtau) * (S[j] + S[j - 1])
                
        print(f"{k} {J0[1]}")
        
    print("\ntau      J0[tau]      Schwarzschild")
    for i in range(Nz):
        print(f"{i * dtau:.5f} {J0[i]:.5f} {schwarzschild(i * dtau):.5f}")
        
getI2() # Uncomment to display  I(tau,mu)

if __name__ == "__main__":
    main()