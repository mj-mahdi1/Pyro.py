"""AVEC Solveur Odeint au lieu de solve_ivp"""
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt  # Pour les tracés
import pandas as pd

# Toutes les données sont en unités SI
"""Objectif: Vcata mini, (Lmax, N, Uv0 Dr) à rajouter kp apparent, """
# Paramètres invariants
Rpore = 10e-9      # Rayon de pore en m (fixé à 10 nm)
Rhop = 1100        # Masse volumique du catalyseur en kg/m³
Mug = 17.5e-6      # Viscosité dynamique du gaz supposée constante en Pa·s
R = 8.314          # Constante des gaz parfaits en J/(mol·K)
P0 = 2*10e5        # Pression initiale en Pa en entrée de tube
T0 = 900           # Température initiale en K
Dp = 0.003         # Diamètre des particules en m (ajusté à 3 mm)
Metoh = 0.046      # Masse molaire de l'éthanol en kg/mol
Malcene = 0.028    # Masse molaire de l'éthène en kg/mol
Meau = 0.018       # Masse molaire de l'eau en kg/mol
Eps = 0.4          # Porosité du lit
Ea = 150*10e3      # Énergie d'activation en J/mol
DrH = 49.1e3       # Enthalpie de réaction à 900 K en J/mol
Tp = 900           # Température de paroi considérée fixe, en K 
""""A varier Tp"""
Cpetoh = 2953      # Capacité calorifique massique de l'éthanol en J/(kg·K)
Cpalcene = 3186    # Capacité calorifique massique de l'éthène en J/(kg·K)
Cpeau = 2217       # Capacité calorifique massique de l'eau en J/(kg·K)
Lambdag = 0.025    # Conductivité thermique en W/(m·K)
K0 = 0.66          # Constante cinétique à 400°C (ou 673.15 K) en s⁻¹
"""A définir A,B,C,D"""
# Paramètres fixes ajustés
Uv0 = 0.5    # Vitesse superficielle en m/s 
"""si Uv0 est moins que 0.45, il n y a plus de conservation de matière dans ces conditions"""
Dr = 0.1     # Diamètre de tube en m 
"""pour l'instant, aucun effet de changement de diamètre de tube sur la simulation"""
Qm0 = 27.77     # Débit massique initial d'éthanol en kg/s #Fixé

#Paramètres correction cinétique
Dm0 = 1e-5  # Diffusivité moléculaire de référence en m²/s (à ajuster selon le système)
Kd = 0.001  # Coefficient de transfert de masse externe (à ajuster selon besoin)

# Fractions massiques initiales
Wetoh0 = 1.0        # Fraction massique initiale d'éthanol 
Walcene0 = 0.0      # Fraction massique initiale d'éthène 
Weau0 = 0.0         # Fraction massique initiale d'eau 

# Calcul de la masse volumique initiale du gaz
M0 = Metoh * Wetoh0 + Malcene * Walcene0 + Meau * Weau0 # Masse molaire moyenne de mélange de gaz
Rho0 = P0 * M0 / (R * T0)  # Masse volumique initiale du gaz en kg/m³

# Débit volumique initial
Qv0 = Qm0 / Rho0   # m³/s

# Définition de la fonction du solveur
def odestiff(E, z):
    # Extraire les 5 variables d'état dans un vecteur ligne appelé E
    Wetoh, Walcene, Weau, T, P = E #vecteur E

    # Calculer fractions molaires yi
    numerator = np.array([Wetoh / Metoh, Walcene / Malcene, Weau / Meau])
    denominator = np.sum(numerator)
    y_etoh = numerator[0] / denominator
    """fractions molaires y_alcene et y_eau sont toujours égales à l'inverse pour chaque terme"""

    # Calculer concentration molaire totale
    C_tot = P / (R * T)  # mol/m³
    # Calculer concentration molaire EtOH
    C_etoh = y_etoh * C_tot  # mol/m³

    # Calculer Masse molaire moyenne priliminaire
    Mm_inv = (Wetoh / Metoh) + (Walcene / Malcene) + (Weau / Meau)
    Mm = 1 / Mm_inv  # kg/mol

    # Calculer constante de vitesse
    T_ref = 673.15  # Température de référence pour notre K0 en Kelvin
    
    k = K0 * np.exp(-Ea / R * ((1 / T) - (1 / T_ref))) 
    #k = min(k, 0.1)  # Arhenus biponctuel
    #Correction cinetique pour tenir en compte de la diffusion et transfert de matière externe et ses limitations qui nous donne Kapp à la fin

    # Calculer la longueur caractéristique pour le module de Thiele
    L = Dp / 6  # Longeur caractéristique
    # Calculer la diffusivité moléculaire moyenne
    Dm = Dm0 * (T / 273.15)**(3/2)  #depend de la T, Dm0 étant la valeur de base avant variation de T
    # Calculer la diffusivité de Knudsen
    Dk = (4/3) * Rpore * np.sqrt(2 * R * T / (np.pi * Mm))
    # Diffusivité effective
    De = 0.5 * (1 / Dm + 1 / Dk)**-1  #elle combine les effets de Dm et Dk
    # Module de Thiele ordinaire
    Th = np.sqrt(k / De) * L #mesure l'importance relative de cinétique intrinsèque vs diffusion
    # Efficacité de surface pour une réaction du premier ordre
    Nusu = np.tanh(Th) / Th     #représente les effets de limitation par diffusion, réduit l'efficacité de vitesse de rxn
    # Fraction de résistance externe
    Fex = 1 / (1 + Kd / (Nusu * L * k)) #représente la résistance externe de transfert de matière
    # Constante de vitesse apparente
    Ka = k * Nusu * (1 - Fex) 
    # Vitesse de réaction avec la constante de vitesse apparente
    r = Ka * C_etoh  # mol/(m³·s)



    # vitesse de réaction
    r = Ka * C_etoh  # mol/(m³·s)

    # Calculs intermédiaires
    Cpm = Wetoh * Cpetoh + Walcene * Cpalcene + Weau * Cpeau 
    """Cp moyen de phase gazeuse il faut utiliser Cpi(T) = A + BT + CT^2 + DT^3 pour l'évolution de Cpi"""
    
    # Masse molaire moyenne
    Mm_inv = Wetoh / Metoh + Walcene / Malcene + Weau / Meau 
    Mm = 1 / Mm_inv  # kg/mol
    
    # Densité du mélange
    Rho = P * Mm / (R * T)  # kg/m³

    # Reynolds et Prandtl numbers
    Re_p = (Rho * Uv0 * Dp) / Mug
    Pr = Mug * Cpm / Lambdag
    
    # Nusselt number et coefficient de transfert de chaleur
    Nu = 2 + 1.8 * np.sqrt(Re_p) * Pr**(1/3)
    h = (Lambdag / Dp) * Nu  # Coefficient de transfert de chaleur interne 
    #hint = 100 # W/m/K varie de [100-500]
    # coefficient de transfert de matière
    coeff_massique = (1 - Eps) / (Eps * Rho * Uv0) 

    # Initialiser vecteur dérivé
    res = np.zeros(5) 

    # Bilan de matière
    res[0] = -coeff_massique * r * Metoh    # Éthanol
    res[1] = coeff_massique * r * Malcene   # Éthylène
    res[2] = coeff_massique * r * Meau      # Eau

    # Bilan d'énergie
    Q_reaction =  -(1 - Eps) * r * DrH / Mm  # r est la vitesse de réaction explicité précédemment
    Q_heat_transfer = (4 * h / Dr) * (Tp - T)
    res[3] = (Q_heat_transfer + Q_reaction) / (Rho * Uv0 * Cpm)

    # Bilan quantité de mouvement
    friction_factor = 1.75 + 150 * (1 - Eps) * Mug / (Rho * Uv0 * Dp)
    res[4] = - (1 - Eps) * friction_factor * Rho * Uv0**2 / (Eps**3 * Dp)

    

    return res


# Conditions initiales
E0 = [Wetoh0, Walcene0, Weau0, T0, P0]

# Définition de la plage de z pour l'intégration
Lmax = 10  # Longueur maximale du réacteur en mètres 
z_span = np.linspace(0, Lmax, 10000)  # Points de z pour l'intégration génère des points de façon linéairement espacée entre deux valeurs.
#adapter la répartition et la précision des points générés en fonction des besoins de simulation
# Appel au solveur odeint format: [np.linspace(start, stop, num)]

sol = odeint(odestiff, E0, z_span)

# Extraction des résultats
Wetoh_sol = sol[:, 0]
Walcene_sol = sol[:, 1]
Weau_sol = sol[:, 2]
T_sol = sol[:, 3]
P_sol = sol[:, 4]

# Recalcul de Mm_inv et Mm à chaque position z
Mm_inv_sol = (Wetoh_sol / Metoh) + (Walcene_sol / Malcene) + (Weau_sol / Meau)
Mm_sol = 1 / Mm_inv_sol  # kg/mol 

# Recalcul de la densité à chaque position z
Rho_sol = (P_sol * Mm_sol) / (R * T_sol)  # kg/m³

# Calcul des concentrations massiques
C_mass_etoh = Rho_sol * Wetoh_sol  # kg/m³
C_mass_alcene = Rho_sol * Walcene_sol  # kg/m³
C_mass_eau = Rho_sol * Weau_sol  # kg/m³

# Affichage des résultats dans Terminal
df = pd.DataFrame({
    "z": z_span,  
    "Wetoh": Wetoh_sol,
    "Walcene": Walcene_sol,
    "Weau": Weau_sol,
    "T": T_sol,
    "P": P_sol,
    "Masse Molaire": Mm_sol
})

print(df)

# Tracés des concentrations massiques en fonction de z
plt.figure(figsize=(11, 8))
plt.plot(z_span, C_mass_etoh, label='Concentration massique d\'éthanol')
plt.plot(z_span, C_mass_alcene, label='Concentration massique d\'éthène')
plt.plot(z_span, C_mass_eau, label='Concentration massique d\'eau')
plt.xlabel('Position z (m)')
plt.ylabel('Concentration massique (kg/m³)')
plt.title('Profil des concentrations massiques le long du réacteur')
plt.xlim(0, Lmax)  
plt.legend()
plt.grid(True)
plt.show()

# Tracé de la température en fonction de z
plt.figure(figsize=(11, 8))
plt.plot(z_span, T_sol, label='Température du fluide')
plt.xlabel('Position z (m)')
plt.ylabel('Température (K)')
plt.title('Profil de température le long du réacteur')
plt.xlim(0, Lmax) 
plt.ylim(500, 1500)  
plt.legend()
plt.grid(True)
plt.show()

# Tracer fractions massiques au long du tube
plt.figure(figsize=(11, 8))
plt.plot(z_span, Wetoh_sol, label='Ethanol')
plt.plot(z_span, Walcene_sol, label='Ethène')
plt.plot(z_span, Weau_sol, label='Eau')
plt.xlabel('Position z (m)')
plt.ylabel('Fraction massique')
plt.title('Fractions massiques au long du réacteur')
plt.legend()
plt.grid(True)
plt.show()

# Tracé de la pression en fonction de z
plt.figure(figsize=(11, 8))
plt.plot(z_span, P_sol, label='Pression le long du réacteur')
plt.xlabel('Position z (m)')
plt.ylabel('Pression (Pa)')
plt.title('Profil de pression le long du réacteur')
plt.xlim(0, Lmax)  
plt.legend()
plt.grid(True)
plt.show()

# Vérification de la conservation massique
sum_Wi = Wetoh_sol + Walcene_sol + Weau_sol
plt.figure(figsize=(11, 8))
plt.plot(z_span, sum_Wi, label='Somme de fractions massiques')
plt.xlabel('Position z (m)')
plt.ylabel('Somme des fractions massiques')
plt.title('Vérification de conservation de matière')
plt.legend()
plt.grid(True)
plt.show()
