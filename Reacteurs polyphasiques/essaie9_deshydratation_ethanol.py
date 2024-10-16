import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt  # Pour les tracés
import pandas as pd

# Toutes les données sont en unités SI
"""Objectif: Vcata mini, (Lmax, N, Uv0 Dr) à rajouter kp apparent, """
# Paramètres invariants
Rpore = 10e-9      # Rayon de pore en m (fixé à 10 nm)
Rhop = 1100        # Masse volumique du catalyseur en kg/m³
Mug = 17.5e-6      # Viscosité dynamique du gaz supposée constante en Pa·s
R = 8.314          # Constante des gaz parfaits en J/(mol·K)
P0 = 2 * 10e5      # Pression initiale en Pa en entrée de tube
T0 = 900           # Température initiale en K
Metoh = 0.046      # Masse molaire de l'éthanol en kg/mol
Malcene = 0.028    # Masse molaire de l'éthène en kg/mol
Meau = 0.018       # Masse molaire de l'eau en kg/mol
Eps = 0.4          # Porosité du lit
Ea = 150 * 10e3    # Énergie d'activation en J/mol
DrH = 49.1e3       # Enthalpie de réaction à 900 K en J/mol
Tp = 900           # Température de paroi considérée fixe, en K 
Cpetoh = 2953      # Capacité calorifique massique de l'éthanol en J/(kg·K)
Cpalcene = 3186    # Capacité calorifique massique de l'éthène en J/(kg·K)
Cpeau = 2217       # Capacité calorifique massique de l'eau en J/(kg·K)
Lambdag = 0.025    # Conductivité thermique en W/(m·K)
K0 = 0.66          # Constante cinétique à 400°C (ou 673.15 K) en s⁻¹

# Paramètres fixes ajustés
Uv0 = 0.4    # Vitesse superficielle en m/s 
Dr = 0.15     # Diamètre de tube en m 
Qm0 = 27.77  # Débit massique initial d'éthanol en kg/s
Dp = 0.003         # Diamètre des particules en m (ajusté à 3 mm)

# Paramètres correction cinétique
Dm0 = 1e-5  # Diffusivité moléculaire de référence en m²/s
Kd = 0.001  # Coefficient de transfert de masse externe

# Fractions massiques initiales
Wetoh0 = 1.0        # Fraction massique initiale d'éthanol 
Walcene0 = 0.0      # Fraction massique initiale d'éthène 
Weau0 = 0.0         # Fraction massique initiale d'eau 

# Calcul de la masse volumique initiale du gaz
M0 = Metoh * Wetoh0 + Malcene * Walcene0 + Meau * Weau0  # Masse molaire moyenne de mélange de gaz
Rho0 = P0 * M0 / (R * T0)  # Masse volumique initiale du gaz en kg/m³

# Débit volumique initial
Qv0 = Qm0 / Rho0   # m³/s

# Définition de la fonction du solveur
def odestiff(z, E):
    # Extraire les 5 variables d'état dans un vecteur ligne appelé E
    Wetoh, Walcene, Weau, T, P = E  # vecteur E

    # Calculer fractions molaires yi
    numerator = np.array([Wetoh / Metoh, Walcene / Malcene, Weau / Meau])
    denominator = np.sum(numerator)
    y_etoh = numerator[0] / denominator

    # Calculer concentration molaire totale
    C_tot = P / (R * T)  # mol/m³
    # Calculer concentration molaire EtOH
    C_etoh = y_etoh * C_tot  # mol/m³

    # Calculer Masse molaire moyenne
    Mm_inv = (Wetoh / Metoh) + (Walcene / Malcene) + (Weau / Meau)
    Mm = 1 / Mm_inv  # kg/mol

    # Calculer constante de vitesse (cinétique Arrhenius)
    T_ref = 673.15  # Température de référence pour notre K0 en Kelvin
    k = K0 * np.exp(-Ea / R * ((1 / T) - (1 / T_ref)))

    # Correction cinétique pour diffusion et transfert de matière externe
    L = Dp / 6  # Longueur caractéristique
    Dm = Dm0 * (T / 273.15)**(3/2)  # Diffusivité moléculaire moyenne
    Dk = (4/3) * Rpore * np.sqrt(2 * R * T / (np.pi * Mm))  # Diffusivité de Knudsen
    De = 0.5 * (1 / Dm + 1 / Dk)**-1  # Diffusivité effective
    Th = np.sqrt(k / De) * L  # Module de Thiele ordinaire
    Nusu = np.tanh(Th) / Th  # Efficacité de surface pour une réaction du 1er ordre
    Fex = 1 / (1 + Kd / (Nusu * L * k))  # Fraction de résistance externe
    Ka = k * Nusu * (1 - Fex)  # Constante de vitesse apparente

    # Vitesse de réaction avec la constante de vitesse apparente
    r = Ka * C_etoh  # mol/(m³·s)

    # Calculs intermédiaires
    Cpm = Wetoh * Cpetoh + Walcene * Cpalcene + Weau * Cpeau  # Capacité calorifique moyenne

    # Densité du mélange
    Rho = P * Mm / (R * T)  # kg/m³

    # Reynolds et Prandtl numbers
    Re_p = (Rho * Uv0 * Dp) / Mug
    Pr = Mug * Cpm / Lambdag

    # Nusselt number et coefficient de transfert de chaleur
    Nu = 2 + 1.8 * np.sqrt(Re_p) * Pr**(1/3)
    
    
    h = 150  # Coefficient de transfert fixé

    # Bilan de matière
    coeff_massique = (1 - Eps) / (Eps * Rho * Uv0)

    # Initialiser vecteur dérivé
    res = np.zeros(5)

    # Bilan de matière
    res[0] = -coeff_massique * r * Metoh    # Éthanol
    res[1] = coeff_massique * r * Malcene   # Éthène
    res[2] = coeff_massique * r * Meau      # Eau

    # Bilan d'énergie
    Q_reaction = -(1 - Eps) * r * DrH / Mm  # Réaction endothermique
    Q_heat_transfer = (4 * h / Dr) * (Tp - T)  # Bilan d'échange thermique avec la paroi
    res[3] = (Q_heat_transfer + Q_reaction) / (Rho * Uv0 * Cpm)

    # Bilan quantité de mouvement
    friction_factor = 1.75 + 150 * (1 - Eps) * Mug / (Rho * Uv0 * Dp)
    res[4]
    res[4] = - (1 - Eps) * friction_factor * Rho * Uv0**2 / (Eps**3 * Dp)

    return res

# Conditions initiales
E0 = [Wetoh0, Walcene0, Weau0, T0, P0]

# Définition de la plage de z pour l'intégration
Lmax = 30  # Longueur maximale du réacteur en mètres 
z_span = np.linspace(0, Lmax, 10000)  # Points de z pour l'intégration

# Appel au solveur solve_ivp avec la méthode BDF (adaptée aux problèmes rigides)
sol = solve_ivp(
    fun=odestiff, 
    t_span=(0, Lmax), 
    y0=E0, 
    method='BDF',  # Méthode pour les problèmes raides
    t_eval=z_span,  # Points où évaluer la solution
    rtol=1e-6, 
    atol=1e-8
)

# Extraction des résultats
Wetoh_sol = sol.y[0]
Walcene_sol = sol.y[1]
Weau_sol = sol.y[2]
T_sol = sol.y[3]
P_sol = sol.y[4]

# Recalcul de Mm_inv et Mm à chaque position z
Mm_inv_sol = (Wetoh_sol / Metoh) + (Walcene_sol / Malcene) + (Weau_sol / Meau)
Mm_sol = 1 / Mm_inv_sol  # kg/mol

# Recalcul de la densité à chaque position z
Rho_sol = (P_sol * Mm_sol) / (R * T_sol)  # kg/m³

# Calcul des concentrations massiques
C_mass_etoh = Rho_sol * Wetoh_sol  # kg/m³
C_mass_alcene = Rho_sol * Walcene_sol  # kg/m³
C_mass_eau = Rho_sol * Weau_sol  # kg/m³

# Affichage des résultats dans un DataFrame
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
plt.plot(z_span, Wetoh_sol, label='Éthanol')
plt.plot(z_span, Walcene_sol, label='Éthène')
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
