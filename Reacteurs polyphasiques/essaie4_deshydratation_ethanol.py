import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt  # Pour les tracés

# Toutes les données sont en unités SI

# Paramètres invariants
Rpore = 10e-9      # Rayon de pore en m (fixé à 10 nm)
Rhop = 1100        # Masse volumique du catalyseur en kg/m³
Mug = 17.5e-6      # Viscosité dynamique du gaz supposée constante en Pa·s
R = 8.314          # Constante des gaz parfaits en J/(mol·K)
P0 = 1e5           # Pression initiale en Pa (ajustée à 1 bar)
T0 = 900           # Température initiale en K
Dp = 0.003         # Diamètre des particules en m (ajusté à 3 mm)
Metoh = 0.046      # Masse molaire de l'éthanol en kg/mol
Malcene = 0.028    # Masse molaire de l'éthène en kg/mol
Meau = 0.018       # Masse molaire de l'eau en kg/mol
Eps = 0.4          # Porosité du lit
Ea = 150e3         # Énergie d'activation en J/mol
DrH = 49.1e3       # Enthalpie de réaction à 900 K en J/mol
Tp = 900           # Température de paroi considérée fixe, en K
Cpetoh = 2953      # Capacité calorifique massique de l'éthanol en J/(kg·K)
Cpalcene = 3186    # Capacité calorifique massique de l'éthène en J/(kg·K)
Cpeau = 2217       # Capacité calorifique massique de l'eau en J/(kg·K)
Lambdag = 0.025    # Conductivité thermique en W/(m·K)
K0 = 0.66          # Constante cinétique à 400°C (ou 673.15 K) en s⁻¹

# Paramètres fixes ajustés
Uv0 = 0.2    # Vitesse superficielle en m/s 
Dr = 0.1     # Diamètre de tube en m 
Qm0 = 27.77     # Débit massique initial d'éthanol en kg/s 

# Fractions massiques initiales
Wetoh0 = 1.0        # Fraction massique initiale d'éthanol 
Walcene0 = 0.0      # Fraction massique initiale d'éthène 
Weau0 = 0.0         # Fraction massique initiale d'eau 

# Calcul de la masse volumique initiale du gaz
M0 = Metoh * Wetoh0 + Malcene * Walcene0 + Meau * Weau0
Rho0 = P0 * M0 / (R * T0)  # Masse volumique initiale du gaz en kg/m³

# Débit volumique initial
Qv0 = Qm0 / Rho0   # m³/s

# Définition de la fonction solveur
def solveur(z, E): #coordonnées spaciales z et vecteur E de parametres
    # Extraction des variables d'état (definition de vecteur E)
    Wetoh, Walcene, Weau, T, P = E

    # Calcul des fractions molaires
   
    numerator = np.array([Wetoh / Metoh, Walcene / Malcene, Weau / Meau])
    denominator = np.sum(numerator)
    y_etoh = numerator[0] / denominator  # Fraction molaire de l'éthanol

    # Calcul de la concentration totale
    C_tot = P / (R * T)  # mol/m³ avec P étant la pression totale
    # Concentration molaire de l'éthanol
    C_etoh = y_etoh * C_tot  # mol/m³ concentration molaire EtOH

    # Calcul de k avec la bonne température de référence
    T_ref = 673.15  # Température de référence pour K0 en K (400°C)
    k = K0 * np.exp(-Ea / R * ((1 / T) - (1 / T_ref)))  # Constante cinétique à la température T

    # Taux de réaction
    r = k * C_etoh  # mol/(m³·s)

    # Calculs intermédiaires
    Cpm = Wetoh * Cpetoh + Walcene * Cpalcene + Weau * Cpeau #W est la masse 
    # Correction du calcul de Mm
    Mm_inv = Wetoh / Metoh + Walcene / Malcene + Weau / Meau
    Mm = 1 / Mm_inv  # kg/mol
    # Densité du mélange
    Rho = P * Mm / (R * T)  # kg/m³
    # Utilisation de Rho dans Re_p
    Re_p = Rho * Uv0 * Dp / Mug
    Pr = Mug * Cpm / Lambdag
    Nu =(100*Dp)/Lambdag  #= 2 + 1.8 * np.sqrt(Re_p) * Pr**(1/3) #une des deux equations est utilisable et h=h_int =100 W/m"/K
    h =100 # =(Lambdag / Dp) * Nu     #on a remplacé h par h_int
    Sfm = numerator.sum() #somme de wi sur Mi
    # Coefficient massique
    coeff_massique = (1 - Eps) / (Eps * Rho * Uv0)

    # Initialisation du vecteur de dérivées
    res = np.zeros(5)

    # Calcul des dérivées
    res[0] = -coeff_massique * r * Metoh / (Rho * Sfm)  # Bilan massique sur l'éthanol
    res[1] = coeff_massique * r * Malcene / (Rho * Sfm)  # Bilan massique sur l'éthène
    res[2] = coeff_massique * r * Meau / (Rho * Sfm)     # Bilan massique sur l'eau
    Q_reaction = (1 - Eps) * r * DrH / (Eps * Rho * Uv0)  # Chaleur de réaction (endo)
    Q_heat_transfer = (4 * h / Dr) * (Tp - T)
    res[3] = (Q_heat_transfer + Q_reaction) / (Rho * Uv0 * Cpm)  # Bilan énergétique
    # Perte de charge 
    friction_factor = (1.75 + 150 * (1 - Eps) * Mug / (Rho * Uv0 * Dp))
    res[4] = - (1 - Eps) / (Eps**3 * Dp) * (Rho * Uv0**2) * friction_factor  # Perte de charge

    return res

# Conditions initiales
E0 = [Wetoh0, Walcene0, Weau0, T0, P0]

# Définition de la plage de z pour l'intégration
Lmax = 1  # Longueur maximale du réacteur en mètres 
z_span = (0, Lmax)

# Appel au solveur d'EDO sans événement (comme avant)
sol = solve_ivp(
    fun=solveur,
    t_span=z_span,
    y0=E0,
    method='BDF',  # Méthode pour les problèmes raides
    rtol=1e-6,
    atol=1e-8
)

# Extraction des résultats
z_sol = sol.t
Wetoh_sol = sol.y[0]
Walcene_sol = sol.y[1]
Weau_sol = sol.y[2]
T_sol = sol.y[3]
P_sol = sol.y[4]  # Pression variable si la perte de charge est calculée

# Recalcul de Mm_inv et Mm à chaque position z
Mm_inv_sol = Wetoh_sol / Metoh + Walcene_sol / Malcene + Weau_sol / Meau
Mm_sol = 1 / Mm_inv_sol  # kg/mol

# Recalcul de la densité à chaque position z
Rho_sol = (P_sol * Mm_sol) / (R * T_sol)  # kg/m³

# Calcul des concentrations massiques
C_mass_etoh = Rho_sol * Wetoh_sol  # kg/m³
C_mass_alcene = Rho_sol * Walcene_sol  # kg/m³
C_mass_eau = Rho_sol * Weau_sol  # kg/m³

# Tracé des concentrations massiques en fonction de z
plt.figure(figsize=(10, 7))
plt.plot(z_sol, C_mass_etoh, label='Concentration massique d\'éthanol')
plt.plot(z_sol, C_mass_alcene, label='Concentration massique d\'éthène')
plt.plot(z_sol, C_mass_eau, label='Concentration massique d\'eau')
plt.xlabel('Position z (m)')
plt.ylabel('Concentration massique (kg/m³)')
plt.title('Profil des concentrations massiques le long du réacteur')
plt.xlim(0, 2)  # Limite de l'axe z à 1 m
plt.legend()
plt.grid(True)
plt.show()

# Tracé de la température en fonction de z
plt.figure(figsize=(10, 6))
plt.plot(z_sol, T_sol, label='Température du fluide')
plt.xlabel('Position z (m)')
plt.ylabel('Température (K)')
plt.title('Profil de température le long du réacteur')
plt.xlim(0, 1)      # Limite de l'axe z à 1 m
plt.ylim(None, 1000)  # Limite supérieure de la température à 1000 K
plt.legend()
plt.grid(True)
plt.show()