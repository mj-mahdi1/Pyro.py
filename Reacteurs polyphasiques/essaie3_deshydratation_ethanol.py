import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt  # Pour les tracés

# Toutes les données sont en unités SI

# Paramètres invariants
Rpore = 10e-9      # Rayon de pore en m (fixé à 10 nm)
Rhop = 1100        # Masse volumique du catalyseur en kg/m3
Mug = 17.5e-6      # Viscosité dynamique du gaz supposée constante en Pa.s
R = 8.314          # Constante des gaz parfaits en J/mol/K
P0 = 2e5           # Pression initiale en Pa
T0 = 900           # Température initiale en K
Dp = 0.001         # Diamètre d'une sphère de particule en m
Metoh = 0.046      # Masse molaire de l'éthanol en kg/mol
Malcene = 0.028    # Masse molaire de l'éthène en kg/mol
Meau = 0.018       # Masse molaire de l'eau en kg/mol
Eps = 0.4          # Porosité du lit
Ea = 150e3         # Énergie d'activation en J/mol
DrH = 49.1e3       # Enthalpie de réaction à 900K en J/mol
Tp = 900           # Température de paroi considérée fixe, en K
Cpetoh = 2953      # Capacité calorifique massique de l'éthanol en J/kg/K
Cpalcene = 3186    # Capacité calorifique massique de l'éthène en J/kg/K
Cpeau = 2217       # Capacité calorifique massique de l'eau en J/kg/K
Lambdag = 0.025    # Conductivité thermique en W/m/K
K0 = 0.66          # Constante cinétique à 400°C (ou 673,15 K) en s-1

# Paramètres fixes
Uv0 = 1.0   # Vitesse superficielle en m/s
Dr = 0.01   # Diamètre de tube en m

# Fractions massiques initiales
Wetoh0 = 1.0        # Fraction massique initiale d'éthanol (1)
Walcene0 = 0.0      # Fraction massique initiale d'éthène (0)
Weau0 = 0.0         # Fraction massique initiale d'eau (0)

# Calcul de la masse volumique initiale du gaz
M0 = Metoh * Wetoh0 + Malcene * Walcene0 + Meau * Weau0
Rho0 = P0 * M0 / (R * T0)  # Masse volumique initiale du gaz en kg/m3

# Définition de la fonction solveur
def solveur(z, E):
    # Extraction des variables d'état
    Wetoh, Walcene, Weau, T, P = E

    # Calcul des fractions molaires
    numerator = np.array([Wetoh / Metoh, Walcene / Malcene, Weau / Meau])
    denominator = np.sum(numerator)
    y_etoh = numerator[0] / denominator  # Fraction molaire de l'éthanol
    y_alcene = numerator[1] / denominator  # Fraction molaire de l'éthène
    y_eau = numerator[2] / denominator     # Fraction molaire de l'eau

    # Calcul de la concentration totale
    C_tot = P / (R * T)  # mol/m3 avec P étant la pression totale
    # Concentrations molaires des composants
    C_etoh = y_etoh * C_tot    # mol/m3 concentration molaire EtOH
    C_alcene = y_alcene * C_tot  # mol/m3 concentration molaire Éthène
    C_eau = y_eau * C_tot      # mol/m3 concentration molaire Eau

    # Calcul de k avec la bonne température de référence
    T_ref = 673.15  # Température de référence pour K0 en K (400 °C)
    k = K0 * np.exp(-Ea / R * (1 / T - 1 / T_ref))  # Constante cinétique à la température T

    # Taux de réaction
    r = k * C_etoh  # mol/(m3·s)

    # Calculs intermédiaires
    Cpm = Wetoh * Cpetoh + Walcene * Cpalcene + Weau * Cpeau
    # Correction du calcul de Mm
    Mm_inv = Wetoh / Metoh + Walcene / Malcene + Weau / Meau
    Mm = 1 / Mm_inv
    # Densité du mélange
    Rho = P * Mm / (R * T)
    # Utilisation de Rho dans Re_p
    Re_p = Rho * Uv0 * Dp / Mug
    Pr = Mug * Cpm / Lambdag
    Nu = 2 + 1.8 * np.sqrt(Re_p) * (Pr)**(1/3)
    h = (Lambdag / Dp) * Nu
    Sfm = numerator.sum()
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
    # Perte de charge (non nécessaire si P est supposé constant)
    res[4] = 0  # On suppose la pression constante

    return res

# Conditions initiales
E0 = [Wetoh0, Walcene0, Weau0, T0, P0]

# Définition de la plage de z pour l'intégration
Lmax = 10  # Longueur maximale du réacteur en mètres
z_span = (0, Lmax)

# Appel au solveur d'EDO sans événement
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

# Calcul des fractions molaires à chaque point
numerator = np.array([Wetoh_sol / Metoh, Walcene_sol / Malcene, Weau_sol / Meau])
denominator = np.sum(numerator, axis=0)
y_etoh_sol = numerator[0] / denominator
y_alcene_sol = numerator[1] / denominator
y_eau_sol = numerator[2] / denominator

# Calcul de la concentration totale à chaque point
C_tot_sol = P0 / (R * T_sol)  # mol/m3, en supposant P constant

# Calcul des concentrations molaires à chaque point
C_etoh_sol = y_etoh_sol * C_tot_sol
C_alcene_sol = y_alcene_sol * C_tot_sol
C_eau_sol = y_eau_sol * C_tot_sol

# Tracé des concentrations en fonction de z
plt.figure(figsize=(10, 6))
plt.plot(z_sol, C_etoh_sol, label='Concentration d\'éthanol')
plt.plot(z_sol, C_alcene_sol, label='Concentration d\'éthène')
plt.plot(z_sol, C_eau_sol, label='Concentration d\'eau')
plt.xlabel('Position z (m)')
plt.ylabel('Concentration molaire (mol/m³)')
plt.title('Profil des concentrations molaires le long du réacteur')
plt.legend()
plt.grid(True)
plt.show()

# Tracé de la température en fonction de z
plt.figure(figsize=(10, 6))
plt.plot(z_sol, T_sol, label='Température du fluide')
plt.xlabel('Position z (m)')
plt.ylabel('Température (K)')
plt.title('Profil de température le long du réacteur')
plt.legend()
plt.grid(True)
plt.show()
