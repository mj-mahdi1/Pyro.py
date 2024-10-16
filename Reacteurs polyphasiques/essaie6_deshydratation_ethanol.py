
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
"""si Uv0 est de 3m/s la perte de charge devient tres grande et deltaP est de 1,8bar, au dessous de 0,4 la conservation de masse n est pas respectée"""
Dr = 0.1     # Diamètre de tube en m #0,1
"""pour l'instant, aucun effet de changement de diamètre de tube sur la simulation"""
Qm0 = 27.77     # Débit massique initial d'éthanol en kg/s #Fixé

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
def odestiff(z, E):
    # Extraire les 5 variables d'état dans un vecteur ligne appelé E
    Wetoh, Walcene, Weau, T, P = E

    # Calculer fractions molaires yi
    numerator = np.array([Wetoh / Metoh, Walcene / Malcene, Weau / Meau])
    denominator = np.sum(numerator)
    y_etoh = numerator[0] / denominator
    """fractions molaires y_alcene et y_eau sont toujours égales à l'inverse pour chaque terme"""

    # Calculer concentration molaire totale
    C_tot = P / (R * T)  # mol/m³
    # Calculer concentration molaire EtOH
    C_etoh = y_etoh * C_tot  # mol/m³
    # Calculer constante de vitesse
    T_ref = 673.15  # Température de référence pour notre K0 en Kelvin
    
    k = K0 * np.exp(-Ea / R * ((1 / T) - (1 / T_ref))) 
    k = min(k, 1000)  # Arhenus biponctuel ; on mets une limite supérieur pour éviter les problemes de calcul

    # vitesse de réaction
    r = k * C_etoh  # mol/(m³·s)

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
z_span = (0, Lmax)

# Appel au solveur d'EDOs
sol = solve_ivp(
    fun=odestiff,
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
Mm_inv_sol = (Wetoh_sol / Metoh) + (Walcene_sol / Malcene) + (Weau_sol / Meau)
Mm_sol = 1 / Mm_inv_sol  # kg/mol 

# Recalcul de la densité à chaque position z
Rho_sol = (P_sol * Mm_sol) / (R * T_sol)  # kg/m³

# Calcul des concentrations massiques
C_mass_etoh = Rho_sol * Wetoh_sol  # kg/m³
C_mass_alcene = Rho_sol * Walcene_sol  # kg/m³
C_mass_eau = Rho_sol * Weau_sol  # kg/m³

# Affichage des résultats
df = pd.DataFrame({
    "z": z_sol,  
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
plt.plot(z_sol, C_mass_etoh, label='Concentration massique d\'éthanol')
plt.plot(z_sol, C_mass_alcene, label='Concentration massique d\'éthène')
plt.plot(z_sol, C_mass_eau, label='Concentration massique d\'eau')
plt.xlabel('Position z (m)')
plt.ylabel('Concentration massique (kg/m³)')
plt.title('Profil des concentrations massiques le long du réacteur')
plt.xlim(0, Lmax)  
plt.legend()
plt.grid(True)
plt.show()

# Tracé de la température en fonction de z
plt.figure(figsize=(11, 8))
plt.plot(z_sol, T_sol, label='Température du fluide')
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
plt.plot(z_sol, Wetoh_sol, label='Ethanol')
plt.plot(z_sol, Walcene_sol, label='Ethène')
plt.plot(z_sol, Weau_sol, label='Eau')
plt.xlabel('Position z (m)')
plt.ylabel('Fraction massique')
plt.title('Fractions massiques au long du réacteur')
plt.legend()
plt.grid(True)
plt.show()

# Tracé de la pression en fonction de z
plt.figure(figsize=(11, 8))
plt.plot(z_sol, P_sol, label='Pression le long du réacteur')
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
plt.plot(z_sol, sum_Wi, label='Somme de fractions massiques')
plt.xlabel('Position z (m)')
plt.ylabel('Somme des fractions massiques')
plt.title('Vérification de conservation de matière')
plt.legend()
plt.grid(True)
plt.show()


















"""import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt  # Pour les tracés
import pandas as pd

# Constantes physiques
R = 8.314  # Constante des gaz parfaits en J/(mol·K)

# Paramètres invariants
Rpore = 10e-9      # Rayon de pore en m (fixé à 10 nm)
Rhop = 1100        # Masse volumique du catalyseur en kg/m³
Mug = 17.5e-6      # Viscosité dynamique du gaz supposée constante en Pa·s
P0 = 2e5           # Pression initiale en Pa en entrée de tube
T0 = 900           # Température initiale en K
Dp = 0.003         # Diamètre des particules en m (ajusté à 3 mm)
Metoh = 0.046      # Masse molaire de l'éthanol en kg/mol
Malcene = 0.028    # Masse molaire de l'éthène en kg/mol
Meau = 0.018       # Masse molaire de l'eau en kg/mol
Eps = 0.4          # Porosité du lit
Ea = 150e3         # Énergie d'activation en J/mol
DrH = 49.1e3       # Enthalpie de réaction à 900 K en J/mol
Tp = 900           # Température de paroi considérée fixe, en K

# Capacités calorifiques massiques (en J/kg·K)
Cpetoh = 2953      # Ethanol
Cpalcene = 3186    # Ethylène
Cpeau = 2217       # Eau

Lambdag = 0.025    # Conductivité thermique en W/(m·K)
K0 = 0.66          # Constante cinétique à 400°C (ou 673.15 K) en s⁻¹

# Paramètres ajustés
Uv0 = 0.45         # Vitesse superficielle en m/s
Dr = 0.1           # Diamètre de tube en m
Qm0 = 27.77        # Débit massique initial d'éthanol en kg/s

# Fractions massiques initiales
Wetoh0 = 1.0
Walcene0 = 0.0
Weau0 = 0.0

# Calcul de la masse molaire moyenne initiale du mélange
M0 = Wetoh0 * Metoh + Walcene0 * Malcene + Weau0 * Meau
Rho0 = P0 * M0 / (R * T0)  # Masse volumique initiale du gaz en kg/m³
Qv0 = Qm0 / Rho0           # Débit volumique initial en m³/s

# Définition de la fonction pour le solveur
def odestiff(z, E):
    Wetoh, Walcene, Weau, T, P = E

    # Calcul des fractions molaires
    numerator = np.array([Wetoh / Metoh, Walcene / Malcene, Weau / Meau])
    denominator = np.sum(numerator)
    y_etoh = numerator[0] / denominator

    # Concentration molaire totale
    C_tot = P / (R * T)  # mol/m³
    C_etoh = y_etoh * C_tot  # mol/m³

    # Constante cinétique selon Arrhenius biponctuel
    T_ref = 673.15  # Température de référence
    k = K0 * np.exp(-Ea / R * ((1 / T) - (1 / T_ref)))
    k = min(k, 0.1)

    # Vitesse de réaction
    r = k * C_etoh  # mol/(m³·s)

    # Capacité calorifique massique moyenne
    Cpm = Wetoh * Cpetoh + Walcene * Cpalcene + Weau * Cpeau

    # Masse molaire moyenne
    Mm_inv = Wetoh / Metoh + Walcene / Malcene + Weau / Meau
    Mm = 1 / Mm_inv  # kg/mol

    # Densité du mélange
    Rho = P * Mm / (R * T)  # kg/m³

    # Nombre de Reynolds et Prandtl
    Re_p = (Rho * Uv0 * Dp) / Mug
    Pr = Mug * Cpm / Lambdag

    # Nombre de Nusselt et coefficient de transfert de chaleur
    Nu = 2 + 1.8 * np.sqrt(Re_p) * Pr**(1/3)
    h = (Lambdag / Dp) * Nu  # Coefficient de transfert de chaleur interne

    # Bilan de matière
    res = np.zeros(5)
    coeff_massique = (1 - Eps) / (Eps * Rho * Uv0)
    res[0] = -coeff_massique * r * Metoh  # Ethanol
    res[1] = coeff_massique * r * Malcene  # Ethylene
    res[2] = coeff_massique * r * Meau  # Eau

    # Bilan d'énergie
    Q_reaction = -(1 - Eps) * r * DrH / Mm
    Q_heat_transfer = (4 * h / Dr) * (Tp - T)
    res[3] = (Q_heat_transfer + Q_reaction) / (Rho * Uv0 * Cpm)

    # Bilan quantité de mouvement (perte de charge)
    friction_factor = 1.75 + 150 * (1 - Eps) * Mug / (Rho * Uv0 * Dp)
    res[4] = -(1 - Eps) * friction_factor * Rho * Uv0**2 / (Eps**3 * Dp)

    return res

# Conditions initiales
E0 = [Wetoh0, Walcene0, Weau0, T0, P0]

# Plage de z pour l'intégration
Lmax = 10  # Longueur maximale du réacteur en mètres
z_span = (0, Lmax)

# Résolution de l'ODE
sol = solve_ivp(
    fun=odestiff,
    t_span=z_span,
    y0=E0,
    method='BDF',
    rtol=1e-6,
    atol=1e-8
)

# Extraction des résultats
z_sol = sol.t
Wetoh_sol = sol.y[0]
Walcene_sol = sol.y[1]
Weau_sol = sol.y[2]
T_sol = sol.y[3]
P_sol = sol.y[4]

# Masse molaire moyenne à chaque position z
Mm_inv_sol = (Wetoh_sol / Metoh) + (Walcene_sol / Malcene) + (Weau_sol / Meau)
Mm_sol = 1 / Mm_inv_sol  # kg/mol

# Densité à chaque position z
Rho_sol = (P_sol * Mm_sol) / (R * T_sol)  # kg/m³

# Concentrations massiques
C_mass_etoh = Rho_sol * Wetoh_sol  # kg/m³
C_mass_alcene = Rho_sol * Walcene_sol  # kg/m³
C_mass_eau = Rho_sol * Weau_sol  # kg/m³

# Affichage des résultats sous forme de tableau
df = pd.DataFrame({
    "z": z_sol,
    "Wetoh": Wetoh_sol,
    "Walcene": Walcene_sol,
    "Weau": Weau_sol,
    "T": T_sol,
    "P": P_sol,
    "Masse Molaire": Mm_sol
})
print(df)

# Tracés
plt.figure(figsize=(11, 8))
plt.plot(z_sol, C_mass_etoh, label='Concentration massique d\'éthanol')
plt.plot(z_sol, C_mass_alcene, label='Concentration massique d\'éthène')
plt.plot(z_sol, C_mass_eau, label='Concentration massique d\'eau')
plt.xlabel('Position z (m)')
plt.ylabel('Concentration massique (kg/m³)')
plt.title('Profil des concentrations massiques le long du réacteur')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(11, 8))
plt.plot(z_sol, T_sol, label='Température du fluide')
plt.xlabel('Position z (m)')
plt.ylabel('Température (K)')
plt.title('Profil de température le long du réacteur')
plt.legend()
plt.grid(True)
plt.show()

plt.figure(figsize=(11, 8))
plt.plot(z_sol, Wetoh_sol, label='Ethanol')
plt.plot(z_sol, Walcene_sol, label='Ethene')
plt.plot(z_sol, Weau_sol, label='Eau')
plt.xlabel('Position z (m)')
plt.ylabel('Fraction massique')
plt.title('Fractions massiques le long du réacteur')
plt.legend()
plt.grid(True)
"""