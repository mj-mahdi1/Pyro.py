import numpy as np
from scipy.integrate import solve_ivp

# Toutes les données sont en unités SI

# Paramètres invariants
Rpore = 10e-9      # Rayon de pore en m (fixé à 10 nm)
Rhop = 1100        # Masse volumique du catalyseur en kg/m3
Mug = 17.5e-6      # Viscosité dynamique du gaz supposée constante en Pa.s
R = 8.314          # Constante des gaz parfaits en J/mol/K
M0 = 0.046         # Masse molaire initiale du gaz en kg/mol
P0 = 2e5           # Pression initiale en Pa
T0 = 900           # Température initiale en K
Rho0 = P0 * M0 / (R * T0)  # Masse volumique initiale du gaz en kg/m3
Dp = 0.001         # Diamètre d'une sphère de particule en m
Metoh = 0.046      # Masse molaire de l'éthanol en kg/mol
Malcene = 0.028    # Masse molaire de l'éthène en kg/mol
Meau = 0.018       # Masse molaire de l'eau en kg/mol
Qm0 = 10e4 / 3600  # Débit massique d'éthanol en kg/s
Qv0 = Qm0 / Rho0   # Débit volumique d'éthanol en m3/s
Eps = 0.4          # Porosité du lit
Ea = 150e3         # Énergie d'activation en J/mol
Wetoh0 = 1         # Fraction massique initiale d'éthanol
Walcene0 = 0       # Fraction massique initiale d'éthène
Weau0 = 0          # Fraction massique initiale d'eau
DrH = 49.1         # Enthalpie de réaction à 900K en J/mol
Tp = 900           # Température de paroi considérée fixe, en K
Cpetoh = 2953      # Capacité calorifique massique de l'éthanol en J/kg/K
Cpalcene = 3186    # Capacité calorifique massique de l'éthène en J/kg/K
Cpeau = 2217       # Capacité calorifique massique de l'eau en J/kg/K
Lambdag = 0.025    # Conductivité thermique en W/m/K
K0 = 0.66          # Constante cinétique à 400°C (ou 673,15k) en s-1 ____ k calculée à 900k = 566,5 s-1 
#ode23s est le solveur pour ode raides
Uv0 = 1 # en m/s  
Dr = 0.02 #en m
"""# Calcul de la concentration totale
C_tot = P / (R * T)  # mol/m3 avec P etant la pression totale
# Concentration molaire de l'éthanol
C_etoh = y_etoh * C_tot        # mol/m3 concentration molaire EtOH
C_alcene = y_alcene * C_tot    # mol/m3 concentration molaire Alcene
C_eau = y_eau * C_tot          # mol/m3 concentration molaire Eau"""

# Paramètres variables
"""Dr = np.arange(0.02, 0.201, 0.01)  # Diamètres de tube possibles en m
Uv0 = np.arange(0.01, 10.01, 0.01)  # Vitesses superficielles en m/s"""

"""ld = len(Dr)"""
"""lu = len(Uv0)"""

# Nombre de tubes pour chaque combinaison de Uv0 et Dr
N = np.zeros((ld, lu))
"""for i in range(ld):
    for j in range(lu):
        N[i, j] = 4 * Qv0 / (np.pi * Dr[i]**2 * Uv0[j])

N = np.ceil(N)  # Arrondi à l'entier supérieur"""

# Initialisation des matrices de résultats
"""H = np.zeros((ld, lu))      # Hauteur du réacteur pour chaque combinaison
DP = np.zeros((ld, lu))     # Perte de charge pour chaque combinaison"""

# Définition de la fonction d'événement pour arrêter l'intégration
"""def event_Wetoh(z, E):
    Wetoh = E[0]
    return Wetoh - 0.05  # Arrêter lorsque Wetoh atteint 0.05"""

"""event_Wetoh.terminal = True       # Arrêter l'intégration
event_Wetoh.direction = -1        # Détection lorsque Wetoh descend en dessous de 0.05"""

# Définition de la fonction solveur
def solveur(z, E):
    # Extraction des variables d'état
    Wetoh, Walcene, Weau, T, P = E

    # Masse molaire des espèces
    # Metoh = 0.046  # kg/mol
    # Malcene = 0.028  # kg/mol
    # Meau = 0.018  # kg/mol

    # Calcul des fractions molaires
    """ numerator = np.array([Wetoh / Metoh, Walcene / Malcene, Weau / Meau])
    denominator = np.sum(numerator)
    y_etoh = numerator[0] / denominator  # Fraction molaire de l'éthanol
    y_alcene = numerator[1] / denominator
    y_eau = numerator[2] / denominator


    # Calcul de la concentration totale
    C_tot = P / (R * T)  # mol/m3 avec P etant la pression totale
    # Concentration molaire de l'éthanol
    C_etoh = y_etoh * C_tot        # mol/m3 concentration molaire EtOH
    C_alcene = y_alcene * C_tot    # mol/m3 concentration molaire Alcene
    C_eau = y_eau * C_tot          # mol/m3 concentration molaire Eau"""  


    # Reinitialisation de k avec la nouvelle relation
    
    T0 = 900  # temperature initiale en K
    k = 566.5 #k = K0 * np.exp((-Ea / R )* (1 / T - 1 / T0))  #constante cinetique à la temperature operatoire en s-1 selon Arhenus biponctuel (hyp que ecart T pas significatif)

    # Taux de réaction
    numerator = np.array([Wetoh / Metoh, Walcene / Malcene, Weau / Meau])
    denominator = np.sum(numerator)

    C_tot = P / (R * T)  # mol/m3 avec P etant la pression totale

    y_etoh = numerator[0] / denominator #calcul de fraction molaire de chaque composant
    C_etoh = y_etoh * C_tot        # mol/m3 concentration molaire EtOH
    r = k * C_etoh # mol/(m3·s)


    # Calculs intermédiaires
    Cpm = Wetoh * Cpetoh + Walcene * Cpalcene + Weau * Cpeau
    Mm = Metoh * Wetoh + Malcene * Walcene + Meau * Weau
    Rho = P * Mm / (R * T)
    Re_p = Rho0 * Uvj * Dp / Mug
    Pr = Mug * Cpm / Lambdag
    Nu = 2 + 1.8 * np.sqrt(Re_p) * (Pr)**(1/3)
    h = (Lambdag / Dp) * Nu
    Sfm = Wetoh / Metoh + Walcene / Malcene + Weau / Meau
    k = K0 * np.exp(-Ea / R * (1 / T - 1 / T0))
    r = k * Wetoh * P  # Taux de réaction

    # Initialisation du vecteur de dérivées
    res = np.zeros(5)

    # Calcul des dérivées
    coeff_massique = (1 - Eps) / (Eps * Rho0 * Uvj)
    res[0] = -coeff_massique * r / (Metoh * Sfm)           # Bilan massique sur l'éthanol
    res[1] = coeff_massique * r / (Malcene * Sfm)          # Bilan massique sur l'éthène
    res[2] = coeff_massique * r / (Meau * Sfm)             # Bilan massique sur l'eau
    Q_reaction = +(1 - Eps) * r * DrH / (Eps * Rho0 * Uvj) #on a changé le signe au positif puisque endo
    Q_heat_transfer = (4 * h / Dri) * (Tp - T)
    res[3] = (Q_heat_transfer + Q_reaction) / (Rho0 * Uvj * Cpm)  # Bilan énergétique
    res[4] = - (1 - Eps) / (Eps**3 * Dp) * (Rho * Uvj**2) * (1.75 + 150 * (1 - Eps) * Mug / (Rho0 * Uvj * Dp))  # Perte de charge

    return res

# Définition de la longueur maximale du réacteur
Lmax = 10  # en mètres

# Parcours de toutes les combinaisons de Dr et Uv0
for i in range(ld):
    for j in range(lu):
        Uvj = Uv0[j]
        Dri = Dr[i]
        E0 = [Wetoh0, Walcene0, Weau0, T0, P0]  # Conditions initiales

        # Appel au solveur d'EDO avec événement
        sol = solve_ivp(
            fun=lambda z, E: solveur(z, E, Uvj, Dri),
            t_span=(0, Lmax),
            y0=E0,
            events=event_Wetoh,
            method='BDF',  # Méthode pour les problèmes raides
            rtol=1e-6,
            atol=1e-8
        )

        # Vérification si l'événement a été détecté
        if sol.status == 1 and len(sol.t_events[0]) > 0:
            H[i, j] = sol.t_events[0][0]  # Hauteur du réacteur nécessaire
            final_P = sol.y[4, sol.t >= H[i, j]][0]
        else:
            # L'événement n'a pas été détecté dans Lmax
            H[i, j] = Lmax
            final_P = sol.y[4, -1]
            print(f"Pour Dr = {Dri:.3f} m et Uv0 = {Uvj:.3f} m/s, Wetoh n'a pas atteint 0.05 dans Lmax")

        DP[i, j] = P0 - final_P  # Perte de charge

# Calcul du volume catalytique Vcata pour chaque combinaison
Vcata = (1 - Eps) * N * H * (np.pi / 4) * (Dr[:, np.newaxis] ** 2)

# Sélection des configurations éligibles (DP < 1 bar)
Vcatae = []

for i in range(ld):
    for j in range(lu):
        if DP[i, j] < 1e5:  # 1 bar en Pa
            Vcatae.append([Vcata[i, j], Dr[i], Uv0[j], DP[i, j]])

if Vcatae:
    Vcatae = np.array(Vcatae).T  # Conversion en array numpy pour faciliter l'accès
    print("\nConfigurations éligibles (Vcata en m³, Dr en m, Uv0 en m/s, DP en Pa) :")
    print(Vcatae)

    # Recherche de la configuration avec le volume catalytique minimal
    Vcatamin_index = np.argmin(Vcatae[0])
    Re = Vcatae[:, Vcatamin_index]
    print("\nConfiguration optimale (Vcata_min en m³, Dr en m, Uv0 en m/s, DP en Pa) :")
    print(Re)
else:
    print("Aucun résultat, veuillez ajuster vos paramètres.")




