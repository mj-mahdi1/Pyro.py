import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt  # Pour les tracés
import pandas as pd
# Toutes les données sont en unités SI
"""Objectif: Vcata mini, (Lmax, N, Uv0 Dr) à rajouter kp apparent, """
"""Objectif: Vcata mini, (Lmax, N, Uv0 Dr) à rajouter kp apparent, ce qui reste a faire le calcul de h en paroi (prendre h constant au long du tube) dp le plus petit pour avoir la vitesse d ecoulement la plus grande, on economise du temps : dp/dr=15 pour eviter l effet de paroi (ca cree des chema preferentiels en paroi qui perturbe (dispersion axiale en 2D par ex :dp/dr=15 c a d 15 particule dans un diametre de tube infintesimal))"""
"""uv plus je limite la vitesse de transport, et hint devient plus efficace, donc ca me permet de travailler"""
"""choisr Dr et Uv est le plus important, tubes trop gros, surface d echange thermique tres faible"""
"""rajouter une boucle d optimisation et fournier les profils et expliquer les causes effets et montrer qu on a utiliser l outil pour optimiser le design. finir """
"""on commonce par : dr:1,5cm , n fixe dp on fixe dr aux plus peits et on trouve uv pour avoir la conversion souhaiter et s assurer que cest compatible avec delta p de 1 bar, tp a 900"""
"""2 eme scenario: d minimiser le nombre de tube est interessant de point de vue realistique les tubes de 1.5mm ca existe pour une reaction endo comme ca par contre, pour une reaction exo cest plus realistic des tubes de 1cm pour que l echange soit possible thermiquement, criteres economiques doivent etre justifier"""

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
"""A definir A,B,C,D"""
# Paramètres fixes ajustés
Uv0 = 0.5    # Vitesse superficielle en m/s 
"""si Uv0 est moins que 0.45, il n y a plus de conservation de matiere dans ces conditions"""
Dr = 0.1     # Diamètre de tube en m 
"""pour l'instant, aucun effet de changement de diametre de tube sur la simulation"""
Qm0 = 27.77     # Débit massique initial d'éthanol en kg/s #Fixé

# Fractions massiques initiales
Wetoh0 = 1.0        # Fraction massique initiale d'éthanol 
Walcene0 = 0.0      # Fraction massique initiale d'éthène 
Weau0 = 0.0         # Fraction massique initiale d'eau 

# Calcul de la masse volumique initiale du gaz
M0 = Metoh * Wetoh0 + Malcene * Walcene0 + Meau * Weau0 #Masse molaire moyenne de melange de gaz
Rho0 = P0 * M0 / (R * T0)  # Masse volumique initiale du gaz en kg/m³

# Débit volumique initial
Qv0 = Qm0 / Rho0   # m³/s

# Définition de la fonction du solveur
def odestiff(z, E):
    # Extraire les 5 variables d'etat dans un vecteur ligne appelé E
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
    T_ref = 673.15  #Temperature de référence pour notre K0 en Kelvin
    
    k = K0 * np.exp(-Ea / R * ((1 / T) - (1 / T_ref))) ; k = min(k, 0.1) #Arhenus biponctuel

    # vitesse de reaction
    r = k * C_etoh  # mol/(m³·s)

    # Calculs intermédiaires

    Cpm = Wetoh * Cpetoh + Walcene * Cpalcene + Weau * Cpeau 
    """Cp moyen de phase gazeuse il faut utiliser Cpi(T) =A+BT+CT^2+DT^3 pour l evolution de Cpi"""
    # Masse molaire moyenne
    Mm_inv = Wetoh / Metoh + Walcene / Malcene + Weau / Meau #formule ici: 1/Mm = (Wetoh/Metoh + Walcene/Malcene + Weau/Meau)
    Mm = 1 / Mm_inv  # kg/mol
    # Mm = ((Wetoh * Metoh)/y_etoh + (Walcene / Malcene)/(1-y_etoh) + (Weau / Meau)/(1-y_etoh))
    # Mm_inv = 1/Mm
    # Densité du mélange
    
    Rho = P * Mm / (R * T)  # kg/m³
    # Reynolds et Prandtl numbers
    Re_p = (Rho * Uv0 * Dp) / Mug
    Pr = Mug * Cpm / Lambdag
    # Nusselt number et coefficient de transfert de chaleur
    Nu = 2 + 1.8 * np.sqrt(Re_p) * Pr**(1/3)
    #h = (Lambdag / Dp) * Nu  #coefficient de transfert de chaleur interne 
    """h =0 pour qu on soit en adiabétique et que la descente en T initiale soit vue"""
    h = 0 
    # coefficient de transfert de matière
    coeff_massique = (1 - Eps) / (Eps * Rho * Uv0) #rho est la seule variable dans cette formule

    # Initializer vecteur dérivé
    res = np.zeros(5) #crée une matrice qui se remplie

    # Bilan de matière
    res[0] = -coeff_massique * r * Metoh    # Ethanol
    res[1] =  coeff_massique * r * Malcene  # Ethylene
    res[2] =  coeff_massique * r * Meau     # Eau

    
    # Bilan d'energie
    Q_reaction =  -(1 - Eps) * r * DrH / ((Wetoh/Metoh)+(Walcene/Malcene)+(Weau/Meau))  # r est la vitesse de reaction explicité précédemment
    Q_heat_transfer = (4 * h / Dr) * (Tp - T)
    res[3] = (Q_heat_transfer + Q_reaction) / (Rho0 * Uv0 * Cpm) # rho0 dans l equation BE 1,2285


    # Bilan quantité de mouvement
    friction_factor = 1.75 + 150 * (1 - Eps) * Mug / (Rho * Uv0 * Dp)
    res[4] = - (1 - Eps) * friction_factor * Rho * Uv0**2 / (Eps**3 * Dp)

    return res


# Conditions initiales
E0 = [Wetoh0, Walcene0, Weau0, T0, P0]

# Définition de la plage de z pour l'intégration
Lmax = 10  # Longueur maximale du réacteur en mètres 
z_span = (0, Lmax)

# Appel au solveur d'EDOs sans événements pour l'instant, Solveur qui utilise un pas adaptatif en fx des points de tolérence
sol = solve_ivp(
    fun=odestiff,
    t_span=z_span, #parce que normalement t_span est utilisée
    y0=E0,
    method='BDF',  # Méthode pour les problèmes raides Voir choix de solveur
    rtol=1e-6,
    atol=1e-8   #parametres de tolérence
)

# Extraction des résultats (paramètres calculés par le solveur)
z_sol = sol.t  #champ sol.t contient tous les points temporelles dont les paramètres ont été évalués (ici on fait z_sol au lieu de t_sol car spacial)
Wetoh_sol = sol.y[0]
Walcene_sol = sol.y[1]
Weau_sol = sol.y[2]
T_sol = sol.y[3]
P_sol = sol.y[4]  # Pression variable si la perte de charge est calculée

# Recalcul de Mm_inv et Mm à chaque position z
Mm_inv_sol = (Wetoh_sol / Metoh) + (Walcene_sol / Malcene) + (Weau_sol / Meau)
Mm_sol = 1 / Mm_inv_sol  # kg/mol #Masse molaire moyenne recalculée

# Recalcul de la densité à chaque position z
Rho_sol = (P_sol * Mm_sol) / (R * T_sol)  # kg/m³

# Calcul des concentrations massiques
C_mass_etoh = Rho_sol * Wetoh_sol  # kg/m³
C_mass_alcene = Rho_sol * Walcene_sol  # kg/m³
C_mass_eau = Rho_sol * Weau_sol  # kg/m³

"""Connaitre la longeur de vecteur sol.t sorti"""
# Afficher la longueur de sol.t
#print(f"Longueur de sol.t : {len(sol.t)}") #le nombre de points de z où les équations différentielles ont été résolues.
# Afficher les valeurs de sol.t
#print(f"Valeurs de sol.t : {sol.t}")       # les points exacts où les calculs ont été effectués.



""" Afficher Tableau de parametres calculés sol."""

# Utiliser les valeurs de z calculées par le solveur solve_ivp
z_values = sol.t  # Cela utilise directement les valeurs de z calculées

# Création du DataFrame avec Pandas
df = pd.DataFrame({
    "z": z_values,  # z_values est maintenant défini comme sol.t
    "Wetoh": Wetoh_sol,
    "Walcene": Walcene_sol,
    "Weau": Weau_sol,
    "T": T_sol,
    "P": P_sol,
    "Masse Molaire": Mm_sol
})

# Affichage du DataFrame pour voir les résultats
print(df)


# Tracé des concentrations massiques en fonction de z
plt.figure(figsize=(11, 8))
plt.plot(z_sol, C_mass_etoh, label='Concentration massique d\'éthanol')
plt.plot(z_sol, C_mass_alcene, label='Concentration massique d\'éthène')
plt.plot(z_sol, C_mass_eau, label='Concentration massique d\'eau')
plt.xlabel('Position z (m)')
plt.ylabel('Concentration massique (kg/m³)')
plt.title('Profil des concentrations massiques le long du réacteur')
plt.xlim(0, 10)  # Limite de l'axe z à 1 m
plt.legend()
plt.grid(True)
plt.show()

# Tracé de la température en fonction de z
plt.figure(figsize=(11, 8))
plt.plot(z_sol, T_sol, label='Température du fluide')
plt.xlabel('Position z (m)')
plt.ylabel('Température (K)')
plt.title('Profil de température le long du réacteur')
plt.xlim(0, 10)      # Limite de l'axe z à 1 m
plt.ylim(500, 1500)  # Limite supérieure de la température à 1000 K
plt.legend()
plt.grid(True)
plt.show()


# tracer fractions massiques au long du tube
plt.figure(figsize=(11, 8))
plt.plot(z_sol, Wetoh_sol, label='Ethanol')
plt.plot(z_sol, Walcene_sol, label='Ethene')
plt.plot(z_sol, Weau_sol, label='eau')
plt.xlabel('Position z (m)')
plt.ylabel('fraction massique')
plt.title('fractions massiques au long du réacteur')
plt.legend()
plt.grid(True)
plt.show()

# Tracé de la pression en fonction de z
plt.figure(figsize=(11, 8))
plt.plot(z_sol, P_sol, label='Pression le long du réacteur')
plt.xlabel('Position z (m)')
plt.ylabel('Pression (Pa)')
plt.title('Profil de pression le long du réacteur')
plt.xlim(0, 10)  # Limite de l'axe z à 10 m
plt.legend()
plt.grid(True)
plt.show()


# Verifier la conservation massique
sum_Wi = Wetoh_sol + Walcene_sol + Weau_sol
plt.figure(figsize=(11, 8))
plt.plot(z_sol, sum_Wi, label='Somme de fractions massiques')
plt.xlabel('Position z (m)')
plt.ylabel('somme des fractions massiques')
plt.title('Verification de conservation de matière')
plt.legend()
plt.grid(True)
plt.show()
