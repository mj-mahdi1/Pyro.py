"""#MODEL DEVELOPMENT: Overview
#Step 1: Generalization for Multiple Reactants: To model more complex systems with multiple reactants and different orders, we need to generalize this to handle a list of reactants and their interactions.
#     a) Define the list of molecules involved.
#     b) Set up a data structure to store reactions.
#     c) Modify the ODE function to handle multiple reactions and species.
#     d) Update initial conditions
#     e) Define Time Points and Temperature
#     f) Solve the ODEs using solver odeint.

#Step 2: Model Accuracy: a) the reaction order for each reactant may not be 1 (as assumed here), so this should be parameterized.
#                        b) Real reactions may involve changes in temperature over time, so this should also be considered for more accurate simulations.
#                        c) Incorporating enthalpy and entropy might help in predicting the final composition of the products under varying conditions.
#Step 3: Adding products not only reactants: This code only tracks the reactant's concentration. We can extend it to also track product formation based on stoichiometry and rate laws.
#Step 4: Result Extraction: we'll want to extract meaningful results, like time to complete conversion, yield of the products, and final composition. We’ll add logic to compute these quantities based on the output from the odeint solver.

#Step 1: Generalize the code to multiple reactants and products: expand the batch_reactor function to accept an arbitrary number of reactants, each with its own concentration and reaction order.

#Note 1: Same format of reactions nomination is used to later get inspired by the ligpy.py model. "stoic-coeff_species,stoic-coeff_species2,…" every molecule (for example C4H8O2 is coded as POGPY specie) where	stoic-coeff	is	negative	for	reactants	and	positive	for	products
#Note 2: The names of groups linked to the solid residue start with P, where radical species start with R (radicals linked to the degrading polymer matrix start with PR).
#Note 3: Molecules are coded (e.g., C₄H₈O₂ is coded as POGPY and has a corresponding structure in the appendix)
#Note 4: k1 : Forward reaction rate constant; K-1 Reverse reaction rate constant; Keq = k1/k-1, if greater than 1 then k1 is predominant forward reaction etc., if no information about Keq or delta_G: reaction is considered irreversible.
"""

import scipy as sp
import numpy as np 
import matplotlib.pyplot as plt
from scipy.integrate import odeint
# Define constants

R = 8.314  # J/(mol*K), Perfect gas constant in standard conditions

# Arrhenius equation
def arrhenius_constant(A, Ea, T):
    return A * np.exp(-Ea / (R * T))


# a) Define the list of molecules involved.
species = ['LIGH', 'RLIGM2A', 'OH', 'C3H6']

# b) Set up a data structure to store reactions: We'll use a list of dictionaries to store reaction details: [Note that this aint the optimal way to store Data, we should use the pickle function perhaps]
# Reaction 1: LIGH → RLIGM2A + OH + C3H6
reactions = [
    {
        'stoichiometry': {'LIGH': -1, 'RLIGM2A': 1, 'OH': 1, 'C3H6': 1},
        'orders': {'LIGH': 1},  # First-order w.r.t LIGH
        'A': 1e13,             # Pre-exponential factor (s^-1)
        'Ea': 163254,          # Activation energy (J/mol)
        'phase': 'solid',      # Reactant phase
        'type': 'initiation',  # Reaction type
        'reversible': False    # Irreversible reaction
    }
    # Additional reactions can be added HERE
]

# Define Initial concentrations (mol/L)
initial_concentrations = {
    'LIGH': 1.0,  # Replace with your actual initial concentration
    'RLIGM2A': 0.0,
    'OH': 0.0,
    'C3H6': 0.0
}

# Convert initial concentrations to a list in the order of 'species'
initial_conc_list = [initial_concentrations[specie] for specie in species]
# Time points (seconds)
time_points = np.linspace(0, 1000000, 1000)  # Adjust as needed (depending on reaction time, adjust to see the whole scope)

# Temperature in Kelvin
T = 200 + 273.15  # Convert from °C to K

# # c) Modify the ODE function to handle multiple reactions and species (Generalize Batch_reactor function)
def batch_reactor(concentrations, t, reactions, T):
    conc_dict = dict(zip(species, concentrations))
    dCdt = {specie: 0.0 for specie in species}
    for reaction in reactions:
        A = reaction['A']
        Ea = reaction['Ea']
        k = arrhenius_constant(A, Ea, T)
        rate = k
        for reactant, order in reaction['orders'].items():
            rate *= conc_dict[reactant] ** order
        for specie, stoich_coeff in reaction['stoichiometry'].items():
            dCdt[specie] += stoich_coeff * rate
    # Return the rates as a list in the order of 'species'
    return [dCdt[specie] for specie in species]

# Solve the ODEs
concentration_profiles = odeint(batch_reactor, initial_conc_list, time_points, args=(reactions, T))

# Extract results
results = {specie: concentration_profiles[:, idx] for idx, specie in enumerate(species)}

# Optionally, plot the results
import matplotlib.pyplot as plt

#plt.figure(figsize=(12, 6))
#for specie in species:
#    plt.plot(time_points, results[specie], label=specie)
#plt.xlabel('Time (s)')
#plt.ylabel('Concentration (mol/L)')
#plt.title('Concentration Profiles')
#plt.legend()
#plt.show()

import matplotlib.pyplot as plt
#Alternative way to plot 3 species:

# First plot: LIGH and RLIGM2A
#plt.figure(figsize=(8, 6))
#plt.plot(time_points, results['LIGH'], label='LIGH (Reactant)')
#plt.plot(time_points, results['RLIGM2A'], label='RLIGM2A (Product)')
#plt.xlabel('Time (s)')
#plt.ylabel('Concentration (mol/L)')
#plt.title('Concentration Profiles of LIGH and RLIGM2A')
#plt.legend()
#plt.grid(True)
#plt.show()

# Second plot: LIGH and OH
#plt.figure(figsize=(8, 6))
#plt.plot(time_points, results['LIGH'], label='LIGH (Reactant)')
#plt.plot(time_points, results['OH'], label='OH (Product)')
#plt.xlabel('Time (s)')
#plt.ylabel('Concentration (mol/L)')
#plt.title('Concentration Profiles of LIGH and OH')
#plt.legend()
#plt.grid(True)
#plt.show()

# Third plot: LIGH and C3H6
#plt.figure(figsize=(8, 6))
#plt.plot(time_points, results['LIGH'], label='LIGH (Reactant)')
#plt.plot(time_points, results['C3H6'], label='C3H6 (Product)')
#plt.xlabel('Time (s)')
#plt.ylabel('Concentration (mol/L)')
#plt.title('Concentration Profiles of LIGH and C3H6')
#plt.legend()
#plt.grid(True)
#plt.show()

import matplotlib.pyplot as plt

fig, axs = plt.subplots(3, 1, figsize=(10, 15), sharex=True)

# First subplot: LIGH and RLIGM2A
axs[0].plot(time_points, results['LIGH'], label='LIGH (Reactant)')
axs[0].plot(time_points, results['RLIGM2A'], label='RLIGM2A (Product)')
axs[0].set_ylabel('Concentration (mol/L)')
axs[0].set_title('LIGH and RLIGM2A')
axs[0].legend()
axs[0].grid(True)

# Second subplot: LIGH and OH
axs[1].plot(time_points, results['LIGH'], label='LIGH (Reactant)')
axs[1].plot(time_points, results['OH'], label='OH (Product)')
axs[1].set_ylabel('Concentration (mol/L)')
axs[1].set_title('LIGH and OH')
axs[1].legend()
axs[1].grid(True)

# Third subplot: LIGH and C3H6
axs[2].plot(time_points, results['LIGH'], label='LIGH (Reactant)')
axs[2].plot(time_points, results['C3H6'], label='C3H6 (Product)')
axs[2].set_xlabel('Time (s)')
axs[2].set_ylabel('Concentration (mol/L)')
axs[2].set_title('LIGH and C3H6')
axs[2].legend()
axs[2].grid(True)

plt.tight_layout()
plt.show()




