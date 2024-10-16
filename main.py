"""This script will be the entry point of the program.
It will import functions and classes from the modules to run the simulation.
Next Steps

IMPORTANT NOTICE: a) In solid-state reactions, the surface area (which might be related to particle size or molar mass)
can also affect the reaction rate. For instance, reactions involving lignin, cellulose, or 
hemicellulose (as in pyrolysis) can be surface-limited rather than concentration-limited.
b) In pyrolysis, different species (e.g., solid and gas phase) interact, and the molar mass helps
determine the mass transfer effects between solid and gas phases, which are relevant for later
stages of model development when we introduce mass and heat transfer effects.
c) its important to implement molar mass to validate the model with mass conservation law (plotteted) also 
in calculating the enthalpy of reaction, the molar mass helps in scaling the heat transfer rate from a molar heat value to a mass-based value showing T evolution

***Rq about Rate constants and Scalability when dealing with very large nb of rxns:
Numerical Instability: High rate constants can cause the solver to take extremely small time steps to maintain stability, leading to long computation times or even solver failures.
Scaling Issues: The large disparity in reaction rates can make it difficult to accurately capture the dynamics of all species in the system.
***Rq about negative concentrations appearing: 
Negative concentrations in your simulation can arise from several factors:
Numerical Instability: Due to stiff equations or inappropriate solver settings.
Incorrect Reaction Stoichiometry: Imbalanced reactions can cause species to appear or disappear unexpectedly.
Improper ODE Implementation: If the species interactions aren't correctly coded, concentrations can go negative. (reaction balance in terms of atoms/stocheo coeff)

- Add more reactions to reactions.csv

2. Incorporate temperature profiles by modifying the T variable in main.py or within the odes function.
3. Consider phase changes and equilibrium constants in future iterations.
4. Error Handling and Validation:

5. Implement error checking in your data loading functions to handle missing data or incorrect formats.
Documentation:
Update the README.md file to document how to run your simulation, the structure of your data files, and any other relevant information

In the end we need to involve temperature, pressure evolution into the model, after adding all the reactions and species, also we should be able to implement mass and heat transfer
equations in parallel so it should become a multilayer model. maybe in the end also implement reactor geometry and phase transfer dynamics (idk how yet) 
and also sensitivity analysis : this final step could be done by a specific solver, also an ai autoevalutation and learning ability should be incorporated
this means: 

1.Introduce Temperature as a Dynamic Variable:
-Modify your ODE system to include temperature as a variable that evolves over time.
-Couple the energy balance with the mass balance equations.

2.Set Up the Energy Balance Equation:
-Formulate an energy balance for the reactor to account for heat generation/consumption by reactions, heat exchange with the surroundings, and any external heating/cooling.

3.Modify the Reaction Rate Calculations:
-Ensure that reaction rates are updated based on the changing temperature at each time step.

4.Update the ODE Solver:
-Adjust your ODE solver to handle the expanded system that now includes temperature AND Sensitivity analysis and correcting along the way.

5.Adjust Initial Conditions and Parameters:
-Include initial temperature and thermal properties such as heat capacities and densities.
Prepare for Future Enhancements:

Design the model structure to facilitate the addition of pressure evolution, mass and heat transfer, reactor geometry, phase dynamics, and sensitivity analysis.
"""

# main.py

import numpy as np
from modules.data_handler import load_species_data, load_reaction_data
from modules.ode_solver import solve_odes
#from modules.utils import plot_concentration_profiles
import matplotlib.pyplot as plt 
from modules.utils import plot_concentration_subplots
#119, 140,161, 182, 203, 224, 245, 266, 308, 329, 350, 371, 392,   3 au lieu de 15
#Rphenol C3h4o2 ph2 Pc2h2 Char???
#133 154 175 196 217 238 259 280  Rphenol: 50 au lieu de 150
#20 : c3h4o2 133 au lieu de 233
#Rphenox : 84 85 : 104 au lieu de 04

def main():
    # Load species data
    species_file = 'data/species.csv'
    species_objects, species_indices = load_species_data(species_file)
    species_list = [species.species_id for species in species_objects]
    
    # Load reaction data
    reactions_file = 'data/reactions.csv'
    reactions_list = load_reaction_data(reactions_file)
    
    # Initial concentrations
    C0 = np.array([species.initial_concentration for species in species_objects])
    
    # Time span and evaluation points
    t_start = 0
    t_end = 180000  #100000# seconds # Increase the end time to see the full transformation
    num_points = 100 #100  # Number of time points
    t_span = (t_start, t_end)
    t_eval = np.linspace(t_start, t_end, num_points)
    
    # Temperature (constant for now)
    T = (230 + 273.15)  # Convert from °C to K
    
    # Solve the ODEs
    solution = solve_odes(t_span, C0, reactions_list, species_indices, T, t_eval)
    
    # Extract results
    concentration_profiles = solution.y
    time_points = solution.t
    
    # Create a dictionary of results
    results = {}
    for idx, specie in enumerate(species_list):
        results[specie] = concentration_profiles[idx]
    
    # Plot the results
    #plot_concentration_profiles(time_points, results, species_list)
    plot_concentration_subplots(time_points, results)
def plot_concentration_profiles(time_points, results, species_list):
    # Check if the data is empty
    if len(time_points) == 0 or not results:
        print("Error: No data available to plot")
        return
if __name__ == '__main__':
    main()




"""
    # Create the plot
    plt.figure(figsize=(10, 6))
    for specie in species_list:
        # Ensure data exists for each specie
        if specie in results and len(results[specie]) > 0:
            plt.plot(time_points, results[specie], label=specie)
        else:
            print(f"Warning: No data available for specie {specie}")

    # Add labels and legend
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration')
    plt.title('Concentration Profiles Over Time')
    plt.legend()
    plt.grid(True)

    # Show the plot
    plt.show()




    # Create the plot
    plt.figure(figsize=(10, 6))
    for specie in species_list:
        # Ensure data exists for each specie
        if specie in results and len(results[specie]) > 0:
            plt.plot(time_points, results[specie], label=specie)
        else:
            print(f"Warning: No data available for specie {specie}")

    # Add labels and legend
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration')
    plt.title('Concentration Profiles Over Time')
    plt.legend()
    plt.grid(True)

    # Show the plot
    plt.show()
"""

