"""This script will be the entry point of the program.
It will import functions and classes from the modules to run the simulation.
Next Steps
1. Implement Additional Reactions:

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
from modules.utils import plot_concentration_profiles

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
    t_end = 1000000  # seconds # Increase the end time to see the full transformation
    num_points = 1000  # Number of time points
    t_span = (t_start, t_end)
    t_eval = np.linspace(t_start, t_end, 1000)
    
    # Temperature (constant for now)
    T = 200 + 273.15  # Convert from °C to K
    
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
    plot_concentration_profiles(time_points, results, species_list)

if __name__ == '__main__':
    main()
