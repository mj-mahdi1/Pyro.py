"""Utility functions (e.g., plotting and logging etc.)
""" 
# modules/utils.py

import matplotlib.pyplot as plt

def plot_concentration_profiles(time_points, results, species):
    plt.figure(figsize=(10, 6))
    for specie in species:
        plt.plot(time_points, results[specie], label=specie)
    plt.xlabel('Time (s)')
    plt.ylabel('Concentration (mol/L)')
    plt.title('Concentration Profiles')
    plt.legend()
    plt.grid(True)
    plt.show()

"""create a new function in utils.py called plot_concentration_subplots that will plot three subplots in the same window, each showing the concentration evolution of one reactant and one product."""