"""Utility functions (e.g., plotting and logging etc.)
""" 
# modules/utils.py
import matplotlib.pyplot as plt

def plot_concentration_subplots(time_points, results):
    # Define species groups

    # Group 1: Initial species + small gas molecules
    group1_species = ['PLIGH', 'PLIGO', 'PLIGC', 'H2', 'CH4', 'CO', 'CO2', 'H2O']

    # Group 2: Initial species + termination products
    termination_species = ['Char', 'KETDM2', 'KETD', 'ETOH', 'VKETM2', 'VKETDM2',
                           'VKET', 'VKETD', 'RKETM2', 'C3H8O2', 'C3H4O', 'PFET']
    group2_species = ['PLIGH', 'PLIGO', 'PLIGC'] + termination_species

    # Group 3: All other intermediates
    all_species = list(results.keys())
    group3_species = [species for species in all_species if species not in set(group1_species + group2_species)]

    # Increase font size for better readability
    plt.rcParams.update({'font.size': 11})

    ### Plot Group 1 in a separate maximized window
    fig1 = plt.figure(figsize=(16, 9))
    ax1 = fig1.add_subplot(111)
    for species in group1_species:
        if species in results:
            ax1.plot(time_points, results[species], label=species)
        else:
            print(f"Warning: Species '{species}' not found in results.")
    ax1.set_title('Group 1: Initial Species and Small Gas Molecules')
    ax1.set_xlabel('Time (s)')
    ax1.set_ylabel('Concentration [Kg/m3]')
    ax1.legend(loc='best', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    ax1.grid(True)
    ax1.set_ylim(bottom=0)
    plt.tight_layout()
    # Maximize the figure window
    mng = plt.get_current_fig_manager()
    try:
        mng.window.state('zoomed')  # For Windows OS
    except AttributeError:
        try:
            mng.window.showMaximized()  # For Linux/Mac OS
        except AttributeError:
            pass  # Cannot maximize window
    plt.show()

    ### Plot Group 2 in a separate maximized window
    fig2 = plt.figure(figsize=(16, 9))
    ax2 = fig2.add_subplot(111)
    for species in group2_species:
        if species in results:
            ax2.plot(time_points, results[species], label=species)
        else:
            print(f"Warning: Species '{species}' not found in results.")
    ax2.set_title('Group 2: Initial Species and Termination Products')
    ax2.set_xlabel('Time (s)')
    ax2.set_ylabel('Concentration [Kg/m3]')
    ax2.legend(loc='best', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    ax2.grid(True)
    ax2.set_ylim(bottom=0)
    plt.tight_layout()
    # Maximize the figure window
    mng = plt.get_current_fig_manager()
    try:
        mng.window.state('zoomed')
    except AttributeError:
        try:
            mng.window.showMaximized()
        except AttributeError:
            pass
    plt.show()

    ### Plot Group 3 in a separate maximized window
    fig3 = plt.figure(figsize=(16, 9))
    ax3 = fig3.add_subplot(111)
    for species in group3_species:
        if species in results:
            ax3.plot(time_points, results[species], label=species)
        else:
            print(f"Warning: Species '{species}' not found in results.")
    ax3.set_title('Group 3: Intermediates Evolution')
    ax3.set_xlabel('Time (s)')
    ax3.set_ylabel('Concentration [kg/m3]')
    ax3.legend(loc='best', bbox_to_anchor=(1.05, 1), borderaxespad=0.)
    ax3.grid(True)
    ax3.set_ylim(bottom=0)
    plt.tight_layout()
    # Maximize the figure window
    mng = plt.get_current_fig_manager()
    try:
        mng.window.state('zoomed')
    except AttributeError:
        try:
            mng.window.showMaximized()
        except AttributeError:
            pass
    plt.show()

"""
import matplotlib.pyplot as plt

def plot_concentration_subplots(time_points, results):
    # Define species groups

    # Group 1: Initial species + small gas molecules
    group1_species = ['PLIGH', 'PLIGO', 'PLIGC', 'H2', 'CH4', 'CO', 'CO2', 'H2O']

    # Group 2: Initial species + termination products
    termination_species = ['Char', 'KETDM2', 'KETD', 'ETOH', 'VKETM2', 'VKETDM2',
                           'VKET', 'VKETD', 'RKETM2', 'C3H8O2', 'C3H4O', 'PFET']
    group2_species = ['PLIGH', 'PLIGO', 'PLIGC'] + termination_species

    # Group 3: All other intermediates
    all_species = list(results.keys())
    group3_species = [species for species in all_species if species not in set(group1_species + group2_species)]

    # Create subplots
    fig, axes = plt.subplots(3, 1, figsize=(12, 18))

    # Plot Group 1
    for species in group1_species:
        if species in results:
            axes[0].plot(time_points, results[species], label=species)
        else:
            print(f"Warning: Species '{species}' not found in results.")
    axes[0].set_title('Group 1: Initial Species and Small Gas Molecules')
    axes[0].set_xlabel('Time (s)')
    axes[0].set_ylabel('Concentration')
    axes[0].legend()
    axes[0].grid(True)

    # Plot Group 2
    for species in group2_species:
        if species in results:
            axes[1].plot(time_points, results[species], label=species)
        else:
            print(f"Warning: Species '{species}' not found in results.")
    axes[1].set_title('Group 2: Initial Species and Termination Products')
    axes[1].set_xlabel('Time (s)')
    axes[1].set_ylabel('Concentration')
    axes[1].legend()
    axes[1].grid(True)

    # Plot Group 3
    for species in group3_species:
        if species in results:
            axes[2].plot(time_points, results[species], label=species)
        else:
            print(f"Warning: Species '{species}' not found in results.")
    axes[2].set_title('Group 3: Intermediates Evolution')
    axes[2].set_xlabel('Time (s)')
    axes[2].set_ylabel('Concentration')
    axes[2].legend()
    axes[2].grid(True)

    plt.tight_layout()
    plt.show()




    


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

#create a new function in utils.py called plot_concentration_subplots that will plot three subplots in the same window, each showing the concentration evolution of one reactant and one product.
"""