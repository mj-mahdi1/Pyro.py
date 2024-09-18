
#Purpose: Compute how the concentration of a reactant changes over time in a batch reactor, considering temperature-dependent reaction rates.
#Parameters calculated :concentrations of reactants at different time points.

import numpy as np 
import matplotlib.pyplot as plt

from scipy.integrate import odeint
 
# Define constants

R = 8.314  # J/(mol*K), universal gas constant
 
# Arrhenius equation to calculate the rate constant

def arrhenius_constant(A, Ea, T):

    return A * np.exp(-Ea / (R * T))
 
# Rate law for the reaction: Calculates the reaction rate based on the rate constant and the concentrations exposant to their respective reaction orders : Rate = k. [conc]^order

def reaction_rate(k, concentrations, orders):

    rate = k

    for conc, order in zip(concentrations, orders): #concentrations: List of reactant concentrations, orders: List of reaction orders for each reactant.

        rate *= conc ** order

    return rate
 
# Differential equations for batch reactor

def batch_reactor(concentrations, t, A, Ea, T, orders): #Defines the differential equations representing the rate of change of concentrations (dCdt) over time.

    k = arrhenius_constant(A, Ea, T) #Arrhenus parameters

    rate = reaction_rate(k, concentrations, orders) #concentrations: Current concentrations of reactants, t: Time variable (required by odeint),

    dCdt = [-rate]  # Assuming first-order reaction for simplicity

    return dCdt

#Set and define values of Initial Conditions and Parameters 
# Initial concentration(s) of reactant(s)

initial_concentration = [1.0]  # mol/L, for examplet
 
# Time points (in seconds)

time_points = np.linspace(0, 1000, 100)  # for example
 
# Arrhenius constants

A = 1e12  # Pre-exponential factor in s^-1

Ea = 80000  # Activation energy in J/mol

T = 298  # Temperature in K   
 
# Reaction orders for each reactant

orders = [1]  # Assuming first-order reaction
 
# Solve the ODE
#Uses odeint Solver to integrate the differential equations and compute concentrations over the specified time points.
#The solver solves the differential equation dCdt = -k. C^order (first order form of ODE)
concentrations = odeint(batch_reactor, initial_concentration, time_points, args=(A, Ea, T, orders))

plt.plot(time_points, concentrations)
plt.xlabel('Time (s)')
plt.ylabel('Concentration (mol/L)')
plt.title('Concentration vs. Time')
plt.show()