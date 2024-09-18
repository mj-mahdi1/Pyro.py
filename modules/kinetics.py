"""Functions for kinetic calculations.
"""
# modules/kinetics.py

import numpy as np

# Universal gas constant
R = 8.314  # J/(mol*K)

def arrhenius_constant(A, Ea, T):
    return A * np.exp(-Ea / (R * T))
