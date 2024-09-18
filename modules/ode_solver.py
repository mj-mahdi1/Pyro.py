"""Functions to set up and solve the ODE system.
"""
# modules/ode_solver.py

import numpy as np
from scipy.integrate import solve_ivp

def odes(t, C, reactions_list, species_indices, T):
    dCdt = np.zeros_like(C)
    for reaction in reactions_list:
        rate = reaction.rate(C, species_indices, T)
        for specie, coeff in reaction.stoichiometry.items():
            index = species_indices[specie]
            dCdt[index] += coeff * rate
    return dCdt

def solve_odes(t_span, C0, reactions_list, species_indices, T, t_eval=None):
    result = solve_ivp(
        fun=lambda t, C: odes(t, C, reactions_list, species_indices, T),
        t_span=t_span,
        y0=C0,
        method='BDF',  # Suitable for stiff systems
        t_eval=t_eval
    )
    return result
