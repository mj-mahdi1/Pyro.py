"""Contains the Reaction class and related functions.
As we add more reactions, simply update the reactions.csv file.
No changes to the code are necessary unless there is need to implement new kinetic models or mechanisms.
"""
# modules/reactions.py
#using the __init__ method to structure data: 

class Reaction:
    def __init__(self, reaction_id, stoichiometry, orders, phase, A, Ea, reversible=False):
        self.reaction_id = reaction_id
        self.stoichiometry = stoichiometry  # Dictionary {species: coefficient}
        self.orders = orders                # Dictionary {species: order}
        self.phase = phase                  # Phase information (e.g., 'solid', 'gas')
        self.A = A                          # Pre-exponential factor
        self.Ea = Ea                        # Activation energy
        self.reversible = reversible

    def rate_constant(self, T):
        from modules.kinetics import arrhenius_constant
        return arrhenius_constant(self.A, self.Ea, T)

    def rate(self, C, species_indices, T):
        k = self.rate_constant(T)
        rate = k
        for specie, order in self.orders.items():
            index = species_indices[specie]
            rate *= C[index] ** order
        return rate
