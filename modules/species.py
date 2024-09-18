"""Contains the Species class and related functions.
"""
# modules/species.py
#using the __init__ method:

class Species:
    def __init__(self, species_id, molecular_weight, phase, initial_concentration):
        self.species_id = species_id
        self.molecular_weight = molecular_weight
        self.phase = phase
        self.initial_concentration = initial_concentration
