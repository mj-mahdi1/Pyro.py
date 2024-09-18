"""Functions to load and process data (to read data from reactions.csv and species.csv)
"""
# modules/data_handler.py

import pandas as pd
from modules.reactions import Reaction
from modules.species import Species
import json

def load_species_data(file_path):
    species_df = pd.read_csv(file_path)
    species_list = species_df['SpeciesID'].tolist()
    species_indices = {specie: idx for idx, specie in enumerate(species_list)}
    species_objects = []
    for idx, row in species_df.iterrows():
        species = Species(
            species_id=row['SpeciesID'],
            molecular_weight=row['MolecularWeight'],
            phase=row['Phase'],
            initial_concentration=row['InitialConcentration']
        )
        species_objects.append(species)
    return species_objects, species_indices

def load_reaction_data(file_path):
    reactions_df = pd.read_csv(file_path)
    reactions_list = []
    for idx, row in reactions_df.iterrows():
        reaction = Reaction(
            reaction_id=row['ReactionID'],
            stoichiometry=json.loads(row['Stoichiometry']),
            orders=json.loads(row['Orders']),
            phase=row['Phase'],             # Added phase here
            A=row['A'],
            Ea=row['Ea'],
            reversible=row['Reversible']
        )
        reactions_list.append(reaction)
    return reactions_list
