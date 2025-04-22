import os
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

from tddft_data_parser import TddftDataParser
from nto_data_parser import NtoDataParser


transition_metal_identifiers = [
    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
    'La', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg'
]

def is_metal(atom):

    """Checks if a given atom is a metal.

    Arguments:
        atom (dict): The atom.

    Returns:
        bool: The flag indicating whether the atom is metal.
    """

    if atom['atom_element'] in transition_metal_identifiers:
        return True

    return False

def build_transition_nature_string(occupied_origin: str, virtual_origin: str):

    """Builds the correct string denoting the given transition nature.

    Arguments:
        occupied_origin (str): The occupied origin type.
        virtual_origin (str): The virtual origin type.

    Returns:
        str: The transition nature string.
    """

    if occupied_origin == 'M' and virtual_origin == 'M':
        return 'ddT'

    if occupied_origin == 'M' and virtual_origin == 'L':
        return 'MLCT'

    if occupied_origin == 'L' and virtual_origin == 'M':
        return 'LMCT'

    if occupied_origin == 'L' and virtual_origin == 'L':
        return 'LLCT'

    return

def get_nto_origin_metal_ligand_ratios(nto_data: list[dict[dict]]):

    """Gets the relative NTO metal and ligand orbital contributions.

    Arguments:
        nto_data (list[dict[dict]]): The NTO data.

    Returns:
        float: The relative metal contribution.
        float: The relative ligand contribution.
    """

    metal_occupations, ligand_occupations = get_metal_ligand_occupations(nto_data)

    metal_ratio = sum(metal_occupations) / (sum(metal_occupations) + sum(ligand_occupations))
    ligand_ratio = sum(ligand_occupations) / (sum(metal_occupations) + sum(ligand_occupations))

    return metal_ratio, ligand_ratio

def get_nto_origin(nto_data: list[dict[dict]], metal_ratio_threshold: float=0.5):

    """Gets the NTO origin.

    Arguments:
        nto_data (list[dict[dict]]): The NTO data.
        metal_ratio_threshold (float): The threhold to determine metal origin.

    Returns:
        str: The metal origin.
    """

    metal_ratio, ligand_ratio = get_nto_origin_metal_ligand_ratios(nto_data)

    if metal_ratio > metal_ratio_threshold:
        return 'M'

    return 'L'

def get_metal_ligand_occupations(nto_data: list[dict[dict]], normalize_by_atoms: bool=False):

    """Gets the absolute NTO metal and ligand orbital occupations.

    Arguments:
        ntos (list[dict[dict]]): The NTO data.
        normalize_by_atoms (bool): Flag to denote whether to normalize by the number of atoms.

    Returns:
        list[float]: The absolute metal occupations.
        list[float]: The absolute ligand occupations.
    """

    metal_occupations = None
    ligand_occupations = [0, 0, 0, 0, 0]

    for atom in nto_data:

        average_absolute_occupations = get_average_absolute_occupations(atom['ntos'])

        if is_metal(atom):
            metal_occupations = average_absolute_occupations
        else:
            for i, _ in enumerate(average_absolute_occupations):
                ligand_occupations[i] += average_absolute_occupations[i]

    # normalize by number of ligand atoms
    if normalize_by_atoms:
        for i, _ in enumerate(ligand_occupations):
            ligand_occupations[i] = ligand_occupations[i] / ( len(nto_data) - 1)

    return metal_occupations, ligand_occupations

def get_average_absolute_occupations(ntos: dict[dict]):

    """Gets the average absolute NTO orbital occupations.

    Arguments:
        ntos (dict[dict]): The NTO data.

    Returns:
        list[float]: The absolute averages.
    """

    absolute_averages = []

    #for _ in ntos.keys():
    for _ in [max(ntos.keys())]:
        absolute_averages.append(_ * get_average_absolute_occupation(ntos[_]))

    return absolute_averages

def get_average_absolute_occupation(nto: dict, normalize_by_orbitals: bool=False):

    """Gets the average absolute orbital occupation of a single NTO.

    Arguments:
        nto (dict): The NTO data.
        normalize_by_orbitals (bool): Flag to denote whether to normalize by the number of orbitals.

    Returns:
        float: The absolute sum.
    """

    absolute_sum = 0
    for _ in nto.keys():
        absolute_sum += np.abs(nto[_])

    # normalize by number of orbitals
    if normalize_by_orbitals:
        return absolute_sum / len(nto.keys())

    return absolute_sum


# read files from input directory
data_dir = sys.argv[1].strip()
out_files = [_ for _ in os.listdir(data_dir) if '.out' in _]
out_files_tddft = sorted([_ for _ in out_files if not 'vis.out' in _])

result_dicts = []
for _ in tqdm(out_files_tddft):

    tddft_result_dict = TddftDataParser(data_dir + _).parse()

    if tddft_result_dict['has_failed']:
        continue

    has_negative_f = False
    for i in range(1, 31, 1):
        if tddft_result_dict['f_' + str(i)] < 0:
            has_negative_f = True
            break
    if has_negative_f:
        continue

    if tddft_result_dict['lambda_max_vis'] is not None:
        nto_result_dict = NtoDataParser(data_dir + _.replace('.out', '-vis.out')).parse()

        if nto_result_dict['has_failed']:
            continue

        #print(nto_result_dict['occupied_nto'][0]['ntos'])

        occupied_origin = get_nto_origin(nto_result_dict['occupied_nto'], 0.2)
        virtual_origin = get_nto_origin(nto_result_dict['virtual_nto'], 0.2)

        occupied_origin_ratio = get_nto_origin_metal_ligand_ratios(nto_result_dict['occupied_nto'])
        virtual_origin_ratio = get_nto_origin_metal_ligand_ratios(nto_result_dict['virtual_nto'])

        tddft_result_dict['transition_nature_vis'] = build_transition_nature_string(occupied_origin, virtual_origin)

        tddft_result_dict['M_contribution_occupied'] = occupied_origin_ratio[0]
        tddft_result_dict['L_contribution_occupied'] = occupied_origin_ratio[1]
        tddft_result_dict['M_contribution_virtual'] = virtual_origin_ratio[0]
        tddft_result_dict['L_contribution_virtual'] = virtual_origin_ratio[1]

    else:
        tddft_result_dict['transition_nature_vis'] = None

    result_dicts.append(tddft_result_dict)

df = pd.DataFrame(result_dicts)
df.to_csv('df_' + data_dir.replace('/', '') + '.csv', index=False)
