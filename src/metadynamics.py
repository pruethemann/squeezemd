from openmmplumed import PlumedForce
import os
from openmm.unit import kelvin
import MDAnalysis as mda
import numpy as np


def add_metadynamics_forces_singledistance(metadynamics_params, T:int, system, args):
    """
    in progress
    """
    # Example: add a distance-based collective variable
    atomId1 = metadynamics_params[0]['d1'][0]['atomId1'] + 1
    atomId2 = metadynamics_params[0]['d1'][1]['atomId2'] + 1

    hills_path = os.path.abspath(args.metadynamics_hills)
    colvar_path = os.path.abspath(args.metadynamics_colvar)

    script = f"""
            d1: DISTANCE ATOMS={atomId1},{atomId2}
            METAD ARG=d1 SIGMA=0.1 HEIGHT=0.3 PACE=50 FILE={hills_path}
            PRINT ARG=d1 STRIDE=50 FILE={colvar_path}
            """
    plumed = PlumedForce(script)
    plumed.setTemperature(T*kelvin)
    system.addForce(plumed)
    print("Metadynamics variable added")
    return system


def add_metadynamics_forces_centerofmass(params, system, args):
    # Metadynamics params
    sigma = params['simulation']['metadynamics']['SIGMA']
    height = params['simulation']['metadynamics']['HEIGHT']
    pace = params['simulation']['metadynamics']['PACE']
    stride = params['simulation']['metadynamics']['STRIDE']

    # Get absolute paths for outputs
    hills_path = os.path.abspath(args.metadynamics_hills)
    colvar_path = os.path.abspath(args.metadynamics_colvar)

    # get relevant atom indexes
    id = extract_atom_indices(args.equilibrated)

    script = f"""
            # get residue and chainID information
            MOLINFO STRUCTURE={args.equilibrated}

            # Define two groups (ligand:Entity0 and receptor:Entity1)
            WHOLEMOLECULES ENTITY0={id['lig_min']}-{id['lig_max']} ENTITY1={id['rec_min']}-{id['rec_max']}

            # Group heavy atoms for contact
            grp_lig: GROUP ATOMS={id['lig_min']}-{id['lig_max']}
            grp_rec: GROUP ATOMS={id['rec_min']}-{id['rec_max']}

            # Define center of mass of the two partners
            lig: COM ATOMS=grp_lig
            rec: COM ATOMS=grp_rec

            # Distance between the two COMs (in nm)
            d1: DISTANCE ATOMS=lig,rec

            METAD ARG=d1 SIGMA={sigma} HEIGHT={height} PACE={pace} FILE={hills_path}
            PRINT ARG=d1 STRIDE={stride} FILE={colvar_path}
            """

    plumed = PlumedForce(script)
    plumed.setTemperature(T*kelvin)
    system.addForce(plumed)
    print("Metadynamics variable added")
    return system

def extract_atom_indices(pdf_file: os.path, cutoff = 5.0):

    u = mda.Universe(pdf_file)

    # get all heavy atoms of lig and rec
    lig_heavy = u.select_atoms("chainID A and not name H*")
    rec_heavy = u.select_atoms("(chainID B or chainID C) and not name H*")

    # get all heavy atoms of lig and rec
    lig = u.select_atoms("chainID A")
    rec = u.select_atoms("chainID B or chainID C")

    # Compute distance matrix between all atoms of the two chains
    dist = mda.lib.distances.distance_array(lig_heavy.positions, rec_heavy.positions)

    # Boolean masks of interface atoms
    lig_interface_mask = np.any(dist < cutoff, axis=1)
    rec_interface_mask = np.any(dist < cutoff, axis=0)

    # Interface atoms selections
    lig_interface = lig_heavy[lig_interface_mask]
    rec_interface = rec_heavy[rec_interface_mask]

    # Print in PLUMED-friendly format. Add +1 because plumed starts at atom id 1 and not 0
    lig_plumed = ",".join(map(str, lig.indices + 1))
    rec_plumed = ",".join(map(str, rec.indices + 1))
    lig_interface_plumed = ",".join(map(str, lig_interface.indices + 1))
    rec_interface_plumed = ",".join(map(str, rec_interface.indices + 1))

    atom_indices = {'lig_index': lig_plumed,
                    'rec_index': rec_plumed,
                    'lig_interface_index':lig_interface_plumed,
                    'rec_interface_index':rec_interface_plumed,
                    'lig_min':lig.indices.min() + 1,
                    'lig_max':lig.indices.max() + 1,
                    'rec_min':rec.indices.min() + 1,
                    'rec_max':rec.indices.max() + 1,
    }

    return atom_indices