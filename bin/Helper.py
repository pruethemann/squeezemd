#!/usr/bin/env python

"""
This Helper modules contains multiple function used by multiple other modules.

"""

import subprocess
import os
import yaml
import MDAnalysis as mda
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


def save_file(content, output_file):
    """
    Saves a string (content) to a text file and closes it.
    :param content:
    :param output_file:
    :return:
    """

    with open(output_file, 'w') as file:
        file.write(content)

def execute(command):
    """
    Executes commands in console
    :param command:
    :return:
    """

    output_text = subprocess.check_output(command, shell=True)
    return output_text

def import_yaml(yaml_path: os.path):
    """
    Opens yaml file containing hyper parameters.

    :param yaml_path: File path to yaml
    :return: dictionary with parameters
    """
    try:
        with open(yaml_path, 'r') as stream:
            return yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


def extract_ligand_sequence(pdb_ligand: os.path):
    # Import pdb file with MDAnalysis
    u = mda.Universe(pdb_ligand)

    # Extract ligand at chain I
    ligand = u.select_atoms('segid I')

    # Return sequence
    return str(ligand.residues.sequence().seq)

def save_yaml(d, filepath):
    """
    Save a yml file
    :param d:
    :param filepath:
    :return:
    """
    with open(filepath, 'w') as file:
        documents = yaml.dump(d, file)


def chain2resid(file_csv):
    """
    Amber preparation removes chain ids and starts renumbering residues from 1.
    This function remaps the numbering based on the amber mapping file.
    :param file_csv:
    :return:
    """
    # Find start and end of chain A

    renum = pd.read_csv(file_csv,
                        delim_whitespace=True,
                        names=['resname', 'chainID', 'resid', 'resname amber', 'resid amber'])

    del renum['resname']
    del renum['resname amber']

    chain_min = renum.groupby('chainID').min().rename(columns={'resid amber': 'amber_start', 'resid': 'start'})
    chain_max = renum.groupby('chainID').max().rename(columns={'resid amber': 'amber_end', 'resid': 'end'})

    chains = pd.concat([chain_min, chain_max], axis=1)

    # TODO: Add one resname of every chain for start resid

    return chains

def remap_MDAnalysis(u: mda.Universe, topo):
    """
    Remaps the correct residues Ids from the OpenMM topology to
    a MDAnylsis universe.

    TODO:
    - Dicts not necessary
    - Chain remapping not yet implemented. Adapt afterwards chainID selection

    :param u: MDAnalysis Universe
    :param topo: OpenMM Toplogy
    :return: Mapping tables from chainIDs to original Ids
    """

    chainIds = {}
    for chain_cont, chainID in zip(u.segments, topo.topology.chains()):
        chainIds[int(chain_cont.segid)] = chainID.id

    resIds = {}
    for res_cont, resid in zip(u.residues, topo.topology.residues()):
        #print(res_cont.resid, resid.id, resid.name)
        resIds[int(res_cont.resid)] = int(resid.id)

        resid_sele = u.select_atoms(f"resid {int(res_cont.resid)}")
        resid_sele.residues.resids = int(resid.id)



    return u

def remap_amber(mapping_file, u):
    """
    Amber preparation removes chain ids and starts renumbering residues from 1.
    This function remaps the numbering based on the amber mapping file.
    :param args:
    :param u:
    :return:
    """

    print("2. Import reside information")
    # Get chain information since amber deleted all chain IDs
    chains = chain2resid(mapping_file)

    # IMPORTANT: chain ID A needs to be at the end, otherwise it leads to misnumbering
    chains.sort_index(ascending=False, inplace=True)

    # 1. Assign chain Ids to amber resids
    for chainID, r in chains.iterrows():
        # Assign chain ID to amber resid numbering
        chain = u.select_atoms(f"resid {r.amber_start} to {r.amber_end}")
        chain.atoms.chainIDs = chainID

    # 2. Renumber resids
    for chainID, r in chains.iterrows():

        # Assign chain ID to amber resid numbering
        chain = u.select_atoms(f"chainID {chainID}")

        # Shift resid numbering to old numbering
        shift_factor = r.start - int(r.amber_start)
        chain.residues.resids += shift_factor

    return u
