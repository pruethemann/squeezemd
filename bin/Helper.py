#!/usr/bin/env python

"""
This Helper modules contains multiple function used by multiple other modules.

"""

import subprocess
import os
import MDAnalysis as mda
import pandas as pd



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
    import yaml
    try:
        with open(yaml_path, 'r') as stream:
            return yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)


def extract_ligand_sequence(pdb_ligand: os.path):
    # Import pdb file with MDAnalysis
    u = mda.Universe(pdb_ligand)

    # Extract ligand at chain A
    ligand = u.select_atoms('chainID A')

    # Return sequence
    return str(ligand.residues.sequence().seq)

def save_yaml(d, filepath):
    """
    Save a yml file
    :param d:
    :param filepath:
    :return:
    """
    import yaml
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

def is_numeric(character):
    """
    This function checks if a given character is numeric.

    :param character: A single character (string) to check.
    :return: True if the character is numeric, False otherwise.
    """
    if len(character) != 1:
        raise ValueError("Input must be a single character.")
    return character.isdigit()

def remap_MDAnalysis_V1(u: mda.Universe, topo):
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

    # Currently resets the segment ID to the original chainID
    for chain_cont, chainID in zip(u.segments, topo.topology.chains()):
        if chainID.id == '1': continue
        if chainID.id == '2': continue
        if chainID.id == '3': continue
        if chainID.id == '4': continue
        if chainID.id == '5': continue
        selected_segid = u.select_atoms(f"segid {chain_cont.segid}")
        selected_segid.segments.segids = chainID.id

    for res_cont, resid in zip(u.residues, topo.topology.residues()):

        if is_numeric(resid.chain.id):
            continue
        resid_sele = u.select_atoms(f"resid {int(res_cont.resid)}")
        resid_sele.residues.resids = int(resid.id)

    return u

import openmm.app as app

def remap_MDAnalysis(u: mda.Universe, topo: app.PDBxFile):
    """
    Remaps the correct residue and chain IDs from the OpenMM PDBxFile topology
    to an MDAnalysis universe.

    :param u: MDAnalysis Universe
    :param topo: openmm.app.PDBxFile object
    :return: updated MDAnalysis Universe
    """
    chains = list(topo.topology.chains())
    residues = list(topo.topology.residues())

    if len(u.segments) != len(chains):
        raise ValueError("Mismatch in number of segments and chains")

    for mda_seg, omm_chain in zip(u.segments, chains):
        mda_seg.segid = omm_chain.id  # Safe: assigns chainID

    if len(u.residues) != len(residues):
        raise ValueError("Mismatch in number of residues between MDAnalysis and OpenMM topology")

    for mda_res, omm_res in zip(u.residues, residues):
        # Optional: Only remap if different
        mda_res.resid = int(omm_res.id)
        mda_res.resname = omm_res.name

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
