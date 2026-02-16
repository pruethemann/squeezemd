from openmmplumed import PlumedForce
import os
from openmm.unit import kelvin
import MDAnalysis as mda

def add_metadynamics_forces_centerofmass(params, system, args, T=300):
    """
    Currently active version
    Only concider C alpha to reduce computational intensity
    """
    print("temperature", T)
    # Metadynamics params
    sigma = params['simulation']['metadynamics']['SIGMA_COM']
    height = params['simulation']['metadynamics']['HEIGHT']
    pace = params['simulation']['metadynamics']['PACE']
    stride = params['simulation']['metadynamics']['STRIDE']

    # Get absolute paths for outputs
    hills_path = os.path.abspath(args.metadynamics_hills)
    colvar_path = os.path.abspath(args.metadynamics_colvar)

    # Get relevant atom indices (Cα only for COM CV)
    idx = extract_atom_indices(args.equilibrated)

    script = f"""
            # Define two groups (ligand:Entity0 and receptor:Entity1)
            WHOLEMOLECULES ENTITY0={idx['lig_min']}-{idx['lig_max']} \
                        ENTITY1={idx['rec_min']}-{idx['rec_max']}

            # Cα-only groups (explicit indices)
            grp_lig: GROUP ATOMS={idx['lig_backbone']}
            grp_rec: GROUP ATOMS={idx['rec_backbone']}

            # Define center of mass of the two partners
            lig: COM ATOMS=grp_lig
            rec: COM ATOMS=grp_rec

            # Distance between the two COMs (in nm) PBC not required because already handled in OpenMM
            d1: DISTANCE ATOMS=lig,rec NOPBC

            METAD ARG=d1 SIGMA={sigma} HEIGHT={height} PACE={pace} FILE={hills_path}
            PRINT ARG=d1 STRIDE={stride} FILE={colvar_path}
            """

    plumed = PlumedForce(script)
    plumed.setTemperature(T*kelvin)
    system.addForce(plumed)
    print("Metadynamics variable added")
    return system

def extract_atom_indices(pdf_file: os.path):
    """Extract ligand/receptor atom indices and Cα subsets for PLUMED."""

    u = mda.Universe(pdf_file)

    # Get all atoms
    lig = u.select_atoms("chainID A") 
    rec = u.select_atoms("chainID B or chainID C")

    # get only C alphas to reduce computational cost of COM calculation
    lig_backbone = u.select_atoms("(chainID A) and backbone")
    rec_backbone = u.select_atoms("((chainID B) or (chainID C)) and backbone")

    # Print in PLUMED-friendly format. Add +1 because plumed starts at atom id 1 and not 0
    lig_plumed_backbone = ",".join(map(str, lig_backbone.indices + 1))
    rec_plumed_backbone = ",".join(map(str, rec_backbone.indices + 1))

    atom_indices = {'lig_backbone': lig_plumed_backbone,
                    'rec_backbone': rec_plumed_backbone,
                    'lig_min':lig.indices.min() + 1,
                    'lig_max':lig.indices.max() + 1,
                    'rec_min':rec.indices.min() + 1,
                    'rec_max':rec.indices.max() + 1,
    }

    return atom_indices


def add_metadynamics_forces_centerofmass_contacts(params, system, args, T=300):
    """
    Collective variable 1: Center of mass between proteins
    Collective variable 2: Contacts
    """
    print("temperature", T)

    # Metadynamics params
    meta = params['simulation']['metadynamics']
    sigma_com = meta['SIGMA_COM']            # kept for COM CV (nm)
    sigma_contacts = meta['SIGMA_CONTACTS']
    height = meta['HEIGHT']
    pace = meta['PACE']
    stride = meta['STRIDE']

    # Optional: contact-switch parameters
    # R_0 is in nm; NN controls sharpness (larger => sharper)
    r0 = meta.get('CONTACT_R0', 0.45)     # ~4.5 Å is a common start for heavy-atom contacts
    nn = meta.get('CONTACT_NN', 6)

    # Get absolute paths for outputs
    hills_path = os.path.abspath(args.metadynamics_hills)
    colvar_path = os.path.abspath(args.metadynamics_colvar)

    # Get relevant atom indices
    id = extract_atom_indices(args.equilibrated)

    script = f"""
            # get residue and chainID information
            MOLINFO STRUCTURE={args.equilibrated}

            # Define two groups (ligand:Entity0 and receptor:Entity1)
            WHOLEMOLECULES ENTITY0={id['lig_min']}-{id['lig_max']} ENTITY1={id['rec_min']}-{id['rec_max']}

            # Group atoms (adjust to heavy atoms if your index ranges include hydrogens)
            grp_lig: GROUP ATOMS={id['lig_min']}-{id['lig_max']}
            grp_rec: GROUP ATOMS={id['rec_min']}-{id['rec_max']}

            # Define center of mass of the two partners
            lig: COM ATOMS=grp_lig
            rec: COM ATOMS=grp_rec

            # CV1: Distance between the two COMs (in nm)
            d1: DISTANCE ATOMS=lig,rec

            # CV2: Interface contacts as coordination number (dimensionless)
            # This counts (smoothly) how many ligand atoms are within ~R_0 of receptor atoms.
            # If you want "native contacts only", we can instead use CONTACTMAP with a reference.
            c1: COORDINATION GROUPA=grp_lig GROUPB=grp_rec R_0={r0} NN={nn} MM=0

            # Bias both CVs
            METAD ARG=d1,c1 SIGMA={sigma_com},{sigma_contacts} HEIGHT={height} PACE={pace} FILE={hills_path}

            # Print both CVs
            PRINT ARG=d1,c1 STRIDE={stride} FILE={colvar_path}
            """

    plumed = PlumedForce(script)
    plumed.setTemperature(T * kelvin)
    system.addForce(plumed)
    print("Metadynamics variables added: COM distance (d1) + interface contacts (c1)")
    return system

