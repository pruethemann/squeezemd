# metadynamics_auxillary.py

from openmmplumed import PlumedForce
import os
from openmm.unit import kelvin
import MDAnalysis as mda


def add_metadynamics_forces_centerofmass_contacts(params, system, args, T=300):
    """
    Plain metadynamics with 2 CVs:
      CV1: COM distance between ligand and receptor (nm)
      CV2: Interface contacts (dimensionless coordination number)

    Notes:
    - This is meant to prevent "cheap" dissociation pathways where only a floppy terminus peels off.
    - For performance, contacts are computed on BACKBONE atoms by default (editable below).
    """
    print("temperature", T)

    meta = params["simulation"]["metadynamics"]
    sigma_com = meta["SIGMA_COM"]                 # nm
    sigma_contacts = meta["SIGMA_CONTACTS"]       # dimensionless (contacts)
    height = meta["HEIGHT"]
    pace = meta["PACE"]
    stride = meta["STRIDE"]

    # Smooth contact switching function parameters
    # R_0 in nm (0.45 nm ~ 4.5 Å is a common start); NN controls steepness
    r0 = meta.get("CONTACT_R0", 0.45)
    nn = meta.get("CONTACT_NN", 6)

    # Output files
    hills_path = os.path.abspath(args.metadynamics_hills)
    colvar_path = os.path.abspath(args.metadynamics_colvar)

    # Atom indices (you already have this helper)
    idx = extract_atom_indices(args.equilibrated)

    # ---- Choose groups for COM and contacts ----
    # COM: keep your backbone-based COM (fast)
    # Contacts: by default also backbone-only to reduce cost (can be made heavier below)
    #
    # If you want *more sensitive* contacts, switch grp_lig_cnt/grp_rec_cnt to all atoms
    # or heavy atoms (recommended only if performance is ok).
    #
    # IMPORTANT: MDAnalysis "backbone" includes N,CA,C,O. If your N-terminus is the issue,
    # contacts CV (c1) will force interface breakage rather than just peeling.
    script = f"""
            # Optional but helpful for residue/chain info (not strictly required for GROUP-based CVs)
            MOLINFO STRUCTURE={args.equilibrated}

            # Keep molecules whole across PBC
            WHOLEMOLECULES ENTITY0={idx['lig_min']}-{idx['lig_max']} ENTITY1={idx['rec_min']}-{idx['rec_max']}

            # --- Groups for COM (backbone only, already precomputed as explicit atom lists) ---
            grp_lig_com: GROUP ATOMS={idx['lig_backbone']}
            grp_rec_com: GROUP ATOMS={idx['rec_backbone']}

            lig: COM ATOMS=grp_lig_com
            rec: COM ATOMS=grp_rec_com
            d1: DISTANCE ATOMS=lig,rec NOPBC

            # --- Groups for CONTACTS ---
            # Default: backbone-only (fast, less noisy)
            grp_lig_cnt: GROUP ATOMS={idx['lig_backbone']}
            grp_rec_cnt: GROUP ATOMS={idx['rec_backbone']}

            # CV2: Smooth coordination number (interface contacts)
            c1: COORDINATION GROUPA=grp_lig_cnt GROUPB=grp_rec_cnt R_0={r0} NN={nn} MM=0

            # Plain metadynamics bias on both CVs
            METAD ARG=d1,c1 SIGMA={sigma_com},{sigma_contacts} HEIGHT={height} PACE={pace} FILE={hills_path}

            # Output CVs
            PRINT ARG=d1,c1 STRIDE={stride} FILE={colvar_path}
            """

    plumed = PlumedForce(script)
    plumed.setTemperature(T * kelvin)
    system.addForce(plumed)
    print("Metadynamics variables added: d1 (COM distance) + c1 (contacts)")
    return system


def extract_atom_indices(pdf_file: os.path):
    """Extract ligand/receptor atom indices and backbone subsets for PLUMED."""
    u = mda.Universe(pdf_file)

    lig = u.select_atoms("chainID A")
    rec = u.select_atoms("chainID B or chainID C")

    lig_backbone = u.select_atoms("(chainID A) and backbone")
    rec_backbone = u.select_atoms("((chainID B) or (chainID C)) and backbone")

    # PLUMED atom indexing starts at 1
    lig_plumed_backbone = ",".join(map(str, lig_backbone.indices + 1))
    rec_plumed_backbone = ",".join(map(str, rec_backbone.indices + 1))

    return {
        "lig_backbone": lig_plumed_backbone,
        "rec_backbone": rec_plumed_backbone,
        "lig_min": int(lig.indices.min() + 1),
        "lig_max": int(lig.indices.max() + 1),
        "rec_min": int(rec.indices.min() + 1),
        "rec_max": int(rec.indices.max() + 1),
    }
