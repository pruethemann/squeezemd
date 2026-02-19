# metadynamics_auxillary.py

"""
Import observation

lig = u.select_atoms("(chainID A) and name CA")
print(lig)

This will print the atom numbers which pdb numbering. Taking this number for pymol gives the CA (1 based)

lig.indices
This will however give numbers which are -1 (0 based)

lig.ids
Those are 1 based: used for plumed and pymol

"""

from openmmplumed import PlumedForce
import os
from openmm.unit import kelvin
import MDAnalysis as mda
import numpy as np


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
    sigma_com = float(meta.get("SIGMA_COM", meta.get("SIGMA", 0.2)))
    sigma_contacts = float(meta.get("SIGMA_CONTACTS", 8.0))
    height = meta["HEIGHT"]
    pace = meta["PACE"]
    stride = meta["STRIDE"]

    # Smooth contact switching function parameters
    # R_0 in nm (0.45 nm ~ 4.5 Å is a common start); NN controls steepness
    r0 = meta.get("CONTACT_R0", 0.45)
    nn = meta.get("CONTACT_NN", 6)
    contact_atom_mode = meta.get("CONTACT_ATOM_MODE", "ca")
    contact_interface_cutoff = meta.get("CONTACT_INTERFACE_CUTOFF", 0.8)
    contact_max_atoms = int(meta.get("CONTACT_MAX_ATOMS_PER_PARTNER", 120))

    # Output files
    hills_path = os.path.abspath(args.metadynamics_hills)
    colvar_path = os.path.abspath(args.metadynamics_colvar)

    # Atom indices (you already have this helper)
    idx = extract_atom_indices(
        args.equilibrated,
        contact_atom_mode=contact_atom_mode,
        contact_interface_cutoff_nm=contact_interface_cutoff,
        contact_max_atoms_per_partner=contact_max_atoms,
    )

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
            # Interface-focused subsets to keep contact CV cheap for protein-protein systems.
            grp_lig_cnt: GROUP ATOMS={idx['lig_contacts']}
            grp_rec_cnt: GROUP ATOMS={idx['rec_contacts']}

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
    print(
        "Metadynamics contact groups:",
        f"lig={idx['lig_contacts_count']} atoms, rec={idx['rec_contacts_count']} atoms",
        f"(mode={contact_atom_mode}, cutoff={contact_interface_cutoff} nm, cap={contact_max_atoms})",
    )
    print("Metadynamics variables added: d1 (COM distance) + c1 (contacts)")
    return system


def extract_atom_indices(
    pdf_file: os.path,
    contact_atom_mode: str = "ca",
    contact_interface_cutoff_nm: float = 0.8,
    contact_max_atoms_per_partner: int = 120,
):
    """Extract ligand/receptor atom indices plus compact interface contact subsets for PLUMED."""
    u = mda.Universe(pdf_file)

    lig = u.select_atoms("chainID A")
    rec = u.select_atoms("chainID B or chainID C")

    lig_backbone = u.select_atoms("(chainID A) and backbone")
    rec_backbone = u.select_atoms("((chainID B) or (chainID C)) and backbone")

    lig_contacts, rec_contacts = select_interface_contact_atoms(
        u,
        lig_sel="chainID A",
        rec_sel="chainID B or chainID C",
        atom_mode=contact_atom_mode,
        interface_cutoff_nm=contact_interface_cutoff_nm,
        max_atoms_per_partner=contact_max_atoms_per_partner,
    )

    # PLUMED atom indexing is 1-based and corresponds to the .ids attribute in MDAnalysis AtomGroups, which is what we want for both PLUMED and PyMOL.
    lig_plumed_backbone = ",".join(map(str, lig_backbone.indices + 1))
    rec_plumed_backbone = ",".join(map(str, rec_backbone.indices + 1))
    lig_plumed_contacts = ",".join(map(str, lig_contacts.indices + 1))
    rec_plumed_contacts = ",".join(map(str, rec_contacts.indices + 1))

    return {
        "lig_backbone": lig_plumed_backbone,
        "rec_backbone": rec_plumed_backbone,
        "lig_contacts": lig_plumed_contacts,
        "rec_contacts": rec_plumed_contacts,
        "lig_contacts_count": int(lig_contacts.n_atoms),
        "rec_contacts_count": int(rec_contacts.n_atoms),
        "lig_min": int(lig.indices.min() + 1),
        "lig_max": int(lig.indices.max() + 1),
        "rec_min": int(rec.indices.min() + 1),
        "rec_max": int(rec.indices.max() + 1),
    }


def _select_contact_mode_atoms(u, chain_selection: str, atom_mode: str):
    atom_mode = atom_mode.lower()
    if atom_mode == "all":
        return u.select_atoms(chain_selection)
    if atom_mode == "heavy":
        return u.select_atoms(f"({chain_selection}) and not name H*")
    if atom_mode == "backbone":
        return u.select_atoms(f"({chain_selection}) and backbone")
    if atom_mode == "ca":
        return u.select_atoms(f"({chain_selection}) and name CA")
    raise ValueError(f"Unknown CONTACT_ATOM_MODE='{atom_mode}'. Use one of: ca, backbone, heavy, all")


def _cap_by_smallest_min_distance(atom_group, min_distances, max_atoms):
    if max_atoms <= 0 or atom_group.n_atoms <= max_atoms:
        return atom_group
    order = np.argsort(min_distances)
    keep = np.sort(order[:max_atoms])
    return atom_group[keep]


def select_interface_contact_atoms(
    u,
    lig_sel: str,
    rec_sel: str,
    atom_mode: str = "ca",
    interface_cutoff_nm: float = 0.8,
    max_atoms_per_partner: int = 120,
):
    """Pick compact ligand/receptor atom groups for contact CV based on interface proximity."""
    lig_mode = _select_contact_mode_atoms(u, lig_sel, atom_mode)
    rec_mode = _select_contact_mode_atoms(u, rec_sel, atom_mode)

    if lig_mode.n_atoms == 0 or rec_mode.n_atoms == 0:
        lig_mode = u.select_atoms(f"({lig_sel}) and backbone")
        rec_mode = u.select_atoms(f"({rec_sel}) and backbone")

    if lig_mode.n_atoms == 0 or rec_mode.n_atoms == 0:
        lig_mode = u.select_atoms(lig_sel)
        rec_mode = u.select_atoms(rec_sel)

    distances = mda.lib.distances.distance_array(lig_mode.positions, rec_mode.positions)
    min_lig = distances.min(axis=1)
    min_rec = distances.min(axis=0)

    # Only take distances within a generous cutoff to avoid keeping too many atoms for large systems.
    interface_cutoff_a = float(interface_cutoff_nm) * 10.0
    lig_mask = min_lig <= interface_cutoff_a
    rec_mask = min_rec <= interface_cutoff_a
    lig_interface = lig_mode[lig_mask]
    rec_interface = rec_mode[rec_mask]

    # If no atoms pass the cutoff, take the closest ones up to the max_atoms_per_partner limit.
    # Happens if no real complex formed during equilibration or if cutoff is too tight. This ensures the CVs are always defined.
    if lig_interface.n_atoms == 0:
        lig_interface = lig_mode[np.sort(np.argsort(min_lig)[: min(max_atoms_per_partner, lig_mode.n_atoms)])]
    if rec_interface.n_atoms == 0:
        rec_interface = rec_mode[np.sort(np.argsort(min_rec)[: min(max_atoms_per_partner, rec_mode.n_atoms)])]

    # Cap the number of atoms by the smallest minimum distance to avoid exceeding the max_atoms_per_partner limit.
    if lig_mask.any():
        lig_interface = _cap_by_smallest_min_distance(lig_interface, min_lig[lig_mask], max_atoms_per_partner)
    else:
        lig_interface = _cap_by_smallest_min_distance(lig_interface, min_lig[np.argsort(min_lig)[:lig_interface.n_atoms]], max_atoms_per_partner)

    if rec_mask.any():
        rec_interface = _cap_by_smallest_min_distance(rec_interface, min_rec[rec_mask], max_atoms_per_partner)
    else:
        rec_interface = _cap_by_smallest_min_distance(rec_interface, min_rec[np.argsort(min_rec)[:rec_interface.n_atoms]], max_atoms_per_partner)

    if lig_interface.n_atoms == 0 or rec_interface.n_atoms == 0:
        lig_interface = _cap_by_smallest_min_distance(lig_mode, min_lig, max_atoms_per_partner)
        rec_interface = _cap_by_smallest_min_distance(rec_mode, min_rec, max_atoms_per_partner)

    return lig_interface, rec_interface
