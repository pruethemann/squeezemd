from openmmplumed import PlumedForce
import os
from openmm.unit import kelvin
import MDAnalysis as mda
import numpy as np

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

    # get relevant atom indexes
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


from typing import Optional
from openmm import (
    Context, System,
    NonbondedForce, CustomNonbondedForce, CustomBondForce,
    CustomExternalForce, CustomAngleForce, CustomTorsionForce,
    HarmonicBondForce, HarmonicAngleForce, PeriodicTorsionForce,
    CMMotionRemover, MonteCarloBarostat
)


def save_active_forces(
    system: System,
    context: Optional[Context] = None,
    logfile: str = "log.txt",
    max_examples: int = 3
) -> None:
    """
    Save all active OpenMM forces and parameters to a log file.

    - If `context` is provided, active global parameter values are printed.
    - Otherwise, default parameter values are printed.
    """

    def _ctx_value(name: str):
        if context is None:
            return None
        try:
            return context.getParameter(name)
        except Exception:
            return None

    with open(logfile, "w") as f:
        def write(msg=""):
            f.write(msg + "\n")

        n_forces = system.getNumForces()
        write(f"=== OpenMM Forces in System: {n_forces} ===")

        for i in range(n_forces):
            force = system.getForce(i)
            cname = force.__class__.__name__
            try:
                fg = force.getForceGroup()
            except Exception:
                fg = "n/a"

            write(f"\n[{i:02d}] {cname} (forceGroup={fg})")

            # ---------- Global parameters (Custom* forces) ----------
            if hasattr(force, "getNumGlobalParameters"):
                ng = force.getNumGlobalParameters()
                if ng > 0:
                    write(f"  Global parameters ({ng}):")
                    for gi in range(ng):
                        name = force.getGlobalParameterName(gi)
                        default = None
                        if hasattr(force, "getGlobalParameterDefaultValue"):
                            default = force.getGlobalParameterDefaultValue(gi)
                        active = _ctx_value(name)
                        if active is None:
                            write(f"    - {name}: default={default}")
                        else:
                            write(f"    - {name}: active={active} (default={default})")

            # ---------- Force-specific details ----------
            if isinstance(force, NonbondedForce):
                write(f"  Nonbonded method: {force.getNonbondedMethod()}")
                write(f"  Cutoff: {force.getCutoffDistance()}")
                write(f"  Ewald error tol: {force.getEwaldErrorTolerance()}")
                write(f"  Dispersion correction: {force.getUseDispersionCorrection()}")
                write(f"  Num particles: {force.getNumParticles()}")

                for p in range(min(max_examples, force.getNumParticles())):
                    q, sig, eps = force.getParticleParameters(p)
                    write(f"    particle[{p}] q={q} sigma={sig} epsilon={eps}")

            elif isinstance(force, CustomExternalForce):
                write(f"  Energy: {force.getEnergyFunction()}")
                write(f"  Per-particle parameters: "
                      f"{[force.getPerParticleParameterName(j) for j in range(force.getNumPerParticleParameters())]}")
                write(f"  Restrained particles: {force.getNumParticles()}")

                for p in range(min(max_examples, force.getNumParticles())):
                    idx, params = force.getParticleParameters(p)
                    write(f"    particle[{p}] atomIndex={idx} params={params}")

            elif isinstance(force, CustomNonbondedForce):
                write(f"  Energy: {force.getEnergyFunction()}")
                write(f"  Num particles: {force.getNumParticles()}")

            elif isinstance(force, HarmonicBondForce):
                write(f"  Num bonds: {force.getNumBonds()}")
                for b in range(min(max_examples, force.getNumBonds())):
                    a1, a2, length, k = force.getBondParameters(b)
                    write(f"    bond[{b}] ({a1},{a2}) length={length} k={k}")

            elif isinstance(force, HarmonicAngleForce):
                write(f"  Num angles: {force.getNumAngles()}")

            elif isinstance(force, PeriodicTorsionForce):
                write(f"  Num torsions: {force.getNumTorsions()}")

            elif isinstance(force, MonteCarloBarostat):
                write(f"  Pressure: {force.getDefaultPressure()}")
                write(f"  Temperature: {force.getDefaultTemperature()}")
                write(f"  Frequency: {force.getFrequency()} steps")

            elif isinstance(force, CMMotionRemover):
                write(f"  Frequency: {force.getFrequency()} steps")

            else:
                # Fallback
                if hasattr(force, "getEnergyFunction"):
                    try:
                        write(f"  Energy: {force.getEnergyFunction()}")
                    except Exception:
                        pass

        write("\n=== End force listing ===")

##### Inactive or depracted functions for metadynamics

def add_metadynamics_forces_centerofmass_old(params, system, args, T=300):
    """
    This version does work but conciders all atoms in COM. I will try to reduce to C_alphas only
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

def add_metadynamics_forces_centerofmass_contacts(params, system, args, T=300):

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

    # get relevant atom indexes
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

def extract_atom_indices_extended(pdf_file: os.path, cutoff = 5.0):

    u = mda.Universe(pdf_file)

    # Get all atoms
    lig = u.select_atoms("chainID A") 
    rec = u.select_atoms("chainID B or chainID C")

    # get all heavy atoms of lig and rec
    lig_heavy = u.select_atoms("chainID A and not name H*")
    rec_heavy = u.select_atoms("(chainID B or chainID C) and not name H*")

    # get only C alphas to reduce computational cost of COM calculation
    lig_ca = u.select_atoms("(chainID A) and (name CA)")
    rec_ca = u.select_atoms("((chainID B) or (chainID C)) and (name CA)")

    # Compute distance matrix between all atoms of the two chains
    dist = mda.lib.distances.distance_array(lig_heavy.positions, rec_heavy.positions)

    # Boolean masks of interface atoms
    lig_interface_mask = np.any(dist < cutoff, axis=1)
    rec_interface_mask = np.any(dist < cutoff, axis=0)

    # Interface atoms selections
    lig_interface = lig_heavy[lig_interface_mask]
    rec_interface = rec_heavy[rec_interface_mask]

    # Print in PLUMED-friendly format. Add +1 because plumed starts at atom id 1 and not 0
    lig_plumed_ca = ",".join(map(str, lig_ca.indices + 1))
    rec_plumed_ca = ",".join(map(str, rec_ca.indices + 1))
    lig_interface_plumed = ",".join(map(str, lig_interface.indices + 1))
    rec_interface_plumed = ",".join(map(str, rec_interface.indices + 1))

    atom_indices = {'lig_ca': lig_plumed_ca,
                    'rec_ca': rec_plumed_ca,
                    'lig_interface_index':lig_interface_plumed,
                    'rec_interface_index':rec_interface_plumed,
                    'lig_min':lig.indices.min() + 1,
                    'lig_max':lig.indices.max() + 1,
                    'rec_min':rec.indices.min() + 1,
                    'rec_max':rec.indices.max() + 1,
    }

    return atom_indices