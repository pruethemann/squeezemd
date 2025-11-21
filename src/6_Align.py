#!/usr/bin/env python


from pymol import cmd
import argparse

def align_structures(input_structures, output, n_structures=10, cutoff=3.0):
    """
    Load, clean, color, and align multiple structures (1–N),
    keeping only surface water and ions *per protein* within a distance cutoff.
    """
    import os
    # Load all structures
    for struct_file in input_structures:
        # Import every final structure
        obj = os.path.basename(struct_file)[:-3]
        cmd.load(struct_file, obj)     

        # Define selections for this object
        protein_sel = f"({obj} and polymer.protein)"
        surface_water = f"({obj} and resn HOH within {cutoff} of {protein_sel})"
        surface_ions = f"({obj} and (resn NA+ or resn CL-) within {cutoff} of {protein_sel})"

        # Remove all atoms in this object not part of the protein, surface water, or nearby ions
        cmd.remove(f"{obj} and not ({protein_sel} or {surface_water} or {surface_ions})")

    # Adjust van der Waals radius for sodium for visualization
    cmd.alter("elem Na", "vdw=0.7")
    cmd.alter("elem Cl", "vdw=3")

    # Color chains consistently
    cmd.color("aquamarine", "chain A")
    cmd.color("lightblue", "chain B")

    # Align all structures to the first
    for i in range(2, n_structures + 1):
        mobile = f"topo_center_{i}"
        target = "topo_center_1"
        cmd.align(mobile, target)

    # Save aligned session
    cmd.save(output)
    print("✅ Alignment complete. Saved as alignment.pse")


def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', nargs='+', required=False)
    parser.add_argument('--traj', nargs='+', required=False)

    # Output
    parser.add_argument('--output', required=False, default='rmsf.svg', help='')

    return parser.parse_args()

def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--input', nargs='+', required=False)

    # Output
    parser.add_argument('--output', required=False, default='align.pse', help='')

    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()
    n_structures = len(args.input)
    align_structures(args.input, n_structures=n_structures)
