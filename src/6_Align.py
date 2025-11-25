#!/usr/bin/env python

"""
This pymol script is execute within python with "run 6_Align.py --input ... --output align.pse"
or in the terminal with "pymol -c 6_Align.py --input ... --output align.pse
"""

from pymol import cmd
import argparse, os

def align_structures(input_structures, output, cutoff=3.0):
    """
    Load, clean, color, and align multiple structures (1–N),
    keeping only surface water and ions *per protein* within a distance cutoff.
    """

    objects = []
    
    # Load all structures
    for struct_file in input_structures:
        # Import every final structure
        obj = os.path.basename(struct_file)[:-3]
        cmd.load(struct_file, obj) 
        objects.append(obj)    

        # Define selections for this object
        protein_sel = f"({obj} and polymer.protein)"
        surface_water = f"({obj} and resn HOH within {cutoff} of {protein_sel})"
        surface_ions = f"({obj} and (resn NA+ or resn CL-) within {cutoff} of {protein_sel})"
        ligand = f"({obj} and chain X)"

        # Remove all atoms in this object not part of the protein, surface water, or nearby ions
        cmd.remove(f"{obj} and not ({protein_sel} or {surface_water} or {surface_ions} or {ligand})")

    # Adjust van der Waals radius for sodium for visualization
    cmd.alter("elem Na", "vdw=0.7")
    cmd.alter("elem Cl", "vdw=3")

    # Color chains consistently/home/peter/Dropbox/code/squeezemd/src/6_Align.py
    cmd.color("aquamarine", "chain A")
    cmd.color("lightblue", "chain B")

    # Align all structures to the first
    reference = objects[0]

    for mobile in objects[1:]:
        print(reference, mobile)
        cmd.align(mobile, reference)

    # Save aligned session
    cmd.save(output)
    print("✅ Alignment complete. Saved as alignment.pse")

def parse_arguments():
    parser = argparse.ArgumentParser()
    # Input
    parser.add_argument('--input', nargs='+', required=False)
    # Output
    parser.add_argument('--output', required=False, default='align.pse', help='')
    return parser.parse_args()

# Execute the script in pymol
args = parse_arguments()
align_structures(args.input, args.output)
