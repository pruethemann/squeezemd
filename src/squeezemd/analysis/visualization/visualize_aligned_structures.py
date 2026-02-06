#!/usr/bin/env python

"""PyMOL alignment helper for final structure comparison.

Run inside PyMOL (batch) to load multiple structures, filter to
protein + nearby solvent/ions, align to a reference, and save a session.
"""

from pymol import cmd
import argparse, os

def align_structures(input_structures, output, cutoff=3.0):
    """
    Load, clean, color, and align multiple structures (1–N),
    keeping only surface water and ions *per protein* within a distance cutoff.
    """

    objects = []

    print(input_structures)
   
    # Load all structures
    for struct_file in input_structures:
        mutation = struct_file.split("/")[-5]
        seed = struct_file.split("/")[-4]
        complex = struct_file.split("/")[-6]

        # Import every final structure
        obj = complex + '_' + '_' + seed + '_' + mutation
        obj = obj + "_" + mutation
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

    print("OBJECTS")
    print(objects)

    # Only perform alignment if more than 1 structure
    if len(objects) > 1:
        for mobile in objects[1:]:
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


def main():
    print("Hello world")
    args = parse_arguments()
    align_structures(args.input, args.output)



print("Hello world")
args = parse_arguments()
align_structures(args.input, args.output)

