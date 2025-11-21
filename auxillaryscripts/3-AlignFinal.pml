from pymol import cmd

def align_structures(n_structures=10, cutoff=3.0):
    """
    Load, clean, color, and align multiple structures (1–N),
    keeping only surface water and ions *per protein* within a distance cutoff.
    """

    # Load all structures
    for i in range(1, n_structures + 1):
        obj = f"topo_center_{i}"
        filename = f"{obj}.pdb"
        cmd.load(filename, obj)

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
    cmd.save("alignment.pse")
    print("✅ Alignment complete. Saved as alignment.pse")


align_structures(n_structures=3)
