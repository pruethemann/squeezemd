
# Load the pdb file
load {input_pdb};
remove solvent
show cartoon;
color grey;
spectrum b, red blue grey;


# Show sticks of interacting reisudes: ligand / Receptor
show sticks, (resid {ligand_resids}) and chain A;
show sticks, (resid {receptor_resids}) and not chain A;

# Gradient for Inhibitor
spectrum b, red blue grey, chain A

# Gradient for receptor
spectrum b, red blue grey, not chain A

# ChainID
show surface, chain not A
set transparency, 0.1, chain not A
extract ligand, chain A
save {output};
png {output_png};
dele all;
