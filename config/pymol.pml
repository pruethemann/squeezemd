
# Load the pdb file
load INPUT;
remove solvent
show cartoon;
color grey;
spectrum b, red blue grey;


# Show sticks of interacting reisudes
#color dblue, (resid LIGAND_RESIDS) and chain I;
show sticks, (resid LIGAND_RESIDS) and chain I;

#color red, (resid REC_RESIDS) and not chain I;
show sticks, (resid REC_RESIDS) and not chain I;

# Gradient for Inhibitor
spectrum b, red blue grey, chain I

# Gradient for Inhibitor
spectrum b, red blue grey, not chain I

show surface, chain A
set transparency, 0.1, chain A
extract BD001_TARGET, chain I
save OUTPUT;