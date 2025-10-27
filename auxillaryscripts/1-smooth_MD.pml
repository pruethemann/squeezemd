# Script loads dcd trajectory and performs a smoothing of the MD


load topo_center.pdb, md
load_traj traj_center.dcd, md,start=1, interval=5

remove solvent
remove name Cl
color aquamarine, chain A
color lightblue, chain B
alter elem Na, vdw=0.7
smooth md, 30, 50


save C-terminus_full.pse

