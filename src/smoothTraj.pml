load topo_center.pdb, md
load_traj traj_center.dcd, md,start=1, stop=500, interval=10

remove solvent
remove name Cl
color aquamarine, chain A
color lightblue, chain B
alter elem Na, vdw=0.7
smooth md, 30, 100


save C-terminus.pse


####
set auto_zoom, off
set defer_builds_mode, 3
set sphere_scale, 0.7

# assume your object name is "traj"
# (change if needed)
set movie_loop, 1
frame first

# total number of frames
set nframes, count_states md

python
from pymol import cmd

obj = "md"        # your object name
residue = 92        # residue of interest
radius = 10.0       # distance in Å
nframes = cmd.count_states(obj)

# make a new selection to update dynamically
cmd.select("nearNa", "none")

for i in range(1, nframes + 1):
    cmd.frame(i)
    # select Na+ ions within radius of residue 92
    cmd.select("nearNa", f"(name NA within {radius} of (resi {residue} and {obj}))")
    # hide all sodium first
    cmd.hide("spheres", f"name NA and {obj}")
    # show only the close ones
    cmd.show("spheres", "nearNa")
    # optional: color the visible ones
    cmd.color("yellow", "nearNa")
    cmd.refresh()
python end



####