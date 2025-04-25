
import mdtraj as md
traj = md.load('output/traj.h5')
traj.save_dcd('output/traj.dcd')