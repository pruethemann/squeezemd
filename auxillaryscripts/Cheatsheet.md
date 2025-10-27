# Copy all topo-center

mkdir -p collected_topos && n=1 && find . -type f -name "topo_center.pdb" -print0 | while IFS= read -r -d '' file; do
  cp "$file" "collected_topos/topo_center_${n}.pdb"
  ((n++))
done


# Adjust ion size in pymol
alter elem Na, vdw=0.7
rebuild

# Show helostasin Args
sele arg, (resid 578,863,504,630)


select C_term, resid 80:100

select near_cterm, (byres (polymer within 10 of C_term)) or (name NA within 10 of C_term)