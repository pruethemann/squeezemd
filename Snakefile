import sys
import numpy as np
import pandas as pd
from os import path

# ======================
# Configuration Handling
# ======================

# Define paths to key files and directories
simulation_data = path.join('config', 'simulations.csv')
configfile: path.join('config', 'params.yml')       # Configuration file for Snakemake
free_energy_settings =  path.join('config', 'mmgbsa.in')

# Retrieve essential parameters from the config
number_frames = config.get("number_frames", 100)  # Default to 100 if not set
replicates = config.get("replicates", 1)  # Default to 1 if not set

# ==================================
# Random Seed Generation for Replicates
# ==================================
np.random.seed(23)  # Ensure reproducibility
seeds = np.random.randint(100, 1000, size=replicates)
representative_seed = seeds[0] # First seeds is used as representation for surface interaction analysis

# ==================================
# Simulation Data Preprocessing
# ==================================
# Load simulation conditions and preprocess data

simulations_df = pd.read_csv(simulation_data)
simulations_df['complex'] = simulations_df['target'] + '_' + simulations_df['ligand']
simulations_df['name'] = simulations_df['complex'] + '_' + simulations_df['mutation_all']
simulations_df.set_index('name', inplace=True)

# Define global variables for easy access
complexes = simulations_df.complex.unique()
mutations = simulations_df.mutation_all.unique()

rule proteinInteraction:
    input:
        'results/martin/interactions.parquet',
        #'results/fingerprints/interactions.parquet',
        expand('results/interactionSurface/{complex}.{mutation}.interaction.pdb', complex=complexes, mutation=mutations),
        #expand('{complex}/{mutation}/{seed}/fingerprint/fingerprint.parquet',complex=complexes,mutation=mutations,seed=seeds),
        expand('{complex}/{mutation}/mutation.pdb', complex=complexes, mutation=mutations),
        expand('{complex}/{mutation}/{seed}/MD/trajectory.h5', complex=complexes, mutation=mutations, seed=seeds),
        expand('{complex}/{mutation}/{seed}/analysis/RMSF.svg', complex=complexes, mutation=mutations, seed=seeds),

rule protein:
    input:
        expand('{complex}/{mutation}/{seed}/MD/traj_center.dcd', complex=complexes, mutation=mutations, seed=seeds),
        expand('{complex}/{mutation}/{seed}/analysis/RMSF.svg', complex=complexes, mutation=mutations, seed=seeds),

rule Mutagensis:
    input:
        'config/params.yml',
        'config/simulations.csv'
    output:
        '{complex}/{mutation}/mutation.pdb'
    params:
        out_dir=directory('{complex}/{mutation}'),
        pdb= lambda wildcards: simulations_df.loc[f'{wildcards.complex}_{wildcards.mutation}']['input'],
        foldX='foldx_20241231'
    log:
        '{complex}/{mutation}/mutation.log'
    shell:
        """
        if [[ {wildcards.mutation} = WT ]]
        then    # No mutagenisis
            cp {params.pdb} {output}
        else    # Mutate
            1_mutation.py --mutation {wildcards.mutation} --ligand {params.pdb} --output {wildcards.complex}/{wildcards.mutation}/mutant_file.txt
            cp {params.pdb} {params.out_dir}/WT.pdb
            {params.foldX} --command=BuildModel --pdb-dir="{params.out_dir}" --pdb=WT.pdb --mutant-file={wildcards.complex}/{wildcards.mutation}/mutant_file.txt --output-dir="{params.out_dir}" --rotabaseLocation ~/tools/foldx/rotabase.txt > {log} || true
            mv {params.out_dir}/WT_1.pdb {output}
        fi
        """

rule MD:
    input:
        md_settings=ancient('config/params.yml'),
        pdb='{complex}/{mutation}/mutation.pdb',
    output:
        topo = '{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj=temp('{complex}/{mutation}/{seed}/MD/trajectory.h5'),
        stats='{complex}/{mutation}/{seed}/MD/MDStats.csv',
        params='{complex}/{mutation}/{seed}/MD/params.yml',
    resources:
        gpu=1
    priority:
        3
    shell:
        """
        2_MD.py --pdb {input.pdb} \
                --topo {output.topo} \
                --traj {output.traj} \
                --md_settings {input.md_settings} \
                --params {output.params} \
                --seed {wildcards.seed} \
                --stats {output.stats}
        """


rule centerMDTraj:
    input:
        topo = '{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj='{complex}/{mutation}/{seed}/MD/trajectory.h5',
    output:
        topo_center = '{complex}/{mutation}/{seed}/MD/topo_center.pdb',
        traj_center='{complex}/{mutation}/{seed}/MD/traj_center.dcd',
    priority:
        3
    shell:
        """
        4_centerMDTraj.py   --topo {input.topo} \
                            --traj {input.traj} \
                            --traj_center {output.traj_center} \
                            --topo_center {output.topo_center}
        """

rule DescriptiveTrajAnalysis:
    input:
        topo='{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj='{complex}/{mutation}/{seed}/MD/traj_center.dcd',
        stats='{complex}/{mutation}/{seed}/MD/MDStats.csv'
    output:
        rmsf= report('{complex}/{mutation}/{seed}/analysis/RMSF.svg',caption="RMSF.rst",category="RMSF",labels={"Complex": "{complex}", "Mutation": "{mutation}"}),
        rmsd= report('{complex}/{mutation}/{seed}/analysis/RMSD.svg',caption="RMSF.rst",category="RMSD",labels={"Complex": "{complex}", "Mutation": "{mutation}"}),
        bfactors='{complex}/{mutation}/{seed}/analysis/bfactors.pdb',
        stats='{complex}/{mutation}/{seed}/analysis/Stats.svg',
    params:
        number_frames = config.get('number_frames')
    shell:
        """
        3_ExplorativeTrajectoryAnalysis.py --topo {input.topo} \
                                                   --traj {input.traj} \
                                                   --stats {input.stats} \
                                                   --rmsf {output.rmsf} \
                                                   --bfactors {output.bfactors} \
                                                   --rmsd {output.rmsd} \
                                                   --fig_stats {output.stats}
        """

rule InteractionAnalyzerMartin:
    input:
        topo = '{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj='{complex}/{mutation}/{seed}/MD/traj_center.dcd',
        pdb='{complex}/{mutation}/mutation.pdb',
    output:
        frame_pdb='{complex}/{mutation}/{seed}/frames/frame_1.pdb',
        ligand_csv=ensure('{complex}/{mutation}/{seed}/frames/lig/1.csv', non_empty=True),
        receptor_csv=ensure('{complex}/{mutation}/{seed}/frames/rec/1.csv', non_empty=True),
        dir=directory('{complex}/{mutation}/{seed}/frames/'),
        final='{complex}/{mutation}/{seed}/MD/final.pdb',
    log:
        '{complex}/{mutation}/{seed}/MD/center.log'
    threads: 4
    shell:
        """
        5_Martin_analyzer.py  --topo {input.topo} \
                              --traj {input.traj} \
                              --pdb {input.pdb} \
                              --n_frames {number_frames} \
                              --final {output.final} \
                              --cpu {threads} \
                              --dir {output.dir}  > {log}
        """


rule GlobalMartinAnalysis:
    input:
        dirs=expand('{complex}/{mutation}/{seed}/frames/', complex=complexes, mutation=mutations, seed=seeds),
    output:
        interactions = report('results/martin/interactions.parquet', caption="posco.rst",category="Posco", labels=({"Name": "All Interactions", "Type": "List"})),
        residueEnergy= report('results/martin/residueEnergy.svg', caption="posco.rst",category="Posco", labels=({"Name": "Per Residue Interaction Energy", "Type": "Plot"})),
        totalEnergy= report('results/martin/totalEnergy.svg', caption="posco.rst",category="Posco", labels=({"Name": "Total Energy", "Type": "Plot"}))
    shell:
        """
        6_GlobalMartinInteractions.py  --dirs {input.dirs} \
                                       --interactions {output.interactions} \
                                       --residueEnergy {output.residueEnergy} \
                                       --totalEnergy {output.totalEnergy}
        """

rule InteractionSurface:
    input:
        final_frame = f'{{complex}}/{{mutation}}/{representative_seed}/MD/topo_center.pdb',
        interactions= 'results/martin/interactions.parquet',
    output:
        bfactor_pdbs = 'results/interactionSurface/{complex}.{mutation}.interaction.pdb',
        pymol_cmd = 'results/interactionSurface/{complex}.{mutation}.pml',
        pymol = report('results/interactionSurface/{complex}.{mutation}.final.pse',caption="RMSF.rst",category="PyMol",labels=({"Complex": "{complex}", "Mutation": "{mutation}", "Type": "PyMol"})),
        surface= report('results/interactionSurface/{complex}.{mutation}.png',caption="RMSF.rst",category="PyMol", labels=({"Complex": "{complex}", "Mutation": "{mutation}", "Type": "Image"}))
    params:
        representative_seed = seeds[0],
    shell:
        """
        10_InteractionSurface.py --interactions {input.interactions} \
                                 --seed {params.representative_seed} \
                                 --mutation {wildcards.mutation} \
                                 --frames {input.final_frame} \
                                 --complex {wildcards.complex}
        pymol -cQ {output.pymol_cmd}
        """


rule interactionFingerprint:
    input:
        topo='{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj='{complex}/{mutation}/{seed}/MD/traj_center.dcd',
    output:
        report('{complex}/{mutation}/{seed}/fingerprint/fingerprint.parquet', labels=({"Name": "Interaction Analysis Prolif", "Type": "List"}),caption="RMSF.rst",category="Interaction Fingerprint")
    params:
        number_frames = config.get('number_frames'),
    threads: 4
    shell:
        """
        7_interactionFingerprint.py --topo {input.topo} \
                                    --traj  {input.traj} \
                                    --output {output} \
                                    --n_frames {params.number_frames} \
                                    --threads {threads} \
                                    --complex {wildcards.complex} \
                                    --mutation {wildcards.mutation} \
                                    --seed {wildcards.seed}
        """

rule GlobalFingerprintAnalysis:
    input:
        fingerprints=expand('{complex}/{mutation}/{seed}/fingerprint/fingerprint.parquet', complex=complexes, mutation=mutations, seed=seeds),
    output:
        interactions = report('results/fingerprints/interactions.parquet',labels=({"Name": "Interaction Analysis Prolif", "Type": "List"}),caption="RMSF.rst",category="Interaction Fingerprint"),
    params:
        n_frames = config.get('number_frames')
    shell:
        """
        8_GlobalFinterprintAnalysis.py  --fingerprints {input.fingerprints} \
                                        --interactions {output.interactions} \
                                        --n_frames {params.n_frames}
        """

onsuccess:
    # If only single protein has been analyzed
    if 'protein' in sys.argv:
        print("Workflow finished, no error. Report for single protein will be generated")
        shell("squeeze protein --report report.html")
        print("YES")
    else:
        print("Workflow finished, no error. Report for protein-protein Interaciton will be generated")
        shell("squeeze --report report.html")
