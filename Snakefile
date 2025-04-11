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
        'results/fingerprints/interactions.parquet',
        'results/metaReport.html',
        #expand('results/interactionSurface/{complex}.{mutation}.interaction.pdb', complex=complexes, mutation=mutations),
        expand('{complex}/{mutation}/{seed}/fingerprint/fingerprint.parquet',complex=complexes,mutation=mutations,seed=seeds),
        expand('{complex}/{mutation}/{seed}/analysis/RMSF.html', complex=complexes, mutation=mutations, seed=seeds),
        expand('{complex}/{mutation}/{seed}/frames/lig_{i}.pdb', i=range(number_frames), complex=complexes, mutation=mutations, seed=seeds),
        expand('{complex}/{mutation}/{seed}/frames/rec_{i}.pdb', i=range(number_frames), complex=complexes, mutation=mutations, seed=seeds),
        expand('{complex}/{mutation}/{seed}/po-sco/{i}.txt', i=range(number_frames), complex=complexes, mutation=mutations, seed=seeds),
        'results/posco/posco_interactions.parquet',
        'results/posco/posco.html',
        expand('{complex}/{mutation}/{seed}/aquaduct/aquaduct.pse',complex=complexes, mutation=mutations, seed=seeds)

rule protein:
    input:
        expand('{complex}/{mutation}/{seed}/MD/traj_center.dcd', complex=complexes, mutation=mutations, seed=seeds),
        expand('{complex}/{mutation}/{seed}/analysis/RMSF.html', complex=complexes, mutation=mutations, seed=seeds),
        expand('{complex}/{mutation}/{seed}/aquaduct/aquaduct.pse',complex=complexes, mutation=mutations, seed=seeds)


rule metaReport:
    input:
        params='config/params.yml',
        sims='config/simulations.csv'
    output:
        report('results/metaReport.html',caption="RMSF.rst",category="Meta"),
    shell:
        """
        metaReport.py  --param {input.params}  \
                       --sims {input.sims}  \
                       --output {output}  \
        """

rule Mutagensis:
    input:
        'config/params.yml',
        'config/simulations.csv'
    output:
        '{complex}/{mutation}/mutation.pdb'
    params:
        out_dir=directory('{complex}/{mutation}'),
        pdb= lambda wildcards: simulations_df.loc[f'{wildcards.complex}_{wildcards.mutation}']['input'],
        foldX='foldx_20251231'              # Needs to be updated every year
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


rule aquaduct:
    input:
        topo_center = '{complex}/{mutation}/{seed}/MD/topo_center.pdb',
        traj_center='{complex}/{mutation}/{seed}/MD/traj_center.dcd',
        aquaduct_template='config/TEMPLATE_aquaduct.txt',
    output:
        aquaduct='{complex}/{mutation}/{seed}/aquaduct/aquaduct.txt',
        aquaduct_folder=directory('{complex}/{mutation}/{seed}/aquaduct/'),
        pymol_cmd='{complex}/{mutation}/{seed}/aquaduct/6_visualize_results.py',
        pymol_output='{complex}/{mutation}/{seed}/aquaduct/aquaduct.pse',
    priority:
        10
    params:
        directory('{complex}/{mutation}/{seed}')

    shell:
        """
        eval "$(micromamba shell hook --shell bash)" &&
        micromamba activate aquaduct_3 &&
        sed  's|PREFIXINPUT|{params}|g' {input.aquaduct_template} > {output.aquaduct} &&
        sed  -i 's|PREFIXOUT|{output.aquaduct_folder}|g' {output.aquaduct} &&
        valve_run -c {output.aquaduct} &&
        /home/anna/micromamba/envs/pymol/bin/python3 {output.pymol_cmd} --save-session aquaduct.pse
        """


rule DescriptiveTrajAnalysis:
    input:
        topo='{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj='{complex}/{mutation}/{seed}/MD/traj_center.dcd',
        stats='{complex}/{mutation}/{seed}/MD/MDStats.csv'
    output:
        rmsf= report('{complex}/{mutation}/{seed}/analysis/RMSF.html',caption="RMSF.rst",category="RMSF",labels={"Complex": "{complex}", "Mutation": "{mutation}"}),
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


rule ExtractithFrameFromTrajectory:
    input:
        topo='{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj='{complex}/{mutation}/{seed}/MD/traj_center.dcd'
    output:
        lig_frame='{complex}/{mutation}/{seed}/frames/lig_{i}.pdb',
        rec_frame='{complex}/{mutation}/{seed}/frames/rec_{i}.pdb'
    threads: 1
    shell:
        """
        5.1_Posco_ExtractLastFrames.py  --topo {input.topo} \
                                        --traj {input.traj} \
                                        --frame {wildcards.i} \
                                        --lig_frame {output.lig_frame} \
                                        --rec_frame {output.rec_frame}
        """

#         touch {output.lig_frame} && touch {output.rec_frame}

# PosCo
rule posco:
    input:
        lig_frames='{complex}/{mutation}/{seed}/frames/lig_{i}.pdb',
        rec_frames='{complex}/{mutation}/{seed}/frames/rec_{i}.pdb'
    output:
        '{complex}/{mutation}/{seed}/po-sco/{i}.txt'
    threads: 1
    shell:
        """
        po-sco {input.rec_frames} {input.lig_frames} -b  > {output}
        """

# PosCo
rule concat:
    input:
        expand('{complex}/{mutation}/{seed}/po-sco/{i}.txt', i=range(number_frames), complex=complexes, mutation=mutations, seed=seeds),
    output:
        'results/posco/posco_interactions.parquet'
    threads: 1
    shell:
        """
        5.2_Posco_TransformDF.py --input {input} --output {output}
        """

# PosCo
rule PoScoAnalysis:
    input:
        'results/posco/posco_interactions.parquet'
    output:
        report('results/posco/posco.html', caption="posco.rst",category="Pose Scoring", labels=({"Name": "Interactions", "Type": "Plot"})),
    threads: 1
    shell:
        """
        5.3_Posco_Analysis.py --input {input} --output {output}
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
        figure = report('results/fingerprints/fingerprints.html',labels=({"Name": "Interaction Analysis Prolif", "Type": "Plot"}),caption="RMSF.rst",category="Interaction Fingerprint"),
    params:
        n_frames = config.get('number_frames')
    shell:
        """
        8_GlobalFinterprintAnalysis.py  --fingerprints {input.fingerprints} \
                                        --interactions {output.interactions} \
                                        --n_frames {params.n_frames} \
                                        --figure {output.figure}
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
