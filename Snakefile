import numpy as np
import pandas as pd
from os import path

### USER INPUT
if 'job' not in config:
    print("ATTENTION: No job defined. Demo will be executed!")
    config["job"] = 'demo'
job = path.join('jobs', config["job"])

# 1. Define all paths
simulation_data = path.join(job, 'simulations.csv')
configfile: path.join(job, 'params.yml')       # Configuration file for Snakemake
md_settings = path.join(job, 'params.yml')
free_energy_settings =  path.join(job, 'mmgbsa.in')

number_frames = config['number_frames']
replicates = config['replicates']

# 2. Determine seed
np.random.seed(23)
seeds = np.random.randint(100,1000, size=config['replicates'])
#seeds = [767,943]
print(seeds)

# Import all simulation conditions
simulations_df = pd.read_csv(simulation_data)
print(simulations_df)
simulations_df['complex'] = simulations_df['target'] + '_' + simulations_df['ligand']
simulations_df['name'] = simulations_df['complex'] + '_' + simulations_df['mutation_all']
complexes = simulations_df.complex
mutations = simulations_df.mutation_all
frames = list(range(number_frames))

simulations_df['complex'] = simulations_df['target'] + '_' + simulations_df['ligand']
simulations_df['name'] = simulations_df['complex'] + '_' + simulations_df['mutation_all']
simulations_df.set_index('name', inplace=True)

rule proteinInteraction:
    input:
        expand('output/{job_id}/{complex}/{mutation}/mutation.pdb', job_id=config["job"], complex=complexes, mutation=mutations),
        #expand('output/{job_id}/{complex}/{mutation}/amber/amber.pdb',job_id=config["job"], complex=complexes, mutation=mutations),
        expand('output/{job_id}/{complex}/{mutation}/{seed}/MD/traj_center.dcd', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        #expand('output/{job_id}/{complex}/{mutation}/{seed}/MD/center/topo.pdb', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        expand('output/{job_id}/{complex}/{mutation}/{seed}/analysis/RMSF.svg', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        expand('output/{job_id}/{complex}/{mutation}/{seed}/frames/lig/1.csv', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds, frame_id=frames),
        expand('output/{job_id}/results/martin/interactions.csv', job_id=config["job"]),
        expand('output/{job_id}/{complex}/{mutation}/{seed}/fingerprint/fingerprint.feather', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        expand('output/{job_id}/results/fingerprints/interactions.csv', job_id=config["job"]),
        expand('output/{job_id}/results/interactions/.checkpoint', job_id=config["job"])





rule SetupMutagesis:
    input:
        csv=simulation_data,
        default=ancient(md_settings)
    output:
        'output/{job_id}/{complex}/{mutation}/mutant_file.txt'
    params:
        job_id=config['job'],
    shell:
        "1_mutation.py --mutation {wildcards.mutation} --output {output}"

rule Mutagensis:
    input:
        'output/{job_id}/{complex}/{mutation}/mutant_file.txt'
    output:
        'output/{job_id}/{complex}/{mutation}/mutation.pdb'
    params:
        out_dir=directory('output/{job_id}/{complex}/{mutation}'),
        pdb= lambda wildcards: simulations_df.loc[f'{wildcards.complex}_{wildcards.mutation}']['input'],
        foldX='foldx_20241231'
    log:
        'output/{job_id}/{complex}/{mutation}/mutation.log'
    shell:
        """
        if [[ {wildcards.mutation} = WT ]]
        then    # No mutagenisis
            cp {params.pdb} {output}
        else    # Mutate
            cp {params.pdb} {params.out_dir}/WT.pdb
            {params.foldX} --command=BuildModel --pdb-dir="{params.out_dir}" --pdb=WT.pdb --mutant-file="{input}" --output-dir="{params.out_dir}" --rotabaseLocation ~/tools/foldx/rotabase.txt > {log} || true
            mv {params.out_dir}/WT_1.pdb {output}
        fi
        """

rule MD:
    input:
        md_settings=ancient(md_settings),
        input_pdb='output/{job_id}/{complex}/{mutation}/mutation.pdb',
    output:
        topo = 'output/{job_id}/{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj=temp('output/{job_id}/{complex}/{mutation}/{seed}/MD/trajectory.h5'),
        stats='output/{job_id}/{complex}/{mutation}/{seed}/MD/MDStats.csv',
        params='output/{job_id}/{complex}/{mutation}/{seed}/MD/params.yml',
    params:
        job_id=config['job']
    resources:
        gpu=1
    priority:
        3
    shell:
        """
        2_MD.py --input_pdb {input.input_pdb} \
                --topo {output.topo} \
                --traj {output.traj} \
                --md_settings {input.md_settings} \
                --params {output.params} \
                --seed {wildcards.seed} \
                --stats {output.stats}
        """


rule centerMDTraj:
    input:
        topo = 'output/{job_id}/{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj='output/{job_id}/{complex}/{mutation}/{seed}/MD/trajectory.h5',
    output:
        topo_center = 'output/{job_id}/{complex}/{mutation}/{seed}/MD/topo_center.pdb',
        traj_center='output/{job_id}/{complex}/{mutation}/{seed}/MD/traj_center.dcd',
    params:
        job_id=config['job']
    resources:
        gpu=1
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
        topo = 'output/{job_id}/{complex}/{mutation}/{seed}/MD/topo_center.pdb',
        traj='output/{job_id}/{complex}/{mutation}/{seed}/MD/traj_center.dcd',
        stats='output/{job_id}/{complex}/{mutation}/{seed}/MD/MDStats.csv'
    output:
        rmsf=report('output/{job_id}/{complex}/{mutation}/{seed}/analysis/RMSF.svg',caption="config/RMSF.rst", category="Trajectory", subcategory='RMSF'),
        bfactors='output/{job_id}/{complex}/{mutation}/{seed}/analysis/bfactors.pdb',
        rmsd=report('output/{job_id}/{complex}/{mutation}/{seed}/analysis/RMSD.svg',caption="config/RMSD.rst",category="Trajectory", subcategory='RMSD'),
        stats='output/{job_id}/{complex}/{mutation}/{seed}/analysis/Stats.svg',
    params:
        number_frames = config['number_frames']
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
        topo='output/{job_id}/{complex}/{mutation}/{seed}/MD/topo_center.pdb',
        traj='output/{job_id}/{complex}/{mutation}/{seed}/MD/traj_center.dcd',
    output:
        frame_pdb='output/{job_id}/{complex}/{mutation}/{seed}/frames/frame_1.pdb',
        ligand_csv='output/{job_id}/{complex}/{mutation}/{seed}/frames/lig/1.csv',
        receptor_csv='output/{job_id}/{complex}/{mutation}/{seed}/frames/rec/1.csv',
        dir=directory('output/{job_id}/{complex}/{mutation}/{seed}/frames/'),
        final='output/{job_id}/{complex}/{mutation}/{seed}/MD/final.pdb',
        checkpoint='output/{job_id}/{complex}/{mutation}/{seed}/frames/.done'
    log:
        'output/{job_id}/{complex}/{mutation}/{seed}/MD/center.log'
    shell:
        """
        5_Martin_analyzer.py  --topo {input.topo} \
                              --traj {input.traj} \
                              --n_frames {number_frames} \
                              --final {output.final} \
                              --dir {output.dir}  > {log}
        touch {output.checkpoint}
        """

rule GlobalMartinAnalysis:
    input:
        dirs=expand('output/{job_id}/{complex}/{mutation}/{seed}/frames/',job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        checkpoint=expand('output/{job_id}/{complex}/{mutation}/{seed}/frames/.done',job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds)
    output:
        interactions = report('output/{job_id}/results/martin/interactions.csv',caption="config/RMSF.rst",category="Interaction Martin"),
        #energy_diff = report('output/{job_id}/results/martin/energy_diff.csv',caption="report/RMSF.rst",category="Energy Differences"),
        dir=directory('output/{job_id}/results/martin/')
    shell:
        """
        6_GlobalMartinInteractions.py  --dirs {input.dirs} \
                                --output {output.dir} \
                                --interactions {output.interactions}
        """

rule interactionFingerprint:
    input:
        topo = 'output/{job_id}/{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj='output/{job_id}/{complex}/{mutation}/{seed}/MD/traj_center.dcd',
    output:
        'output/{job_id}/{complex}/{mutation}/{seed}/fingerprint/fingerprint.feather',
    params:
        number_frames = config['number_frames']
    shell:
        """
        7_interactionFingerprint.py --topo {input.topo} \
                                    --traj  {input.traj} \
                                    --output {output} \
                                    --n_frames {params.number_frames}
        """

rule GlobalFingerprintAnalysis:
    input:
        fingerprints=expand('output/{job_id}/{complex}/{mutation}/{seed}/fingerprint/fingerprint.feather',job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
    output:
        interactions = report('output/{job_id}/results/fingerprints/interactions.csv',caption="config/RMSF.rst",category="Interaction Martin"),
        figure='output/{job_id}/results/fingerprints/interactions.svg'
    shell:
        """
        8_GlobalFinterprintAnalysis.py  --fingerprints {input.fingerprints} \
                              --interactions {output.interactions} \
                              --figure {output.figure}
        """

rule InteractionSurface:
    input:
        interactions = report('output/{job_id}/results/martin/interactions.csv',caption="config/RMSF.rst",category="Interaction Martin"),
    output:
        dir=directory('output/{job_id}/results/interactionSurface/'),
        checkpoint='output/{job_id}/results/interactions/.checkpoint'
    params:
        job=config["job"],
        simulations='output/demo/simulations.csv'
    shell:
        """
        10_InteractionAminoAcids.py --output {output.dir} \
                                            --interactions {input.interactions} \
                                            --simulations {params.simulations} \
                                            --checkpoint {output.checkpoint}
        """
