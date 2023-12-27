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
        expand('output/{job_id}/{complex}/{mutation}/amber/amber.pdb',job_id=config["job"], complex=complexes, mutation=mutations),
        expand('output/{job_id}/{complex}/{mutation}/{seed}/MD/trajectory.dcd', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        expand('output/{job_id}/{complex}/{mutation}/{seed}/MD/center/topo.pdb', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        expand('output/{job_id}/{complex}/{mutation}/{seed}/analysis/RMSF.svg', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        expand('output/{job_id}/results/martin/interactions.csv', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        #expand('output/{job_id}/{complex}/{mutation}/{seed}/frames/lig/1.csv', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds, frame_id=frames),
        #expand('output/{job_id}/results/martin/interactions.csv', job_id=config["job"]),
        #expand('output/{job_id}/{complex}/{mutation}/{seed}/analysis/fingerprint.csv', job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),


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
        foldX="install/foldX/foldx_20231231"
    log:
        'output/{job_id}/{complex}/{mutation}/mutation.log'
    shell:
        """
        if [[ {wildcards.mutation} = WT ]]
        then    # No mutagenisis
            cp {params.pdb} {output}
        else    # Mutate
            cp {params.pdb} {params.out_dir}/WT.pdb
            {params.foldX} --command=BuildModel --pdb-dir="{params.out_dir}" --pdb=WT.pdb --mutant-file="{input}" --output-dir="{params.out_dir}" --rotabaseLocation install/foldX/rotabase.txt > {log} || true
            mv {params.out_dir}/WT_1.pdb {output}
        fi
        """

rule Amber:
    input:
        'output/{job_id}/{complex}/{mutation}/mutation.pdb'
    output:
        amber='output/{job_id}/{complex}/{mutation}/amber/amber.pdb',
        prmtop='output/{job_id}/{complex}/{mutation}/amber/{complex}.prmtop',
        inpcrd='output/{job_id}/{complex}/{mutation}/amber/{complex}.inpcrd',
        leap='output/{job_id}/{complex}/{mutation}/amber/tleap_nosolvent.in',
        mapping='output/{job_id}/{complex}/{mutation}/amber/amber_renum.txt',
        tleappdb='output/{job_id}/{complex}/{mutation}/amber/{complex}.tleap.pdb',
    log:
        amber='output/{job_id}/{complex}/{mutation}/amber/amber.log',
        tleap='output/{job_id}/{complex}/{mutation}/amber/tleap.log'
    priority: 5
    shell:
        """
        pdb4amber -i {input} -o {output.amber} --dry --logfile {log.amber}  &&
        2_createleap.py --pdb {output.amber} \
                                --prmtop {output.prmtop} \
                                --inpcrd {output.inpcrd} \
                                --leap {output.leap} \
                                --tleappdb {output.tleappdb}
        tleap -f {output.leap} > {log.tleap}
        """


rule MD:
    input:
        prmtop='output/{job_id}/{complex}/{mutation}/amber/{complex}.prmtop',
        inpcrd='output/{job_id}/{complex}/{mutation}/amber/{complex}.inpcrd',
        default=ancient(md_settings),
        amber='output/{job_id}/{complex}/{mutation}/amber/amber.pdb',
        tleappdb='output/{job_id}/{complex}/{mutation}/amber/{complex}.tleap.pdb',
    output:
        traj='output/{job_id}/{complex}/{mutation}/{seed}/MD/trajectory.dcd',
        stats='output/{job_id}/{complex}/{mutation}/{seed}/MD/MDStats.csv',
        last='output/{job_id}/{complex}/{mutation}/{seed}/MD/frame_end.cif',
        pdb_last='output/{job_id}/{complex}/{mutation}/{seed}/MD/last.pdb',
        pdbx='output/{job_id}/{complex}/{mutation}/{seed}/MD/last.pdbx',
        state='output/{job_id}/{complex}/{mutation}/{seed}/MD/state.xml'
    params:
        directory=directory('output/{job_id}/{complex}/{mutation}/{seed}/MD/'),
        job_id=config['job']
    resources:
        gpu=1
    priority:
        3
    shell:
        """
        3_MD.py --topo {input.prmtop} \
                        --crd {input.inpcrd} \
                        --directory {params.directory} \
                        --traj {output.traj} \
                        --config {input.default} \
                        --simulation_name {wildcards.complex} \
                        --seed 23 \
                        --stats {output.stats} \
                        --job_id {params.job_id} \
                        --amber {input.amber} \
                        --tleap {input.tleappdb} \
                        --xml {output.state} \
                        --last {output.last} \
                        --pdb_last {output.pdb_last} \
                        --pdbx {output.pdbx} \
        """

rule centerTraj:
    input:
        topo='output/{job_id}/{complex}/{mutation}/{seed}/MD/frame_end.cif',
        traj='output/{job_id}/{complex}/{mutation}/{seed}/MD/trajectory.dcd',
        mapping='output/{job_id}/{complex}/{mutation}/amber/amber_renum.txt'
    #topo='output/{job_id}/{complex}/{mutation}/{seed}/MD/last.pdb',
    output:
        topo_center='output/{job_id}/{complex}/{mutation}/{seed}/MD/center/topo.pdb',
        traj_center='output/{job_id}/{complex}/{mutation}/{seed}/MD/center/traj.dcd',
        frame_pdb='output/{job_id}/{complex}/{mutation}/{seed}/frames/frame_1.pdb',
        ligand_csv='output/{job_id}/{complex}/{mutation}/{seed}/frames/lig/1.csv',
        receptor_csv='output/{job_id}/{complex}/{mutation}/{seed}/frames/rec/1.csv',
        dir=directory('output/{job_id}/{complex}/{mutation}/{seed}/frames/'),
        checkpoint='output/{job_id}/{complex}/{mutation}/{seed}/frames/.done'
    log:
        'output/{job_id}/{complex}/{mutation}/{seed}/MD/center.log'
    shell:
        """
        4_center.py    --topo {input.topo} \
                               --traj {input.traj} \
                               --topo_center {output.topo_center} \
                               --mapping {input.mapping} \
                               --frame_id {number_frames} \
                               --dir {output.dir} \
                               --traj_center {output.traj_center} > {log}

        touch {output.checkpoint}
        """

rule DescriptiveTrajAnalysis:
    input:
        topo='output/{job_id}/{complex}/{mutation}/{seed}/MD/center/topo.pdb',
        traj='output/{job_id}/{complex}/{mutation}/{seed}/MD/center/traj.dcd',
        mapping='output/{job_id}/{complex}/{mutation}/amber/amber_renum.txt',
        stats='output/{job_id}/{complex}/{mutation}/{seed}/MD/MDStats.csv'
    output:
        rmsf=report('output/{job_id}/{complex}/{mutation}/{seed}/analysis/RMSF.svg',caption="config/RMSF.rst", category="Trajectory", subcategory='RMSF'),
        rmsd=report('output/{job_id}/{complex}/{mutation}/{seed}/analysis/RMSD.svg',caption="config/RMSD.rst",category="Trajectory", subcategory='RMSD'),
        checkpoint='output/{job_id}/{complex}/{mutation}/{seed}/analysis/interactions/.checkpoint',
        ana_dir=directory('output/{job_id}/{complex}/{mutation}/{seed}/analysis/')
    params:
        dir = directory('output/{job_id}/{complex}/{mutation}/{seed}/'),
        number_frames = config['number_frames']
    shell:
        """
        6_ExplorativeTrajectoryAnalysis.py --topo {input.topo} \
                                                   --traj {input.traj} \
                                                   --stats {input.stats} \
                                                   --dir {params.dir} \
                                                   --mapping {input.mapping} \
                                                   --ana {output.ana_dir} \
                                                   --checkpoint {output.checkpoint}
        """




rule GlobalStatsMartin:
    input:
        dirs=expand('output/{job_id}/{complex}/{mutation}/{seed}/frames/',job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds),
        checkpoint=expand('output/{job_id}/{complex}/{mutation}/{seed}/frames/.done',job_id=config["job"], complex=complexes, mutation=mutations, seed=seeds)
    output:
        interactions = report('output/{job_id}/results/martin/interactions.csv',caption="config/RMSF.rst",category="Interaction Martin"),
        #energy_diff = report('output/{job_id}/results/martin/energy_diff.csv',caption="report/RMSF.rst",category="Energy Differences"),
        dir=directory('output/{job_id}/results/martin/')
    shell:
        """
        8_InteractionMartin.py  --dirs {input.dirs} \
                                        --output {output.dir} \
                                        --interactions {output.interactions}
        """


rule interactionFingerprint:
    input:
        traj = ('output/{job_id}/{complex}/{mutation}/{seed}/MD/trajectory.dcd'),
        topo = 'output/{job_id}/{complex}/{mutation}/{seed}/MD/frame_end.pdb',
        #topo='output/{job_id}/{complex}/{mutation}/{seed}/MD/center/topo.pdb',
        #traj='output/{job_id}/{complex}/{mutation}/{seed}/MD/center/traj.dcd',
        mapping='output/{job_id}/{complex}/{mutation}/amber/amber_renum.txt',
    output:
        'output/{job_id}/{complex}/{mutation}/{seed}/analysis/fingerprint.csv',
    shell:
        """
        5_interactionFingerprint.py --topo {input.topo} \
                                            --traj  {input.traj} \
                                            --residMapping {input.mapping} \
                                            --output {output}
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
