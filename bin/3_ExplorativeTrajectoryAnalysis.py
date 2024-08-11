#!/usr/bin/env python

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import rms
from Helper import remap_MDAnalysis
import openmm.app as app
from MDAnalysis.analysis.dssp import DSSP, translate
import argparse
import matplotlib.pyplot as plt
import seaborn as sns

def visualize_MDStats(stats_file, output_graph):
    data = pd.read_csv(stats_file, sep='\t')

    data['time (ns)'] = data['Time (ps)'] / 1000

    # Total Energy
    plt.subplot(2,2,1)
    sns.lineplot(data=data,
                 x='time (ns)',
                 y='Total Energy (kJ/mole)')

    plt.title("Total Energy")

    # Potential Energy
    plt.subplot(2,2,2)
    sns.lineplot(data=data,
                 x='time (ns)',
                 y='Potential Energy (kJ/mole)')

    plt.title("Potential Energy (kJ/mole)")

    # Temperature
    plt.subplot(2,2,3)
    sns.lineplot(data=data,
                 x='time (ns)',
                 y='Temperature (K)')

    plt.title("Temperature (K) Mean: " + str(data['Temperature (K)'].mean().round(2)))

    # Volume
    plt.subplot(2,2,4)
    sns.lineplot(data=data,
                 x='time (ns)',
                 y='Box Volume (nm^3)')

    plt.title("Box Volume (nm^3) Mean: " + str(data['Box Volume (nm^3)'].mean().round(2)))

    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
    plt.savefig(output_graph)


def calculate_RMSF_and_secondary_structure(u: mda.Universe, args):
    # Extract all unique chainIDs from the Universe object
    chains = list(set(atom.chainID for atom in u.atoms))
    # Exclude numeric values which correspond to salts and solvents
    chains = [item for item in chains if not item.isdigit()]
    
    rmsf_data = {}
    secondary_structure_data = {}

    for chain in chains:
        c_alphas = u.select_atoms(f'chainID {chain} and name CA')
        R = rms.RMSF(c_alphas).run()

        # Store RMSF and secondary structure data
        rmsf_data[chain] = (c_alphas.resids, R.results.rmsf)
        secondary_str = predict_secondary_structure(u, chain)
        secondary_structure_data[chain] = secondary_str

    # Visualize combined RMSF and secondary structure for all chains
    visualize_combined_rmsf_and_secondary_structure(rmsf_data, secondary_structure_data, args.rmsf)

    # Calculate bfactors
    c_alphas = u.select_atoms('protein and name CA')
    R = mda.analysis.rms.RMSF(c_alphas).run()
    calculate_bfactors(R)


def calculate_bfactors(R):
    u.add_TopologyAttr('tempfactors')  # add empty attribute for all atoms
    protein = u.select_atoms('protein')  # select protein atoms
    for residue, r_value in zip(protein.residues, R.results.rmsf):
        residue.atoms.tempfactors = r_value

    u.atoms.write(args.bfactors)

def predict_secondary_structure(u: mda.Universe, chainID: str):
    chain = u.select_atoms(f'chainID {chainID}')
    dssp_analysis = DSSP(chain).run()
    mean_secondary_structure = translate(dssp_analysis.results.dssp_ndarray.mean(axis=0))
    secondary = ''.join(mean_secondary_structure)
    print(secondary)
    return secondary


def visualize_combined_rmsf_and_secondary_structure(rmsf_data, secondary_structure_data, output_file):
    num_chains = len(rmsf_data)
    fig = make_subplots(rows=num_chains, cols=1, shared_xaxes=False, vertical_spacing=0.05)

    helix_color = 'rgba(0, 100, 250, 0.3)'
    sheet_color = 'rgba(250, 150, 0, 0.3)'

    row = 1
    for chain, (resids, rmsf_values) in sorted(rmsf_data.items()):
        secondary_str = secondary_structure_data[chain]

        # Plot RMSF
        fig.add_trace(go.Scatter(x=resids, y=rmsf_values, mode='lines', name=f'Chain {chain} RMSF'),
                      row=row, col=1)

        # Add secondary structure as filled areas instead of shapes
        helix_x, helix_y = [], []
        sheet_x, sheet_y = [], []
        y_max = max(rmsf_values)  # maximum value on y-axis
        for i, ss in enumerate(secondary_str):
            if ss == 'H':  # Helix
                helix_x.extend([resids[i], resids[i] + 1, resids[i] + 1, resids[i], resids[i]])  # Loop back to start
                helix_y.extend([0, 0, y_max, y_max, 0])  # Loop back to start
            elif ss == 'E':  # Sheet
                sheet_x.extend([resids[i], resids[i] + 1, resids[i] + 1, resids[i], resids[i]])  # Loop back to start
                sheet_y.extend([0, 0, y_max, y_max, 0])  # Loop back to start

        # Add filled areas for helices
        if helix_x:
            fig.add_trace(go.Scatter(x=helix_x, y=helix_y, fill='toself', mode='lines', line=dict(color=helix_color),
                                     name='Helix', legendgroup='Helix', showlegend=(row == 1)),
                          row=row, col=1)

        # Add filled areas for sheets
        if sheet_x:
            fig.add_trace(go.Scatter(x=sheet_x, y=sheet_y, fill='toself', mode='lines', line=dict(color=sheet_color),
                                     name='Sheet', legendgroup='Sheet', showlegend=(row == 1)),
                          row=row, col=1)

        # Set the x-axis title for each subplot
        fig.update_xaxes(title_text=f"Residue IDs (Chain {chain})", row=row, col=1)

        row += 1

    # Update layout with proper size and a clear title
    fig.update_layout(height=300*num_chains, width=800, title_text="RMSF and Secondary Structure per Chain",
                      legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))

    # Save the figure to the specified output file
    fig.write_image('output_rmsf_secondary_structure_corrected.png')  # Save as an image file
    fig.write_html(output_file)  # Save as an HTML file for interactive viewing
    #fig.show()

def calculate_RMSD(u: mda.Universe, args):
    """
    Calculate RMSD of receptor and ligand
    :param u:
    :param output:
    :return:
    """

    print("Init RMSD analysis")
    # CHAINIDENTIFICAITON
    ligand = u.select_atoms('chainID A')
    receptor = u.select_atoms('chainID B')

    # 3. Compute RMSD for receptor and ligand
    RMSD_ligand = rms.RMSD(ligand, ref_frame=0).run()
    RMSD_receptor = rms.RMSD(receptor, ref_frame=0).run()

    # 4. Save the data in a dataframe
    data = {
        'Time (ns)': RMSD_ligand.times,
        'Ligand': RMSD_ligand.results.rmsd[:, 2],  # Column 2 contains the RMSD values
        'Receptor': RMSD_receptor.results.rmsd[:, 2],  # Column 2 contains the RMSD values
    }
    df = pd.DataFrame(data)

    # Melt the dataframe for seaborn plotting
    df_melted = df.melt(id_vars=["Time (ns)"], var_name="Molecule", value_name="RMSD")

    # 5. Plot the data with seaborn
    sns.lineplot(data=df_melted, x="Time (ns)", y="RMSD", hue="Molecule")
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (Å)')
    plt.title('RMSD over Time')
    plt.legend(title='Molecule')
    plt.tight_layout()

    df_melted.to_csv(args.rmsd[-4:] + '.csv')
    plt.savefig(args.rmsd)
    plt.close()


def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=True, help='Topo file', default='output/WT-simple/topo.pdb')
    parser.add_argument('--traj', required=True, help='Trajectory', default='output/WT-simple/traj_150.dcd')
    parser.add_argument('--stats', required=True, help='MD Stats file')

    # Output
    parser.add_argument('--rmsf', required=False, default='results/rmsf.html', help='')
    parser.add_argument('--bfactors', required=False, help='', default='results/bfactors.pdb')
    parser.add_argument('--rmsd', required=False, help='', default='results/rmsd.png')
    parser.add_argument('--fig_stats', required=False, help='', default='results/stats.png')

    return parser.parse_args()

# Example of running the function
if __name__ == '__main__':
    args = parse_arguments()

    # Import Trajectory
    topo = app.PDBxFile(args.topo)
    u = mda.Universe(topo, args.traj, in_memory=False)

    traj_length = len(u.trajectory)
    print(f'Number of frames: {traj_length}')

    calculate_RMSF_and_secondary_structure(u, args)

    calculate_RMSD(u, args)

    # Visualize Energies, T, ...
    visualize_MDStats(args.stats, args.fig_stats)
