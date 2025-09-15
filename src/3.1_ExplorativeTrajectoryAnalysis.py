#!/usr/bin/env python

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import pandas as pd
from Helper import remap_MDAnalysis
import MDAnalysis as mda
from MDAnalysis.analysis import rms
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

def find_series(sec_str, structure):
    """
    This function finds all consecutive series of 'H' in a given string and returns a list of tuples
    with the start and end positions of each series.

    :param sec_str: A string containing '-' and 'H'.
    :return: A list of tuples, where each tuple contains the start and end positions of a series of 'H'.
    """
    secondary_series = []
    serie_start = None
    
    for i, char in enumerate(sec_str):
        if char == structure and serie_start is None:
            serie_start = i  # Start of a new H series
        elif char != structure and serie_start is not None:
            secondary_series.append((serie_start+2, i))  # End of the current series
            serie_start = None
    
    # Handle case where string ends with an H series
    if serie_start is not None:
        secondary_series.append((serie_start, len(sec_str) - 1))
    
    return secondary_series

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
    visualize_RMSF(rmsf_data, secondary_structure_data, args.rmsf)

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
    return secondary


def visualize_RMSF(rmsf_data, secondary_structure_data, output_file):
    num_chains = len(rmsf_data)
    fig = make_subplots(rows=num_chains, cols=1, shared_xaxes=False, vertical_spacing=0.05)

    helix_color = 'rgba(0, 100, 250, 0.3)'
    sheet_color = 'rgba(250, 150, 0, 0.3)'
    legend = True
    legend_beta = True

    row = 1

    for chain, (resids, rmsf_values) in sorted(rmsf_data.items()):
        secondary_str = secondary_structure_data[chain]

        # Plot RMSF
        fig.add_trace(go.Scatter(x=resids, y=rmsf_values, mode='lines', name=f'Chain {chain} RMSF'),
                      row=row, col=1)
        
        ceil_rmsf = max(rmsf_values)
        floor_rmsf = min(rmsf_values)

        # Add secondary structure as filled areas instead of shapes
        beta_series = find_series(secondary_str, 'E')
        helix_series = find_series(secondary_str, 'H')

        start_residue = min(resids)

        # Draw secondary structure shades
        # 4. Add a blue shade in the background of the graph from x=5 to x=10           

        for helix in helix_series:
            # Second shade from x=15 to x=20
            # TODO: Check if the residue numbers for secondary structures are correct
            helix += start_residue
            fig.add_trace(go.Scatter(
                        x=[helix[0], helix[0], helix[1], helix[1]], y=[floor_rmsf, ceil_rmsf, ceil_rmsf, floor_rmsf], fill='toself',
                        fillcolor=helix_color, line=dict(color="rgba(0, 0, 0, 0)"),
                        name='Alpha', legendgroup='Alpha', showlegend=legend),row=row, col=1)
            # Only show one legend (first)
            legend=False

        for beta in beta_series:
            beta += start_residue
            fig.add_trace(go.Scatter(
                        x=[beta[0], beta[0], beta[1], beta[1]], y=[floor_rmsf, ceil_rmsf, ceil_rmsf, floor_rmsf], fill='toself',
                        fillcolor=sheet_color, line=dict(color="rgba(0, 0, 0, 0)"),
                        name='Beta', legendgroup='Beta', showlegend=legend_beta),row=row, col=1)
            # Only show one legend (first)
            legend_beta=False

        # Set the x-axis title for each subplot
        fig.update_xaxes(title_text=f"Residue IDs (Chain {chain})", row=row, col=1)

        # Update layout to ensure the shading fits the entire y-axis range and show legend
        fig.update_layout(title="Toy Example with Shaded Area", xaxis_title="Residue number", yaxis_title="RMSF", 
                  showlegend=True)
        
        x_axis_range = [min(resids), max(resids)]

        # Inside your function after plotting everything but before saving or showing the figure:
        fig.update_xaxes(range=x_axis_range, row=row, col=1)

        row += 1

    # Update layout with proper size and a clear title
    fig.update_layout(height=300*num_chains, width=800, title_text="RMSF and Secondary Structure per Chain",
                      legend=dict(yanchor="top", y=0.99, xanchor="left", x=0.01))
    
    # Set the x-axis range to the min and max of all residue IDs


    # Save the figure to the specified output file
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

    df_melted.to_csv(args.rmsd_raw)
    plt.savefig(args.rmsd)
    plt.close()


def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--topo', required=False, help='Topo file', default='frame_end.cif')
    parser.add_argument('--traj', required=False, help='Trajectory', default='traj_center.dcd')
    parser.add_argument('--stats', required=False, default='MDStats.csv')

    # Output
    parser.add_argument('--rmsf', required=False, default='results/rmsf.html', help='')
    parser.add_argument('--bfactors', required=False, help='', default='results/bfactors.pdb')
    parser.add_argument('--rmsd', required=False, help='', default='results/rmsd.png')
    parser.add_argument('--rmsd_raw', required=False, help='', default='results/rmsd.png')
    parser.add_argument('--fig_stats', required=False, help='', default='results/stats.png')

    return parser.parse_args()

# Example of running the function
if __name__ == '__main__':
    args = parse_arguments()

    # Import Trajectory
    topo = app.PDBxFile(args.topo)
    u = mda.Universe(topo, args.traj, in_memory=False)
    u = remap_MDAnalysis(u, topo)

    traj_length = len(u.trajectory)
    #print(f'Number of frames: {traj_length}')

    calculate_RMSF_and_secondary_structure(u, args)

    calculate_RMSD(u, args)

    # Visualize Energies, T, ...
    visualize_MDStats(args.stats, args.fig_stats)
