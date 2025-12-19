#!/usr/bin/env python

import pandas as pd
import plotly.express as px
import plotly.subplots as sp
import plotly.graph_objects as go
import argparse

def import_data(fingerprints):
    """
        Imports and processes fingerprint data from a list of file paths.

        Parameters:
        - fingerprints: A list of file paths to fingerprint data files.

        Returns:
        - A combined DataFrame containing processed fingerprint data from all files.

        Note: This function relies on 'data_engineering' to process individual DataFrames.
        """

    ## Initialize an empty list to store processed data
    combined_data = []

    for fp_path in fingerprints:

        # Import data
        try:
            fp = pd.read_parquet(fp_path)
            fp = data_engineering(fp, args.n_frames)
        except FileNotFoundError:
            print("Error with import from: ", fp_path)
            continue

        # Determine metrics lables
        fp['name'] = fp.attrs['complex']
        fp['target'] = fp.attrs['target']
        fp['lig'] = fp.attrs['ligand']
        fp['mutation'] = fp.attrs['mutation']
        fp['seed'] = fp.attrs['seed']

        # Save data
        combined_data.append(fp)

    # Merge all data together
    data = pd.concat(combined_data)
    data.attrs['n_residues'] = fp.attrs['n_residues_ligand']
    return data

def data_engineering(data, n_frames):
    """
    Processes individual fingerprint data.

    Parameters:
    - fp_data: DataFrame containing fingerprint data.
    - n_frames: An integer indicating the number of frames to process.

    Returns:
    - Processed DataFrame.
    """

    # Aggregate all interactions
    data_agg = data.sum()
    data_agg = data_agg.reset_index()

    # rename columns
    data_agg.columns = ['ligand', 'target', 'interaction', 'sum']
    data_agg['sum']  /= n_frames

    # Group interactions
    interaction_map = {
        'Cationic': 'salt bridge',
        'Anionic': 'salt bridge',
        'HBAcceptor': 'H bonds',
        'HBDonor': 'H bonds',
        'PiStacking': 'PiStacking',
        "PiCation": 'Cation-Pi',
        "CationPi": 'Cation-Pi',
        "Hydrophobic": "Hydrophobic"
    }

    data_agg['interaction_type'] = data_agg['interaction'].map(interaction_map)

    # extract resids
    data_agg['resid'] = data_agg['ligand'].str.extract('(\d+)').astype(int)
    return data_agg

def create_fig(fp_df, fig_path):

    # Group by 'interaction_type', 'mutation', and 'resid', and calculate mean and standard deviation of the 'sum'
    df_grouped = fp_df.groupby(['interaction_type', 'mutation', 'resid']).agg(
        mean_sum=('sum', 'mean'),
        std_sum=('sum', 'std')
    ).reset_index()

    # Determine the minimum and maximum residue across all data
    min_resid = fp_df['resid'].min()
    max_resid = fp_df['resid'].max()

    # Define a color map for mutations
    unique_mutations = df_grouped['mutation'].unique()
    color_map = px.colors.qualitative.Plotly[:len(unique_mutations)]  # Use Plotly's color scheme
    mutation_color_mapping = dict(zip(unique_mutations, color_map))

    # Create subplots for each 'interaction_type' with consistent colors for mutations
    interaction_types = df_grouped['interaction_type'].unique()
    fig = sp.make_subplots(rows=len(interaction_types), cols=1, subplot_titles=interaction_types)

    # Add bar plots to the corresponding subplots
    for i, interaction_type in enumerate(interaction_types):
        filtered_data = df_grouped[df_grouped['interaction_type'] == interaction_type]
        
        for mutation in filtered_data['mutation'].unique():
            mutation_data = filtered_data[filtered_data['mutation'] == mutation]
            
            # Create a bar plot for each mutation, ensuring color consistency
            fig.add_trace(
                go.Bar(
                    x=mutation_data['resid'],
                    y=mutation_data['mean_sum'],
                    error_y=dict(type='data', array=mutation_data['std_sum']),
                    name=mutation,
                    legendgroup=mutation,  # Group by mutation name for synchronized legend
                    marker_color=mutation_color_mapping[mutation],  # Ensure consistent color for the mutation
                    showlegend=(i == 0)  # Show legend only in the first subplot
                ),
                row=i+1, col=1
            )

        # Update x-axis to cover the full range of residues for consistency
        fig.update_xaxes(range=[min_resid, max_resid], tickvals=list(range(min_resid, max_resid + 1)), row=i+1, col=1)

    # Update the layout of the figure to have a shared legend and consistent color mapping
    fig.update_layout(
        height=300 * len(interaction_types),  # Adjust height based on number of subplots
        title_text="Interaction Types vs Residue",
        showlegend=True,
        barmode='group'
    )

    # Save the updated figure as an html file
    fig.write_html(fig_path)


def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--fingerprints', nargs='+', help="List of fingerprint parquet files", required=True, default=['fp1.parquet', 'fp2.parquet'])
    parser.add_argument('--n_frames', help='How many frames to be analysed. Only the last n frames will be analyzed. Default 100', type=int, required=False, default=100)

    # Output
    parser.add_argument('--interactions', help="Joined fingerprint csv file path", required=False, default='fingerprint.csv')
    parser.add_argument('--figure', help="Fingerprints figure", required=False, default='fingerprint.html')
    return parser.parse_args()

if __name__ == '__main__':
    args = parse_arguments()

    # Import all fingerprints data
    fingerprints = import_data(args.fingerprints)

    create_fig(fingerprints, args.figure)

    # Export data
    fingerprints.to_parquet(args.interactions)