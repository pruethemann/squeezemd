#!/usr/bin/env python
import prolif as plf
from prolif.plotting.network import LigNetwork
import MDAnalysis as mda
import sys
import argparse
import pandas as pd
import networkx as nx
from pyvis.network import Network
from tqdm.auto import tqdm
from matplotlib import cm, colors
from IPython.display import IFrame
#import py3Dmol

# deactivate import in
# /home/peter/miniconda3/envs/flow/lib/python3.10/sqlite3/__init__.py
#import sqlite3

"""
Visualisation tutorial:
https://prolif.readthedocs.io/en/latest/notebooks/visualisation.html


python3 5_interactionFingerprint.py --topo output/3_HellostasinChimera/C1s_BD001/A81S/amber/C1s_BD001.prmtop \
                                   --traj  output/3_HellostasinChimera/C1s_BD001/A81S/767/MD/trajectory.dcd \
                                   --residMapping output/3_HellostasinChimera/C1s_BD001/A81S/amber/amber_renum.txt \
                                   --output  fingerprint_prmtop.csv

python3 5_interactionFingerprint.py --topo output/3_HellostasinChimera/C1s_BD001/D60S/amber/C1s_BD001.prmtop \
                                   --traj  output/3_HellostasinChimera/C1s_BD001/D60S/767/MD/trajectory.dcd \
                                   --residMapping output/3_HellostasinChimera/C1s_BD001/A81S/amber/amber_renum.txt \
                                   --output  fingerprint_D60S.csv

"""


def import_renum(file_csv):
    """
    Import amber remapping file.
    :param file_csv:
    :return:
    """
    print(file_csv)
    renum = pd.read_csv(file_csv,
                        delim_whitespace=True,
                        names=['resname', 'chainID', 'resid', 'resname amber', 'resid amber'])
    return renum


def chain2resid(renum):
    """
    Find start and end of chain A
    :param renum:
    :return:
    """

    chains = {}

    for chainID in renum.chainID.unique():
        chain = renum[renum.chainID == chainID]

        resID_min = chain['resid amber'].min()
        resID_max = chain['resid amber'].max()

        chains[chainID] = (resID_min, resID_max)

    return chains

def make_graph(values, df=None,
               node_color=["#FFB2AC", "#ACD0FF"], node_shape="dot",
               edge_color="#a9a9a9", width_multiplier=1):
    """Convert a pandas DataFrame to a NetworkX object

    Parameters
    ----------
    values : pandas.Series
        Series with 'ligand' and 'protein' levels, and a unique value for
        each lig-prot residue pair that will be used to set the width and weigth
        of each edge. For example:

            ligand  protein
            LIG1.G  ALA216.A    0.66
                    ALA343.B    0.10

    df : pandas.DataFrame
        DataFrame obtained from the fp.to_dataframe() method
        Used to label each edge with the type of interaction

    node_color : list
        Colors for the ligand and protein residues, respectively

    node_shape : str
        One of ellipse, circle, database, box, text or image, circularImage,
        diamond, dot, star, triangle, triangleDown, square, icon.

    edge_color : str
        Color of the edge between nodes

    width_multiplier : int or float
        Each edge's width is defined as `width_multiplier * value`
    """
    lig_res = values.index.get_level_values("ligand").unique().tolist()
    prot_res = values.index.get_level_values("protein").unique().tolist()

    G = nx.Graph()
    # add nodes
    # https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.network.Network.add_node
    for res in lig_res:
        G.add_node(res, title=res, shape=node_shape,
                   color=node_color[0], dtype="ligand")
    for res in prot_res:
        G.add_node(res, title=res, shape=node_shape,
                   color=node_color[1], dtype="protein")

    for resids, value in values.items():
        label = "{} - {}<br>{}".format(*resids, "<br>".join([f"{k}: {v}"
                                       for k, v in (df.xs(resids,
                                                          level=["ligand", "protein"],
                                                          axis=1)
                                                      .sum()
                                                      .to_dict()
                                                      .items())]))
        # https://pyvis.readthedocs.io/en/latest/documentation.html#pyvis.network.Network.add_edge
        G.add_edge(*resids, title=label, color=edge_color,
                   weight=value, width=value*width_multiplier)

    return G

def create_network(df):
    data = (df.groupby(level=["ligand", "protein"], axis=1, sort=False)
            .sum()
            .astype(bool)
            .mean())

    G = make_graph(data, df, width_multiplier=5)

    # color each node based on its degree
    max_nbr = len(max(G.adj.values(), key=lambda x: len(x)))
    palette = cm.get_cmap('YlGnBu', max_nbr)
    for n, d in G.nodes(data=True):
        n_neighbors = len(G.adj[n])
        d["color"] = colors.to_hex(palette(n_neighbors / max_nbr))

    # convert to pyvis network
    net = Network(width=640, height=500, notebook=True, heading="")
    net.from_nx(G)

    # use specific layout
    layout = nx.circular_layout(G)
    for node in net.nodes:
        node["x"] = layout[node["id"]][0] * 1000
        node["y"] = layout[node["id"]][1] * 1000
    net.toggle_physics(False)

    net.write_html("residue-network_graph.html")
    IFrame("residue-network_graph.html", width=650, height=510)

def fancy_network(df):
    data = (df.groupby(level=["ligand", "protein"], axis=1, sort=False)
            .sum()
            .astype(bool)
            .mean())

    G = make_graph(data, df, width_multiplier=8)

    # color each node based on its degree
    max_nbr = len(max(G.adj.values(), key=lambda x: len(x)))
    blues = cm.get_cmap('Blues', max_nbr)
    reds = cm.get_cmap('Reds', max_nbr)
    for n, d in G.nodes(data=True):
        n_neighbors = len(G.adj[n])
        # show TM3 in red and the rest of the protein in blue
        palette = reds if d["dtype"] == "ligand" else blues
        d["color"] = colors.to_hex(palette(n_neighbors / max_nbr))

    # convert to pyvis network
    net = Network(width=640, height=500, notebook=True, heading="")
    net.from_nx(G)
    net.write_html("prot-prot_graph.html")
    IFrame("prot-prot_graph.html", width=650, height=510)


def create_interactionFingerprint(args):

    # Handle the renaming of residues from amber
    resid_mapping = import_renum(args.residMapping)
    chains = chain2resid(resid_mapping)

    # Import trajectory
    u = mda.Universe(args.topo, args.traj, in_memory=False)

    # Define ligand (gigastasin) and receptor
    ligand = u.select_atoms(f"resid {chains['I'][0]}:{chains['I'][1]}")
    protein = u.select_atoms(f"resid {chains['A'][0]}:{chains['A'][1]}")

    # Run interaction fingerprint analysis
    fp = plf.Fingerprint(["HBDonor", "HBAcceptor", "PiStacking", "PiCation", "CationPi", "Anionic", "Cationic"])
    fp.run(u.trajectory, ligand, protein) # u.trajectory[::10]

    # Export interactions
    interactions_df = fp.to_dataframe()
    interactions_df.to_csv(args.output)

    interactions_df.to_csv('delete.csv')

    create_network(interactions_df)

    fancy_network(interactions_df)

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('--topo', required=False,help='Topo file', default='')
    parser.add_argument('--traj', required=False,help='Trajectory file', default='')
    parser.add_argument('--residMapping', required=False, help='Mapping resid to amber resid', default='')
    parser.add_argument('--output', required=False, help='Mapping resid to amber resid', default='')

    args = parser.parse_args()

    print(args)
    create_interactionFingerprint(args)


"""
Code Snippets
    Interaction graph
    Only Version 2.0

    net = fp.to_ligplot(
        lmol,
        # replace with `kind="frame", frame=0,` for the other depiction
        kind="aggregate",
        threshold=0.3,
        rotation=270,
    )
    net.display()
"""
