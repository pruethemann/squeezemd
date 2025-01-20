import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go

def plot_interaction_fingerprint(data):

    data['ligand_resid'] = data['ligand_resid'].astype(int)
    data = data[data.ligand_resid <= 122]



    fig = px.bar(data,
                x='ligand_resid',
                y='energy_mean',
                color='Interaction Type',
                error_y='energy_std'
                )
    
        # Set the bars to be grouped horizontally rather than stacked
    fig.update_layout(xaxis_title="Residue ID",
                    yaxis_title="Energy",
                    xaxis_tickangle=-90,
                    barmode='group')  # Group bars for different mutations side by side

    fig.update_layout(xaxis_title="Residue ID",
                        yaxis_title="Energy",
                        xaxis_tickangle=-90)
    
    # Save figure as an HTML file
    fig.write_html('test2.html')

df = pd.read_csv('posco_interactions.csv')

# clean up
#df = df[df.seed==695]
#df = df[df.frame==0]
#df = df[df.mutation=='WT']


#del df['frame']
#del df['seed']
#del df['mutation']

del df['Angle (a)']
del df['Unnamed: 0']
del df['Distance (r)']
del df['receptor_resid']

df = df[df.receptor_resname!='HOH']
df = df[df.ligand_resname!='HOH']



# ['Interaction Type', 'target' , 'lig', 'mutation', 'ligand_resid', 'seed']
print(df)
# Sum all energy contributions for every resid
resid_energy =  df.groupby(['Interaction Type', 'target' , 'lig', 'mutation','seed','frame', 'ligand_resid']).sum(numeric_only=True).reset_index()

print("resid energy")
print(resid_energy)

# Take the mean of all residue contributions for all frames
frame_mean =  resid_energy.groupby(['Interaction Type', 'target' , 'lig', 'mutation','seed', 'ligand_resid']).agg(energy_mean=('Energy (e)', 'mean'),energy_std=('Energy (e)', 'std')).reset_index()

print("frame mean energy")
print(frame_mean)

plot_interaction_fingerprint(frame_mean)

