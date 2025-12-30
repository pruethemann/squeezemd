from glob import glob
from os import path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Relevant data paths

P_20 = "/home/peter/caracara/Squeeze/P-20_SAR_meta_L55R/"
P_19 = "/home/peter/caracara/Squeeze/P-19_SAR_metadynamics/"

first = path.join(P_19, "metadynamics/C1s_Gigastasin/**/**/metadynamics/fes.dat")
second = path.join(P_20, "metadynamics/C1s_Gigastasin/**/**/metadynamics/fes.dat")
energy_data_1 = glob(first)
energy_data_2 = glob(second)

energy_data_1.extend(energy_data_2)

energy_data = energy_data_1

data = []

filter = ['WT', 'L55R']

for fes in energy_data:
    df = pd.read_csv(fes, sep='\s+', comment="#", header=None, names=['d1', 'F' ,'der_d1'])
    seed = fes.split('/')[-3]
    mutation = fes.split('/')[-4]
    id = fes.split('/')[-3]
    if 'P-20' in fes:
        id = "P-20"
    else:
        id = "P-19"

    df['seed'] = seed
    df['mutation'] = mutation
    df['id'] = id
    
    data.append(df)

data = pd.concat(data)

data['sim'] = data['id'] + '_' + data['seed'] + '_' + data['mutation']

data.to_parquet('energy.parquet')
print(data)
sims = data.sim.unique()

data_filtered = data[(data.mutation == 'WT')]
for sim in sims:
    sns.lineplot(data=data_filtered[data_filtered.sim==sim],
                x='d1',
                y='F',
                color="black")

plt.show()

data_filtered = data[(data.mutation == 'L55R')]
for sim in sims:
    sns.lineplot(data=data_filtered[data_filtered.sim==sim],
                x='d1',
                y='F',
                color="red")  

plt.show()
