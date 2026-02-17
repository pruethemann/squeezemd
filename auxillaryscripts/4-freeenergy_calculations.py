from glob import glob
from os import path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# Relevant data paths

P_22 = "/home/peter/caracara/Squeeze/P-25_SAR_meta_L55R_contacts/"

first = path.join(P_22, "metadynamics/C1s_Gigastasin/**/**/metadynamics/fes.dat")

energy_data = glob(first)

data = []

filter = ['WT', 'L55R']

for fes in energy_data:
    names_1CV = ['d1', 'F' ,'der_d1']
    columns_2CV = ['d1', 'c1','F' ,'der_d1', 'der_c1']
    df = pd.read_csv(fes, sep='\s+', comment="#", header=None, names=columns_2CV)
    seed = fes.split('/')[-3]
    mutation = fes.split('/')[-4]
    id = fes.split('/')[-3]

    df['seed'] = seed
    df['mutation'] = mutation
    df['id'] = 'P_22'
    
    data.append(df)

data = pd.concat(data)

data['sim'] = data['id'] + '_' + data['seed'] + '_' + data['mutation']

data.to_parquet('energy.parquet')
print(data)
sims = data.sim.unique()

print(data)

data_filtered = data[(data.mutation == 'WT')]
for sim in sims:
    sns.lineplot(data=data_filtered[data_filtered.sim==sim],
                x='d1',
                y='F')
    
    sns.lineplot(data=data_filtered[data_filtered.sim==sim],
                x='c1',
                y='F')

plt.savefig("WT.png")
plt.show()

data_filtered = data[(data.mutation == 'L55R')]
for sim in sims:
    sns.lineplot(data=data_filtered[data_filtered.sim==sim],
                x='d1',
                y='F')

plt.savefig("L55R.png")
plt.show()
