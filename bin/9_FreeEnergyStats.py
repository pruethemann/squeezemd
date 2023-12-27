#!/usr/bin/env python
import argparse
import ast
import pandas as pd
import os
import re
import seaborn as sns
import matplotlib.pyplot as plt


def convert_line_to_numbers(line):
    line = line.strip()[11:]
    line = re.sub(' +', ' ', line).strip().split(' ')
    line = [float(number) for number in line]
    return line

def extract_data(ene_path):
    try:
        with open(ene_path, 'r') as f:
            for line in f.readlines():
                if line.startswith('DELTA TOTAL'):
                    ene_values = convert_line_to_numbers(line)
                    ene_values = {'free_ene_diff': ene_values[0],
                                  'SD': ene_values[1],
                                  'SEM': ene_values[2]}

                    return ene_values
    except FileNotFoundError:
        return {'free_ene_diff': None,
                                  'SD': None,
                                  'SEM': None}





if __name__ == '__main__':
    # Parse Arguments
    parser = argparse.ArgumentParser()

    parser.add_argument('--input', required=False, default='output/demo/simulations.csv')

    parser.add_argument('--csv', required=False, default='output/demo/results/free_energy.csv')
    parser.add_argument('--fig', required=False, default='output/demo/results/free_energy.png')
    parser.add_argument('--job_id', required=False, default='demo')


    args = parser.parse_args()


    sims = pd.read_csv(args.input, index_col='sim_id', converters={"mutations": ast.literal_eval})

    # TODO move to setup
    # define paths free energy stats
    sims['free_ene_path'] = sims['path'] + '/analysis/freeEnergy.dat'

    # Extract simulation data into different columns
    sims['delta free energy'] = sims['free_ene_path'].apply(extract_data)
    df2 = pd.json_normalize(sims['delta free energy']).set_index(sims.index)
    sims = pd.concat([sims, df2], axis='columns')

    sims.to_csv(f'output/{args.job_id}/results/simulation_data.csv')

    # Aggregate over Ligand, receptor and mutation
    #sims_agg = sims.groupby(['ligand', 'target', 'mutation_all']).mean()
    sims_agg = sims[['name', 'free_ene_diff']].groupby(['name']).agg(['mean', 'std'])

    # Visualisation
    print(sims_agg)
    sims_agg.to_csv(args.csv)

    print(sims['free_ene_diff'])

    sims_agg.reset_index(inplace=True)
    sns.boxplot(data=sims, x='name',
                y='free_ene_diff',
                #order=["C1s_BD001_T69R", "C1s_BD001_Q45K", "C1s_BD001_", "C1s_BD001_D18W", "C1s_BD001_D18E_G36R", "C1s_BD001_R65A"]
                )


#  order=sims.sort_values('free_ene_diff').name)

    plt.ylabel('Free energy kcal/mol')
    plt.xlabel('Simulations')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(args.fig, dpi=600)
    plt.show()
    plt.close()


    sns.barplot(data=sims, x='name',
                y='free_ene_diff',
                order=["C1s_BD001_T69R", "C1s_BD001_Q45K", "C1s_BD001_", "C1s_BD001_D18W", "C1s_BD001_D18E_G36R", "C1s_BD001_R65A"],
                errorbar='sd',
                capsize=.3
                )

    plt.ylabel('Free energy kcal/mol')
    plt.xlabel('Simulations')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f'output/{args.job_id}/results/barplot.png', dpi=600)
    plt.show()

# TODO: determine missing values to give feedback about progress!
incompleted_sims = sims['free_ene_diff'].isna().sum()
print(incompleted_sims, 1- incompleted_sims/len(sims))


"""
python3 9_FreeEnergyStats.py --job_id 231222_Trp17 --input output/231222_Trp17/simulations.csv --csv free_energy_trp17.csv --fig fee_energy.png

"""
