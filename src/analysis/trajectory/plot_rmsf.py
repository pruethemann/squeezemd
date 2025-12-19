#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse


def parse_arguments():
    parser = argparse.ArgumentParser()

    # Input
    parser.add_argument('--input', required=False)

    # Output
    parser.add_argument('--output', required=False, default='rmsf.svg', help='')

    return parser.parse_args()

# Example of running the function
if __name__ == '__main__':

    args = parse_arguments()

    rmsf_df = pd.read_parquet(args.input)

    sns.lineplot(data=rmsf_df,
                x='resid',
                y='rmsf',
                errorbar='sd')
    
    plt.savefig(args.output)
    plt.close()
