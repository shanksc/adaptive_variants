import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
from scipy.stats import linregress
import sys

def get_args():
    parser = argparse.ArgumentParser(description='plot scatter between xtx and bfs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tsv', type=str, required=True, help='bfs/xtx tsv')
    parser.add_argument('--contig', type=str, help='contigs to plot', required=True)
    parser.add_argument('--cutoff', type=float, required=False, help='minimum dB', default=10)
    #parser.add_argument('--env', type=str, required=True, help='environmental variable to plot (ex latitude)')
    parser.add_argument('--o', type=str, required=True, help='output filename')

    return parser.parse_args()

def plot_scatter(df, cutoff, out):
    #calculate decibans of BFs
    #df['dB'] = 10 * np.log10(df['BF(dB)'])

    #df = df.sort_values(['#[1]CHROM', 'POSITION'])
    df = df[df['BF(dB)'] > cutoff]
    #df.reset_index(inplace=True, drop=True)
    #df['i'] = df.index
    #sns.regplot(data=df, x='BF(dB)', y='XtXst')
    bfs = df['BF(dB)'].to_numpy()
    xtx = df['XtXst'].to_numpy()

    res = linregress(bfs, xtx)

    #fig, ax = plt.subplots(figsize=(5,5))
    #plt.axis('equal')
    lim = max(max(bfs), max(xtx))
    print(lim)
    plt.ylim(0,lim)
    plt.xlim(0,lim)
    plt.scatter(bfs, xtx, alpha=.75)
    plt.plot(bfs, res.intercept + res.slope*bfs, 'r', label=f'R^2: {round(res.rvalue**2,4)}')
    plt.legend()
    plt.xlabel('BF (dB)')
    plt.ylabel('XtX')
    #pal = ['#66c2a5', '#fc8d62']
    #sns.set_palette([pal])
    '''plot = sns.relplot(data=df, x='i', y='BF(dB)', aspect=6, 
                        hue='#[1]CHROM', palette='colorblind', legend=None) 
    chrom_df=df.groupby('#[1]CHROM')['i'].median()
    plot.ax.set_xlabel('contig')
    plot.ax.set_xticks(chrom_df, labels=chrom_df.index)
    plot.ax.tick_params(axis='x', labelrotation=45)
    '''
    
    plt.tight_layout()
    plt.savefig(out+'.png', bbox_inches='tight')
    '''
    plt.clf()
    sns.scatterplot(data=df, x='Beta_is', y='XtXst')
    plt.tight_layout()
    plt.savefig(out+'_beta.png', bbox_inches='tight')
    '''

def main():
    args = get_args()

    #df = pd.read_csv(args.tsv, sep='\t', names=col_names, dtype=col_types)
    df = pd.read_csv(args.tsv, sep='\t', index_col=False)
    print(df.head())
    print(df.columns)
    #sys.exit()

    df = df[df['#[1]CHROM'] == args.contig]
    print(df.head())
    plot_scatter(df, args.cutoff, args.o)
    #print(df.head())

    #print(df['#[1]CHROM'].unique())

if __name__ == '__main__':
    main()
