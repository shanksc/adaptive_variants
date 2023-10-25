import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import sys

def get_args():
    parser = argparse.ArgumentParser(description='plot manhattan plot', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tsv', type=str, required=True, help='bfs tsv')
    parser.add_argument('--contigs', type=str, help='txt containing qualified contigs', required=False, default=None)
    parser.add_argument('--cutoff', type=float, required=False, help='minimum dB', default=0)
    parser.add_argument('--xtx', required=False, help='plot XtX statistic', action='store_true')
    #parser.add_argument('--env', type=str, required=True, help='environmental variable to plot (ex latitude)')
    parser.add_argument('--o', type=str, required=True, help='output filename')

    return parser.parse_args()


def plot_manhattan(df, cutoff, var, out):
    #calculate decibans of BFs
    #df['dB'] = 10 * np.log10(df['BF(dB)'])
    df = df.sort_values(['#[1]CHROM', '[2]POS'])
    df = df[df[var] > cutoff]
    df.reset_index(inplace=True, drop=True)
    df['i'] = df.index

    #pal = ['#66c2a5', '#fc8d62']
    #sns.set_palette([pal])
    plot = sns.relplot(data=df, x='i', y=var, aspect=6, 
                        hue='#[1]CHROM', palette='colorblind', legend=None) 
    chrom_df=df.groupby('#[1]CHROM')['i'].median()
    plot.ax.set_xlabel('contig')
    plot.ax.set_xticks(chrom_df, labels=chrom_df.index)
    plot.ax.tick_params(axis='x', labelrotation=45)
    
    #add a horizontal line 
    #ax1, ax2 = plot.axes[0]

    #ax1.axhline(20, ls='--')
    #ax2.axhline(30, ls='--')
    if var == 'BF(dB)':
        plt.axhline(y=20, linestyle='--', color='grey', linewidth=1.5)
    plt.tight_layout()
    plt.savefig(out+'.png', bbox_inches='tight')
    '''
    with open(out+'_contig_names.txt', 'w') as f:
        for chrom in chrom_df.index:
    '''

def main():
    args = get_args()
    df = pd.read_csv(args.tsv, sep='\t', index_col=False)
    print(df.head())
    print(df.columns)

    if args.contigs is not None:
        contigs = []
        with open(args.contigs, 'r') as f:
            for line in f.readlines():
                contigs.append(line.strip())
        #contigs = [contig for contig in args.contigs]
        #for contig in args.contigs:
        contigs = set(contigs)
        df = df[df['#[1]CHROM'].isin(contigs)]

    var = 'XtXst' if args.xtx else 'BF(dB)'
    print(var)
    plot_manhattan(df, args.cutoff, var, args.o)
    #print(df.head())

    #print(df['#[1]CHROM'].unique())

if __name__ == '__main__':
    main()
