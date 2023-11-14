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
    parser.add_argument('--no_p', help='do not plot p value', action='store_true')
    parser.add_argument('--o', type=str, required=True, help='output filename')

    return parser.parse_args()

def plot_scatter(df, cutoff, no_p, out):
    #calculate decibans of BFs
    #df['dB'] = 10 * np.log10(df['BF(dB)'])

    #df = df.sort_values(['#[1]CHROM', 'POSITION'])
    df = df[df['BF(dB)'] > cutoff]
    #df.reset_index(inplace=True, drop=True)
    #df['i'] = df.index
    #sns.regplot(data=df, x='BF(dB)', y='XtXst')
    bfs = df['BF(dB)'].to_numpy()
    xtx = df['XtXst'].to_numpy()
    pvals = df['log10(1/pval)'].to_numpy()

    res = linregress(bfs, xtx)

    #fig, ax = plt.subplots(figsize=(5,5))
    #plt.axis('equal')
    lim = max(max(bfs), max(xtx)) + 1
    print(lim)
    plt.ylim(0,lim)
    plt.xlim(0,lim)

    if no_p:
        plt.scatter(x=bfs, y=xtx, alpha=.65)
    else:
        #plot p vals
        plt.scatter(x=bfs, y=xtx, c=pvals, cmap=plt.cm.viridis)
        cbar = plt.colorbar()
        cbar.ax.set_ylabel(r"$10\log_{10}( \frac{1}{pval} )$", rotation=270, labelpad=20)
    
    
   
    
    plt.plot(bfs, res.intercept + res.slope*bfs, 'r', label=f'R^2: {round(res.rvalue**2,4)}')
    plt.legend()
    plt.xlabel('BF (dB)')
    plt.ylabel('XtX')
    plt.savefig(out+'.png', dpi=800, bbox_inches='tight')
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
    plot_scatter(df, args.cutoff, args.no_p, args.o)
    #print(df.head())

    #print(df['#[1]CHROM'].unique())

if __name__ == '__main__':
    main()
