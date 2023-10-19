import pandas as pd
import argparse
import matplotlib
matplotlib.use('Agg')
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns


def get_args():
    parser = argparse.ArgumentParser(description='plot manhattan plot', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--tsv', type=str, required=True, help='bfs tsv')
    #maybe we want to do this earlier when we create the tsv
    parser.add_argument('--contaminent', type=str, required=False, help='remove contaminent chromosome(s)', nargs='*', default=None)
    parser.add_argument('--cutoff', type=float, required=False, help='minimum dB', default=0)
    parser.add_argument('--env', type=str, required=True, help='environmental variable to plot (ex latitude)')
    parser.add_argument('--o', type=str, required=True, help='output filename')

    return parser.parse_args()


def plot_manhattan(df, env, cutoff, out):
    #calculate decibans of BFs
    df['dB'] = 10 * np.log10(df[env])
    df = df.sort_values(['chromosome', 'position'])
    df = df[df['dB'] > cutoff]
    df.reset_index(inplace=True, drop=True)
    df['i'] = df.index

    #pal = ['#66c2a5', '#fc8d62']
    #sns.set_palette([pal])
    plot = sns.relplot(data=df, x='i', y='dB', aspect=6, 
                        hue='chromosome', palette='colorblind', legend=None) 
    chrom_df=df.groupby('chromosome')['i'].median()
    plot.ax.set_xlabel('contig')
    plot.ax.set_xticks(chrom_df, labels=chrom_df.index)
    plot.ax.tick_params(axis='x', labelrotation=45)
    #was chrom_df.index
    #plot.ax.set_xticklabels(chrom_df.index, )
    #plot.ax.set_xticklabels([i for i in range(1, len(chrom_df.index)+1)])
    #fig = plot.get_fig()
    #fig.savefig(out+'.png')
    #print(chrom_df.index)
    plt.savefig(out+'.png', bbox_inches='tight')
    '''
    with open(out+'_contig_names.txt', 'w') as f:
        for chrom in chrom_df.index:
    '''

def main():
    args = get_args()

    df = pd.read_csv(args.tsv, sep='\t')

    if args.contaminent is not None:
        for chrom in args.contaminent:
            df = df[df['chromosome'] != chrom]

    plot_manhattan(df, args.env, args.cutoff, args.o)
    #print(df.head())

    #print(df['chromosome'].unique())

if __name__ == '__main__':
    main()
