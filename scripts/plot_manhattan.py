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
    #parser.add_argument('--env', type=str, required=True, help='environmental variable to plot (ex latitude)')
    parser.add_argument('--o', type=str, required=True, help='output filename')

    return parser.parse_args()


def plot_manhattan(df, cutoff, out):
    #calculate decibans of BFs
    #df['dB'] = 10 * np.log10(df['BF(dB)'])

    df = df.sort_values(['CHROMOSOME', 'POSITION'])
    df = df[df['BF(dB)'] > cutoff]
    df.reset_index(inplace=True, drop=True)
    df['i'] = df.index

    #pal = ['#66c2a5', '#fc8d62']
    #sns.set_palette([pal])
    plot = sns.relplot(data=df, x='i', y='BF(dB)', aspect=6, 
                        hue='CHROMOSOME', palette='colorblind', legend=None) 
    chrom_df=df.groupby('CHROMOSOME')['i'].median()
    plot.ax.set_xlabel('contig')
    plot.ax.set_xticks(chrom_df, labels=chrom_df.index)
    plot.ax.tick_params(axis='x', labelrotation=45)
    
    #add a horizontal line 
    #ax1, ax2 = plot.axes[0]

    #ax1.axhline(20, ls='--')
    #ax2.axhline(30, ls='--')
    plt.axhline(y=20, linestyle='--', color='grey', linewidth=1.5)
    plt.tight_layout()
    plt.savefig(out+'.png', bbox_inches='tight')
    '''
    with open(out+'_contig_names.txt', 'w') as f:
        for chrom in chrom_df.index:
    '''

def main():
    args = get_args()
    col_names = ['CHROMOSOME','POSITION','REF','ALT','COVARIABLE','MRK','M_Pearson',
                 'SD_Pearson','M_Spearman','SD_Spearman','BF(dB)','Beta_is','SD_Beta','_iseBPis']
    str_types = ['CHROMOSOME', 'REF', 'ALT', 'POSITION']
    int_types = ['POSITION']
    col_types = {}
    for name in col_names:
        if name in str_types:
            col_types[name] = str
        elif name in int_types:
            col_types[name] = int
        else:
            col_types[name] = np.float64
    #col_types = {name:str if name in str_types else name:np.float64 for name in col_names}
    df = pd.read_csv(args.tsv, sep='\t', names=col_names, dtype=col_types)
    print(df.head())
    #sys.exit()

    if args.contigs is not None:
        contigs = []
        with open(args.contigs, 'r') as f:
            for line in f.readlines():
                contigs.append(line.strip())
        #contigs = [contig for contig in args.contigs]
        #for contig in args.contigs:
        contigs = set(contigs)
        df = df[df['CHROMOSOME'].isin(contigs)]

    plot_manhattan(df, args.cutoff, args.o)
    #print(df.head())

    #print(df['chromosome'].unique())

if __name__ == '__main__':
    main()
