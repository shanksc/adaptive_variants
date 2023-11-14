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
    parser.add_argument('--annotate', type=str, required=False, help='tsv containing ANP receptors to annotate', default=None)
    parser.add_argument('--o', type=str, required=True, help='output filename')

    return parser.parse_args()


def plot_manhattan(df, cutoff, var, out, annotate):
    #calculate decibans of BFs
    #df['dB'] = 10 * np.log10(df['BF(dB)'])
    df = df.sort_values(['#[1]CHROM', '[2]POS'])
    df = df[df[var] > cutoff]
    df.reset_index(inplace=True, drop=True)
    df['i'] = df.index

    #pal = ['#66c2a5', '#fc8d62']
    #sns.set_palette([pal])
    #change to binary colors
    #plot = sns.relplot(data=df, x='i', y=var, aspect=6, 
    #                    hue='#[1]CHROM', palette='colorblind', legend=None)
    #for alternating between two colors
    '''
    palette = []
    for i, chrom in enumerate(df['#[1]CHROM'].unique()):
        if i % 2 == 0:
            palette.append('b')
        else:
            palette.append('grey')
    '''    
    
    plot = sns.relplot(data=df, x='i', y=var, aspect=3, 
                        hue='#[1]CHROM', palette=sns.color_palette('colorblind'), legend=None, linewidth=0)
    plt.margins(x=0.01, y=0.01)
    #chroms = df['#[1]CHROM'].unique()
    #plt.xticks(ticks=[x for x in range(1, len(chroms))], labels=[x for x in range(1, len(chroms))])
    chrom_df=df.groupby('#[1]CHROM')['i'].median()
    plot.ax.set_xlabel('contig')
    plot.ax.set_xticks(chrom_df, labels=[x for x in range(1, len(chrom_df.index)+1)])
    plot.ax.tick_params(axis='x', labelsize=7)
    
    #plot.ax.set_xticks(chrom_df, labels=)

    if annotate is not None:
        anno = pd.read_csv(annotate, sep='\t')
        anno = anno.sort_values(['#[1]CHROM', '[2]POS'])
        print(anno)
        not_labeled=True
        for idx, row in anno.iterrows():
            #MRK variable
            #i = df.index.get_loc(row[''])
            print(row['#[1]CHROM'],row['[2]POS'])
            #i = df.index.get_loc(row['[2]POS'])
            i = df.loc[df['[2]POS'] == row['[2]POS']]

            print(i)
            #sys.exit()
            #if row.index['MRK'] == df.iloc[i]['MRK']:
            #x = df.iloc[i]['']
            #print(row['MRK'])
            print(f'axes shape: {plot.axes.shape}')
            ax = plot.axes[0,0]
            #if i[0] == i[1]:
            #    i = i[0]
            #ax.text(i['i']+100, row['BF(dB)']+100, 'ANP Receptor')
            #change this in future, should just gather all of the points then plot wiht one func call...
            if not_labeled:
                ax.plot(i['i'], row['BF(dB)'], color='black', marker='X', markersize=40, label='ANP receptor', linestyle='None')
                not_labeled=False
            else:
                ax.plot(i['i'], row['BF(dB)'], color='black', marker='X', markersize=40)
            #zoom inset 
            #max_x = anno['[2]POS'].max() + 1e5
            #min_x = anno['[2]POS'].min() - 1e5
            #axins = ax.inset_axes([0.5, 0.5, 0.47, 0.47])

            #ax.indicate_inset_zoom(axins, edgecolor="black")
            #get df with contig 14 
            #contig_14 = df[df['#[1]CHROM'] == 'JALGQA010000001.1']

            #axins.scatter(contig_14.index, contig_14['BF(dB)'])
            #sns.distplot(d, ax=ax2)
            #plot = sns.relplot(data=df, ax=axins, x='i', y=var, aspect=3, 
            #            hue='#[1]CHROM', palette=sns.color_palette('colorblind'), legend=None, linewidth=0)
            #axins.scatter()
            #axins.set_title('contig 14')
            #all positions containing contig 14 
            #max_x = anno['[2]POS'].max() + 1e5
            #min_x = anno['[2]POS'].min() - 1e5
            #axins.set_xlim([min_x, max_x])
            
        #df = df[df[var] > cutoff]
        #df.reset_index(inplace=True, drop=True)
        #df['i'] = df.index
        #plot = sns.relplot(data=df, x='i', y=var, aspect=3, 
        #                 legend=None, linewidth=0)
        print(df)

    #add a horizontal line 
    #ax1, ax2 = plot.axes[0]

    #ax1.axhline(20, ls='--')
    #ax2.axhline(30, ls='--')
    if var == 'BF(dB)':
        plt.axhline(y=20, linestyle='--', color='black', linewidth=1.5)
    #plt.tight_layout()
    #plt.legend()
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
    plot_manhattan(df, args.cutoff, var, args.o, args.annotate)
    #print(df.head())

    #print(df['#[1]CHROM'].unique())

if __name__ == '__main__':
    main()
