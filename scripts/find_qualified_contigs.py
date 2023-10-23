import argparse
import pandas as pd
import sys

#add a filter for contig size?
def get_args():
    parser = argparse.ArgumentParser(description='create txt containing names of qualified contigs', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    #this can be expanded to accept multiple samples files which can then support different populations
    parser.add_argument('--contig_sizes', type=str, required=True, help='tsv containing contig size')
    parser.add_argument('--min_contig_size', type=int, required=False, help='minimum size for contigs', default=100000)
    parser.add_argument('--o', type=str, required=True, help='output filename')

    return parser.parse_args()


def get_qualified_contigs(contig_sizes_file, min_contig_size):
    
    qualified_contigs = set([])
    
    with open(contig_sizes_file) as f:
        for line in f.readlines():
            #should always be second column
            chrom, size, _, _, _ = line.strip().split('\t')
            size = int(size)
            if size >= min_contig_size:
                qualified_contigs.add(chrom)
    
    return qualified_contigs


def main():
    args = get_args()
    contigs = get_qualified_contigs(args.contig_sizes, args.min_contig_size)
    
    with open(args.o+'.txt', 'w') as f:
        for contig in contigs:
            f.write(f'{contig}\n')
        
            
if __name__ == '__main__':
    main()