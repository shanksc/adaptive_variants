import argparse
import pandas as pd
import sys

#add a filter for contig size?
def get_args():
    parser = argparse.ArgumentParser(description='plot manhattan plot', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--bfs', type=str, required=True, help='bfs file')
    parser.add_argument('--info', type=str, required=True, help='INFO file')
    #this can be expanded to accept multiple samples files which can then support different populations
    parser.add_argument('--contig_sizes', type=str, required=False, help='tsv containing contig size')
    parser.add_argument('--min_contig_size', type=int, required=False, help='minimum size for contigs', default=100000)
    parser.add_argument('--indices', type=str, help='indicies of random sample', default=None)
    parser.add_argument('--env_labels', type=str, help='indicies', nargs='*')
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
    #print(len(qualified_contigs))
    return qualified_contigs
            

def get_info(info_file, indices):

    if indices is not None:
        indices = set(indices)
    #let info be by idx, mapping to chrom and position, so that it corresponds to SNPS file
    info = {}

    with open(info_file, 'r') as f:
        for idx, line in enumerate(f.readlines()):
            chrom, pos, ref, alt = line.strip().split('\t')
            if indices is not None:
                if idx in indices:
                    info[idx] = [chrom, int(pos)]
            else:
                info[idx] = [chrom, int(pos)]

    return info


def get_indices(indicies_file):

    indicies = []

    with open(indicies_file, 'r') as f:
        for line in f.readlines():
            idx = int(line.strip())
            indicies.append(idx)
    
    return indicies


def get_bfs(bfs_file, env_labels):

    label_to_bfs = {label:[] for label in env_labels}

    with open(bfs_file, 'r') as f:
        for line in f.readlines():
            #if indices is not None:
                    #get all bfs, excluding chunk file name/location in first col
            bfs = line.strip().split('\t')[1:]
            #print(bfs)
            #sys.exit()

            for i, label in enumerate(env_labels):
                label_to_bfs[label].append(bfs[i])
            #else:
                #bfs[idx] = [line.strip().split()[1:]]
    #print(len(label_to_bfs['latitude']))
    #sys.exit()
    return label_to_bfs
                    

def main():
    args = get_args()
    #print(args.env_labels)
    indices = get_indices(args.indices)
    label_to_bfs = get_bfs(args.bfs, args.env_labels)
    info = get_info(args.info, indices)

    #print(len(bfs))
    #print(len(info))
    #for idx in bfs:
    #    info[idx].extend(bfs[idx])
    
    #write out tsv
    chroms = []
    positions = []
    
    for idx in info:
        chroms.append(info[idx][0])
        positions.append(info[idx][1])

        #for i in range(2,len(info[idx])):
        #    bfs[i].append(info[idx][i])
    
    #build dict for df
    #print(len(chroms))
    #print(len(positions))
    df_dict = {'chromosome': chroms, 'position': positions}
    
    for label, bfs in label_to_bfs.items():
        df_dict[label] = bfs
        #print(len(df_dict[label]))

    df = pd.DataFrame.from_dict(df_dict)
    
    #now filter on qualified contigs
    if args.contig_sizes is not None:
        qualified_contigs = get_qualified_contigs(args.contig_sizes, args.min_contig_size)
        df = df[df.chromosome.isin(qualified_contigs)]

    df.to_csv(path_or_buf=args.o, sep='\t', index=False)
    

if __name__ == '__main__':
    main()