import argparse
import random


def get_args():
    parser = argparse.ArgumentParser(description='randomly sample baypass SNPS', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--SNPS', type=str, required=True, help='SNPS file')
    #this can be expanded to accept multiple samples files which can then support different populations
    parser.add_argument('--k', type=int,required=True, help='number of SNPS to randomly sample')
    parser.add_argument('--o', type=str, required=True, help='output filename')

    return parser.parse_args()


def read_SNPS(SNPS):
    with open(SNPS, 'r') as f:
        return list(f.readlines())


def write_random_sample(lines, k, out):
    
    #this really shouldn't happen but if it does...
    if k <= len(lines):
        #enumerate to keep track of indicies
        sample = random.sample(list(enumerate(lines)), k)  
        # 1000-200 800 1000      
        #idx = random.randint(0, len(line_pairs) - k)
    else:
        raise ValueError(f'SNPS file not large enough for sample {k} large')

    with open(out+'.SNPS', 'w') as f:
        for _, line in sample:
            f.write(line)
    
    #we want to write out indicies so we can map back to clean.INFO files
    with open(out+'.indices', 'w') as f:
        for idx, _ in sample:
            f.write(f'{idx}\n')

def main():

    args = get_args()
    line_pairs = read_SNPS(args.SNPS)
    write_random_sample(line_pairs, args.k, args.o)


if __name__ == '__main__':
    main()

