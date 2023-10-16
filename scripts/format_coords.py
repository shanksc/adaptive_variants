#DEPRACATED 

import argparse
import numpy as np


def get_args():
    parser = argparse.ArgumentParser(description='format coords file to ENVIRONFILE containing latitude', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--coords', type=str, required=True, help='file containing longitude and latitude corresponding to samples')
    parser.add_argument('--samples', type=str, required=True, help='samples that we need from coords file', nargs='*') 
    parser.add_argument('--o', type=str, required=True, help='output filename')
    
    return parser.parse_args()


def read_samples(samples_file):

    samples = []
    with open(samples_file, 'r') as f:
        for s in f.readlines():
            samples.append(s.strip())

    return set(samples)
            

def parse_coords(coords_file, samples):

    coords = []
    with open(coords_file, 'r') as f:
        for line in f.readlines():
            sample, long, lat = line.split()
            if sample in samples: 
                coords.append(float(lat))

    return np.array(coords)


#tab delimited standardized across populations for BayEnv
def write_coords(coords, out):
    
    with open(out+'.ENVS', 'w') as f:
        line = ''
        for i in range(len(coords)-1):
            line += f'{coords[i]}\t'
        f.write(line+f'{coords[-1]}\n')


def main():

    args = get_args()
    samples = [read_samples(s) for s in args.samples] 
    coords = []
    for s in samples:
        coords.append(np.mean(parse_coords(args.coords, s)))
    
    #now we standardize
    sd = np.std(coords)
    u = np.mean(coords)
    for i, lat in enumerate(coords):
        coords[i] = (lat - u) / sd


    write_coords(coords, args.o) 


if __name__ == '__main__':
    main()




