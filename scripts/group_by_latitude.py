import argparse


#this would need standardized across populations, but we only have one population for now 
def get_args():
    parser = argparse.ArgumentParser(description='create populations by latitude', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--coords', type=str, required=True, help='vcf.gz file')
    parser.add_argument('--o', type=str, required=True, help='output filename')
    
    return parser.parse_args()


def parse_coords(coords_file):

    coords = []
    with open(coords_file, 'r') as f:
        for line in f.readlines():
            sample, long, lat = line.split()
            coords.append((float(long), float(lat)))

    return coords

def main():
    
    args = get_args()
    coords = parse_coords(args.coords)


