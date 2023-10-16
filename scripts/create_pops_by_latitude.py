import argparse


def get_args():
    parser = argparse.ArgumentParser(description='create populations by latitude', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--coords', type=str, required=True, help='file containing longitude and latitude corresponding to samples')
    parser.add_argument('--o', type=str, required=True, help='output filename')
    
    return parser.parse_args()
    
def parse_coords(coords_file):

    coords = {}
    with open(coords_file, 'r') as f:
        next(f)
        for line in f.readlines():
            sample, long, lat = line.strip().split()
            lat = round(float(lat))
            if lat in coords:
                coords[lat].append(sample)
            else:
                coords[lat] = [sample]

    return coords


def write_populations(coords, out):

    for lat in coords:
        with open(out+f'.{lat}.txt', 'w') as f:
            for sample in coords[lat]:
                f.write(f'{sample}\n')
            
def main():

    args = get_args()
    
    coords = parse_coords(args.coords)
    
    write_populations(coords, args.o)


if __name__ == '__main__':
    main()
