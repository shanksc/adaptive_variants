import numpy as np
import xarray as xr
import warnings
import argparse

#supress warnings from xarray
warnings.simplefilter("ignore") 

def get_args():
    parser = argparse.ArgumentParser(description='format ENV file containing latitude min monthly temp and max monthly temp', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--coords', type=str, required=True, help='file containing longitude and latitude corresponding to samples')
    parser.add_argument('--samples', type=str, required=True, help='samples that we need from coords file', nargs='*') 
    parser.add_argument('--sst', type=str, required=True, help='nc file containing long term monthly temperatures')
    parser.add_argument('--o', type=str, required=True, help='output filename')
    
    return parser.parse_args()


def read_samples(samples_file):

    samples = []
    with open(samples_file, 'r') as f:
        for s in f.readlines():
            samples.append(s.strip())

    return set(samples)


def parse_coords(coords_file, samples):

    lats = []
    lons = []
    with open(coords_file, 'r') as f:
        for line in f.readlines():
            sample, lon, lat = line.split()
            if sample in samples: 
                lats.append(float(lat))
                lons.append(float(lon))

    return lats, lons


#tab delimited standardized across populations for BayEnv
def write_out(vars, out):
    
    with open(out+'.ENVS', 'w') as f:
        for i in range(len(vars)):
            standardize(vars[i])
            line = ''
            for n in vars[i][:-1]:
                line += f'{n}\t'
            f.write(line+f'{vars[i][-1]}\n')


def get_min_max_temps(lats, lons, sst):

    sst = xr.open_dataset(sst)

    min_temps = []
    max_temps = []
    for i in range(len(lats)):
        #lon needs to be 0-360 deg format, not [-180, 180]
        sliced = sst['sst'].sel(lat=lats[i], lon=lons[i]%360, method='nearest')
        #handle edge case where the "nearest" is actually land
        #this is hacky but xarray has virtually no good solution and we don't want to interpolate
        #can't even drop na values effectively
        #so we're moving nearest by longitude, and taking the smallest step
        #need a better solution to scale
        arr = sliced.to_numpy()
        if np.isnan(arr).any():
            sliced = sst['sst'].sel(lat=lats[i], lon=(lons[i]%360)-.5, method='nearest')
            arr = sliced.to_numpy()
            if np.isnan(arr).any():
                raise ValueError("Population is not close to valid lat/long coordinates")
        
        min_temps.append(np.amax(arr))
        max_temps.append(np.amin(arr))

    return min_temps, max_temps

#standardize list in place
def standardize(a):
    
    u = np.mean(a)
    sd = np.std(a)
    
    for i, n in enumerate(a):
        a[i] = (n - u) / sd


def main():
    args = get_args()
    samples = [read_samples(s) for s in args.samples] 

    lats = []
    lons = []
    for s in samples:
        lat, lon = parse_coords(args.coords, s)
        lats.append(np.mean(lat))
        lons.append(np.mean(lon))

    min_temps, max_temps = get_min_max_temps(lats, lons, args.sst)

    write_out([lats, min_temps, max_temps], args.o)


if __name__ == "__main__":
    main()
