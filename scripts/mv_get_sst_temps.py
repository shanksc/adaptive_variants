import numpy as np
import pandas as pd
import xarray as xr
import fsspec


MIN_LAT=25
MAX_LAT=50
MIN_LON=-130
MAX_LON=-110

sst = xr.open_zarr('https://mur-sst.s3.us-west-2.amazonaws.com/zarr-v1', consolidated=True)


#we don't need to filter on ice given out latitude
#print(ds_sst)
#we probably won't SST monthly min, monthly max? 
sliced = sst['analysed_sst'].sel(lat=slice(MIN_LAT,MAX_LAT),lon=slice(MIN_LON,MAX_LON))
print(sliced)
print('SLICED')
sst_monthly = sliced.resample(time='1MS').mean('time',keep_attrs=True,skipna=False)
print('MONTHLY')
mean_monthly = sst_monthly.groupby('time.month').mean('time',keep_attrs=True,skipna=False)
print(mean_monthly)
print('MEAN MONTHLY')
#average over monthly variables
mean_monthly.to_netcdf("mean_monthly.nc")


