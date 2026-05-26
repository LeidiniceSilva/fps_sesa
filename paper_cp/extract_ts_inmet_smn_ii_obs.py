# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Feb 23, 2026"
__description__ = "This script extract timeseries"

import os
import argparse
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_smn_ii_stations import smn_ii

parser = argparse.ArgumentParser()
parser.add_argument('--var', required=True, help='Variable')
parser.add_argument('--inst', required=True, help='Institution')
args = parser.parse_args()

var = args.var
inst = args.inst

if var == 'tp':
	nc_var = 'tp'
	unit_var = 'mm'
	name_var = 'Daily total of precipitation'
	std_var = 'precipitation'
elif var == 't2m':
	nc_var = 't2m'
	unit_var = 'degrees C'
	name_var = 'Daily mean of air temperature'
	std_var = 'temperature'
else: # 3
	nc_var = 'ws10' 
	unit_var = 'm.s**-1'
	name_var = 'Daily mean of wind speed'
	std_var = 'wind speed'

path_in = '/home/mda_silv/users/FPS_SESA/database/obs/era5'
path_out = '/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/{0}'.format(nc_var)


def extract_station(ts, station_code, station_name, institution):

	dt = pd.date_range('2018-06-01 00:00', '2021-05-31 23:00', freq='d')
	nc_output = os.path.join(path_out, f'{nc_var}_{station_code}_{station_name}_D_2018-06-01-2021-05-31.nc')
	ds_out = Dataset(nc_output, 'w', format='NETCDF4_CLASSIC')

	# Global attributes
	ds_out.Conventions  = 'CF-1.6'
	ds_out.title        = f'ERA5 extracted for {institution} stations'
	ds_out.institution  = institution
	ds_out.source       = 'ERA5 CSAM-4i reanalysis'
	ds_out.history      = 'Extracted via Python script'
	ds_out.references   = 'https://cds.climate.copernicus.eu/'

	# Dimensions
	ds_out.createDimension('time', None)

	# Time variable
	time_var = ds_out.createVariable('time', float, ('time',))
	time_var.units = f'days since {dt[0]}'
	time_var.calendar = 'standard'
	time_var[:] = np.arange(len(dt))

	# Data variable
	var_nc = ds_out.createVariable(nc_var, float, ('time',))
	var_nc.units = unit_var
	var_nc.long_name = name_var
	var_nc.standard_name = std_var
	var_nc.missing_value = np.nan
	var_nc[:] = ts

	ds_out.close()
	print(f'Done --> {nc_output}')
	
	
# Load ERA5 netCDF
var_ds = Dataset('{0}/{1}_CSAM-4i_ERA5_day_2018060100-2021053123.nc'.format(path_in, nc_var))
lats = var_ds.variables['lat'][:]
lons = var_ds.variables['lon'][:]
time = var_ds.variables['time'][:]

# SMN stations
for station in smn_ii:

	station_name, lat_s, lon_s = smn_ii[station]
	station_code = f'SMN{station:03d}'
	print(f'Processing SMN station {station_code} - {station_name}')
	print(f'Lat: {lat_s} - Lon: {lon_s}')
	
	lat_idx = np.abs(lats - lat_s).argmin()
	lon_idx = np.abs(lons - lon_s).argmin()
	ts = var_ds.variables[nc_var][:, lat_idx, lon_idx]
	extract_station(ts, station_code, station_name, inst)
	
var_ds.close()
exit()


