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
from dict_smn_i_stations import smn_i
from dict_inmet_stations import inmet

parser = argparse.ArgumentParser()
parser.add_argument('--var', required=True, help='Variable')
parser.add_argument('--mdl', required=True, help='Model')
parser.add_argument('--inst', required=True, help='Institution')
args = parser.parse_args()

var = args.var
mdl = args.mdl
inst = args.inst

if var == 'pr':
	nc_var = 'pr'
	unit_var = 'mm'
	name_var = 'Hourly total of precipitation'
	std_var = 'precipitation'
elif var == 'tas':
	nc_var = 'tas'
	unit_var = 'degrees C'
	name_var = 'Hourly mean of air temperature'
	std_var = 'temperature'
else: # 3
	nc_var = 'sfcWind' 
	unit_var = 'm.s**-1'
	name_var = 'Hourly mean of wind speed'
	std_var = 'wind speed'

if mdl == 'reg_ictp':
	exp = 'CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5_v1'
elif mdl == 'reg_ictp_pbl1':
	exp = 'CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0'
elif mdl == 'reg_ictp_pbl2':
	exp = 'CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0'
elif mdl == 'reg_usp':
	exp = 'CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_USP-RegCM471_v2'
elif mdl == 'wrf_ncar':
	exp = 'CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415'
else:
	exp = 'CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1'

# Paths
path_in = '/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/{0}'.format(mdl)
path_out = '/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/{0}/{1}'.format(mdl, nc_var)

# Skip list for INMET stations
skip_inmet = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] # 2018-2021


def extract_station(ts, station_code, station_name, institution):

	dt = pd.date_range('2018-06-01 00:00', '2021-05-31 23:00', freq='h')
	nc_output = os.path.join(path_out, f'{nc_var}_{station_code}_{station_name}_H_2018-06-01-2021-05-31.nc')
	ds_out = Dataset(nc_output, 'w', format='NETCDF4_CLASSIC')

	# Global attributes
	ds_out.Conventions  = 'CF-1.6'
	ds_out.title        = f'Model extracted for {institution} stations'
	ds_out.institution  = institution
	ds_out.source       = exp
	ds_out.history      = 'Extracted via Python script'
	ds_out.references   = 'https://cds.climate.copernicus.eu/'

	# Dimensions
	ds_out.createDimension('time', None)

	# Time variable
	time_var = ds_out.createVariable('time', float, ('time',))
	time_var.units = f'hours since {dt[0]}'
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
var_ds = Dataset('{0}/{1}_{2}_1hr_2018060100-2021053123.nc'.format(path_in, nc_var, exp))
lats = var_ds.variables['lat'][:]
lons = var_ds.variables['lon'][:]
time = var_ds.variables['time'][:]

# INMET stations
if inst == 'INMET':
	for station in range(360, 567):
		if station in skip_inmet:
			continue
		station_code, station_name, lat_s, lon_s, alt_s = inmet[station]
		print(f'Processing INMET station {station_code} - {station_name}')
		lat_idx = np.abs(lats - lat_s).argmin()
		lon_idx = np.abs(lons - lon_s).argmin()
		ts = var_ds.variables[nc_var][:, lat_idx, lon_idx]
		extract_station(ts, station_code, station_name, inst)
	var_ds.close()

# SMN stations
else:
	for station in smn_i:
		station_name, lat_s, lon_s = smn_i[station]
		station_code = f'SMN{station:03d}'
		print(f'Processing SMN station {station_code} - {station_name}')
		lat_idx = np.abs(lats - lat_s).argmin()
		lon_idx = np.abs(lons - lon_s).argmin()
		ts = var_ds.variables[nc_var][:, lat_idx, lon_idx]
		extract_station(ts, station_code, station_name, inst)
	var_ds.close()

exit()


