# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script convert .csv to .nc of weather station"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_smn_ii_stations import smn_ii

var = 'pre'

if var == 'pre':
	nc_var = 'pre'
	unit_var = 'pre'
	name_var = 'Daily total of precipitation'
	std_var = 'precipitation'
elif var == 'tmax':
	nc_var = 'tmax'
	unit_var = 'C'
	name_var = 'Daily mean of maximum temperature'
	std_var = 'maximum temperature'
elif var == 'tmin':
	nc_var = 'tmin'
	unit_var = 'C'
	name_var = 'Daily mean of minimum temperature'
	std_var = 'minimum temperature'
else:
	nc_var = 'uv'
	unit_var = 'C'
	name_var = 'Daily mean of wind speed'
	std_var = 'wind speed'

# create date list
dt = pd.date_range('1979-01-01','2021-12-31', freq='D').strftime('%Y-%m-%d').tolist()

list_st = np.arange(0, 109)

for idx in list_st:

	print('Reading smn station:', idx+1, smn_ii[idx+1][0])
	# Reading smn station
	data = np.loadtxt(os.path.join('/home/mda_silv/users/FPS_SESA/database/obs/smn_ii/', 'precip_BrasUruArgPar_D_19790101-20211231.csv'), dtype='str', delimiter=',', unpack=True)
	file_path = os.path.join('/home/mda_silv/users/FPS_SESA/database/obs/smn_ii/', 'precip_BrasUruArgPar_D_19790101-20211231.csv')
	df = pd.read_csv(file_path, sep=',', na_values=['null'], parse_dates=[0])
	data_values = (df.iloc[:, idx + 1].fillna(-999.).to_numpy(dtype=float))
	
	data_dates  = []
	for i in range(len(dt)):
			
		data_dates.append('{0}'.format(dt[i]))
		print('Date organized:', data_dates[i], data_values[i])
		
	nc_output = '/home/mda_silv/users/FPS_SESA/database/obs/smn_ii/smn_nc/pre/pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[idx+1][0])

	print('Create netcdf')
	# Create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions 	       = 'CF-1.6'
	ds.title 	       = 'Weather stations.'
	ds.institution 	       = 'SMN.'
	ds.source 	       = 'Automatic weather station.'
	ds.history 	       = 'Rewrote via python script.'
	ds.references 	       = 'https://www.smn.gob.ar/.'
	ds.comment 	       = 'This script convert .csv to .nc of weather station'
	
	ds.createDimension('time', None)

	time 		       = ds.createVariable('time', float, ('time'))
	time.axis 	       = 'L'
	time.calendar 	       = 'standard'
	time.units	       = 'days since {}'.format(data_dates[0])
	time[:]		       = range(len(data_dates))

	var 		       = ds.createVariable(nc_var,  float, ('time'))
	var.units 	       = unit_var
	var.long_name 	       = name_var
	var.standard_name      = std_var
	var.missing_value      = -999
	var[:] 		       = data_values[:]

	ds.close()

	if os.path.exists(nc_output): 
		print('Done -->', nc_output)

exit()
