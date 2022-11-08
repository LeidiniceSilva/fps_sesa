# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "10/17/2022"
__description__ = "This script convert .csv to .nc from each Uruguai station"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_stations_urug_smn import urug_smn

# Reading Uruguai weather station
df = pd.read_csv(os.path.join('/home/nice/Documentos/FPS_SESA/urug_smn/', 'PP_URUG_SMN_2018_2021.csv'))
df['Date'] = pd.to_datetime(df['Date'])
df = df.set_index('Date').resample('H').sum()
df.iloc[:,:].to_csv('/home/nice/Documentos/FPS_SESA/urug_smn/PP_URUG_SMN_H_2018-01-01_2021-12-31.csv')

nc_var = 'pre'
unit_var = 'mm'
name_var = 'Hourly total of precipitation'
std_var = 'precipitation'
		 
for idx in range(0, 73):

	print()
	print('Reading Uruguai weather station:', idx+1, urug_smn[idx+1][0])
	# Reading Uruguai weather station
	data = np.loadtxt(os.path.join('/home/nice/Documentos/FPS_SESA/urug_smn', 'PP_URUG_SMN_H_2018-01-01_2021-12-31.csv'), dtype='str', delimiter=',', unpack=True)
	data = data[:,1:-1]
	data_values = np.where(data[idx+1,:] == 'null', -999., data[idx+1,:])
	data_values = np.array(data_values, dtype=float)

	data_dates  = []
	for i in range(len(data[0,:])):
			
		data_dates.append('{0}'.format(data[0,i]))
		
		if i < 24:
			print('Date organized:', data_dates[i], data_values[i])
	
	nc_output = '/home/nice/Documentos/FPS_SESA/urug_smn/urug_smn_nc/pre_{0}_H_2018-01-01_2021-12-31.nc'.format(urug_smn[idx+1][0])

	print('Create netcdf')
	# Create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions 	= 'CF-1.6'
	ds.title 		= 'Uruguai weather station.'
	ds.institution 	= 'SMN.'
	ds.source 		= 'Automatic weather station.'
	ds.history 		= 'Rewrote via python script.'
	ds.references 	= 'https://www.smn.gob.ar/.'
	ds.comment 		= 'This script convert .csv to .nc from each Uruguai weather station'
		
	ds.createDimension('time', None)

	time 				= ds.createVariable('time', float, ('time'))
	time.axis 			= 'L'
	time.calendar 		= 'standard'
	time.units			= 'hours since {}'.format(data_dates[0])
	time[:]				= range(len(data_dates))

	var 				= ds.createVariable(nc_var,  float, ('time'))
	var.units 			= unit_var
	var.long_name 		= name_var
	var.standard_name 	= std_var
	var.missing_value 	= -999
	var[:] 				= data_values[:]
		
	ds.close()

	if os.path.exists(nc_output): 
		print('Done -->', nc_output)

exit()
