# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script convert .csv to .nc of the INMET weather station"

import os
import numpy as np

from netCDF4 import Dataset
from dict_sesa_inmet_stations import inmet

idx=4

if idx == 2:
	nc_var = 'pre'
	unit_var = 'mm'
	name_var = 'Hourly total of precipitation'
	std_var = 'precipitation'
elif idx == 4:
	nc_var = 'rad'
	unit_var = 'kj m**-2'
	name_var = 'Hourly mean of solar radiation'
	std_var = 'solar radiation'
elif idx == 7:
	nc_var = 'tmp'
	unit_var = 'C'
	name_var = 'Hourly mean of temperature'
	std_var = 'temperature'
elif idx == 8:
	nc_var = 'dtmp'
	unit_var = 'C'
	name_var = 'Hourly mean of dew point temperature'
	std_var = 'dew point temperature'
else:
	nc_var = 'uv'
	unit_var = 'm s**-1'
	name_var = 'Hourly mean of wind'
	std_var = 'wind'
	
dt = 'H_2018-01-01_2021-12-31'

for j in range(1, 155):
																																																						
	# This load the file (date, hour, variable)
	data = np.loadtxt(os.path.join('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_used_sesa', 'dados_{0}_{1}.csv'.format(inmet[j][0], dt)), dtype='str', delimiter=',', unpack=True)
	data = data[:,1:]
	data_values = np.where(data[idx,:] == 'null', -999., data[idx,:])
	data_values = np.array(data_values, dtype=float)

	print()
	print('Reading inmet weather station:', j, inmet[j][0], inmet[j][1])
	
	data_dates  = []
	for i in range(len(data[0,:])):
		
		data_dates.append('{0} {1:02d}:{2:02d}'.format(data[0,i], int(int(data[1,i])/100), int(data[1,i][-2:])))

		if i < 24:
			print('Date organized:', data_dates[i], data_values[i])

	nc_output = '/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc_sesa/{0}/{1}_{2}_{3}.nc'.format(nc_var, nc_var, inmet[j][0], dt)

	# create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions 	= 'CF-1.6'
	ds.title 		= 'INMET automatic weather station.'
	ds.institution 	= 'INMET.'
	ds.source 		= 'INMET automatic weather station.'
	ds.history 		= 'Rewrote via python script.'
	ds.references 	= 'https://bdmep.inmet.gov.br/.'
	ds.comment 		= 'This script convert .csv to .nc from inmet automatic weather station'
		
	ds.createDimension('time', None)

	time 				= ds.createVariable('time', float, ('time'))
	time.axis 			= 'L'
	time.calendar 		= 'standard'
	time.units			= 'hours since {}:00'.format(data_dates[0])
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















