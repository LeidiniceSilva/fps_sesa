# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "10/17/2022"
__description__ = "This script convert .csv to .nc from each inmet station"

import os
import numpy as np

from netCDF4 import Dataset
from dict_inmet_stations_code import codes
from dict_inmet_stations_name import names

idx=6

if idx == 2:
	nc_var = 'pre'
	unit_var = 'mm'
	name_var = 'Hourly total of precipitation'
	std_var = 'precipitation'
	
elif idx == 6:
	nc_var = 'tmp'
	unit_var = 'C'
	name_var = 'Hourly mean of temperature'
	std_var = 'temperature'

else:
	nc_var = 'uv'
	unit_var = 'm s**-1'
	name_var = 'Hourly mean of wind'
	std_var = 'wind'

dt = 'H_2018-01-01_2021-12-31'

for j in range(1, 289):

	if j == 4:
		continue
	if j == 12:
		continue
	if j == 45:
		continue
	if j == 55:
		continue
	if j == 77:
		continue
	if j == 98:
		continue
	if j == 99:
		continue
	if j == 118:
		continue
	if j == 122:
		continue
	if j == 130:
		continue
	if j == 135:
		continue
	if j == 151:
		continue
	if j == 155:
		continue
	if j == 159:
		continue
	if j == 160:
		continue
	if j == 163:
		continue
	if j == 164:
		continue
	if j == 181:
		continue
	if j == 183:
		continue
	if j == 186:
		continue
	if j == 187:
		continue
	if j == 188:
		continue
	if j == 209:
		continue
	if j == 216:
		continue
	if j == 228:
		continue
	if j == 236:
		continue
	if j == 268:
		continue
	if j == 287:
		continue
																																																											
	# This load the file (date, hour, variable)
	data = np.loadtxt(os.path.join('/home/nice/Downloads/FPS_SESA/inmet/inmet_used/', 'dados_{0}_{1}.csv'.format(codes[j], dt)), dtype='str', delimiter=',', unpack=True)
	data = data[:,1:]
	data_values = np.where(data[idx,:] == 'null', -999., data[idx,:])
	data_values = np.array(data_values, dtype=float)
	
	print()
	print('Reading inmet weather station:', j, codes[j], names[j][1])

	data_dates  = []
	for i in range(len(data[0,:])):
		
		data_dates.append('{0} {1:02d}:{2:02d}'.format(data[0,i], int(int(data[1,i])/100), int(data[1,i][-2:])))

		if i < 24:
			print('Date organized:', data_dates[i], data_values[i])

	nc_output = '/home/nice/Downloads/FPS_SESA/inmet/inmet_nc/{0}_{1}_{2}.nc'.format(nc_var, codes[j], dt)

	# create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions 	= 'CF-1.6'
	ds.title 		= 'INMET weather station.'
	ds.institution 	= 'INMET.'
	ds.source 		= 'Automatic weather station.'
	ds.history 		= 'Rewrote via python script.'
	ds.references 	= 'https://bdmep.inmet.gov.br/.'
	ds.comment 		= 'This script convert .csv to .nc from each inmet weather station'
		
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















