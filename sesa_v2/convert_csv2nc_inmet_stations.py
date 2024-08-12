# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Aug 12, 2024"
__description__ = "This script convert .csv to .nc from each INMET station"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_inmet_stations import inmet

idx=2
dt = 'H_2012-01-01_2021-12-31'
path = '/marconi/home/userexternal/mdasilva/user/mdasilva/FPS_SESA/database/obs/inmet/2012-2021'

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

# create date list
date = pd.date_range('2012-01-01','2022-01-01', freq='H')
date = date[:-1]
	
for j in range(1, 100):
	if j == 13:
		continue
	if j == 18:
		continue
	if j == 19:
		continue
	if j == 24:
		continue
	if j == 28:
		continue
	if j == 29:
		continue
	if j == 43:
		continue
	if j == 44:
		continue
	if j == 63:
		continue
	if j == 64:
		continue
	if j == 77:
		continue
	if j == 79:
		continue
	if j == 82:
		continue
	if j == 98:
		continue
																													
	# Reading inmet station																																																								
	print('Reading inmet weather station:', j, inmet[j][0])
	data = pd.read_csv(os.path.join('{0}/csv/'.format(path), 'dados_{0}_{1}.csv'.format(inmet[j][0], dt)), skiprows=9, encoding='ISO-8859-1', decimal=',', delimiter=';')
	data_i = data.iloc[:,2]
	data_ii = data_i.replace(-999., np.nan)
	data_values = np.array(data_ii, dtype=float)
	
	data_dates  = []
	for i in range(len(date)):
		data_dates.append('{0}'.format(date[i]))
		print('Date organized:', data_dates[i], data_values[i])
		
	nc_output = '{0}/nc/{1}_{2}_{3}.nc'.format(path, nc_var, inmet[j][0], dt)

	# create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions  	= 'CF-1.6'
	ds.title 		= 'INMET weather station.'
	ds.institution 	        = 'INMET.'
	ds.source 		= 'Automatic weather station.'
	ds.history 		= 'Rewrote via python script.'
	ds.references 	        = 'https://bdmep.inmet.gov.br/.'
	ds.comment 		= 'This script convert .csv to .nc from each inmet weather station'
		
	ds.createDimension('time', None)

	time 			= ds.createVariable('time', float, ('time'))
	time.axis 		= 'L'
	time.calendar 		= 'standard'
	time.units		= 'hours since {}:00'.format(data_dates[0])
	time[:]			= range(len(data_dates))
	var 			= ds.createVariable(nc_var,  float, ('time'))
	var.units 		= unit_var
	var.long_name 		= name_var
	var.standard_name 	= std_var
	var.missing_value 	= -999
	var[:] 			= data_values[:]
		
	ds.close()

	if os.path.exists(nc_output): 
		print('Done -->', nc_output)

exit()
