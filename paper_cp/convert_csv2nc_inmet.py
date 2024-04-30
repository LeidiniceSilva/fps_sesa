# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script convert .csv to .nc of weather station"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_inmet_stations import inmet

path='/afs/ictp.it/home/m/mda_silv/Documents/FPS_SESA/database/obs/inmet/inmet_br'

idx=10 # 1, 4, 5, 6 and 10
if idx == 1:
	nc_var = 'pre'
	unit_var = 'mm'
	name_var = 'Daily total of precipitation'
	std_var = 'precipitation'
elif idx == 4:
	nc_var = 'tmax'
	unit_var = 'C'
	name_var = 'Daily maximum air temperature'
	std_var = 'maximum air temperature'
elif idx == 5:
	nc_var = 'tmp'
	unit_var = 'C'
	name_var = 'Daily mean air temperature'
	std_var = 'mean air temperature'
elif idx == 6:
	nc_var = 'tmp'
	unit_var = 'C'
	name_var = 'Daily minimum air temperature'
	std_var = 'temperature'
else:
	nc_var = 'uv'
	unit_var = 'm s**-1'
	name_var = 'Daily mean wind speed'
	std_var = 'wind'

# create date list
dt = pd.date_range('2018-01-01','2021-12-31', freq='D').strftime('%Y-%m-%d').tolist()
	
for idx in range(1, 151):
																																																						
	print('Reading inmet station:', idx, inmet[idx][0])
	# Reading smn station
	data = pd.read_csv(os.path.join('{0}/inmet_csv/daily'.format(path), 'dados_{0}_D_2018-01-01_2021-12-31.csv'.format(inmet[idx][0])), skiprows=9, encoding='ISO-8859-1', decimal=',', delimiter=';')
	data_i = data.iloc[:,10]
	data_ii = data_i.replace(-999., np.nan)
	data_values = np.array(data_ii, dtype=float)

	data_dates  = []
	for i in range(len(dt)):
			
		data_dates.append('{0}'.format(dt[i]))
		print('Date organized:', data_dates[i], data_values[i])

	nc_output = '{0}/inmet_nc/daily/{1}/{1}_{2}_D_2018-01-01_2021-12-31.nc'.format(path, nc_var, inmet[idx][0])

	# create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions 	= 'CF-1.6'
	ds.title 		= 'Automatic weather stations.'
	ds.institution 	= 'Instituto Nacional de Meteorologia, INMET'
	ds.source 		= 'INMET Meteorological Database'
	ds.history 		= 'Rewrote via python script.'
	ds.references 	= 'https://bdmep.inmet.gov.br/.'
	ds.comment 		= 'This script convert .csv to .nc of weather station'
		
	ds.createDimension('time', None)

	time 				= ds.createVariable('time', float, ('time'))
	time.axis 			= 'L'
	time.calendar 		= 'standard'
	time.units			= 'days since {}'.format(data_dates[0])
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
