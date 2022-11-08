# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "10/17/2022"
__description__ = "This script convert .csv to .nc from each Argentina station"

import os
import numpy as np
import pandas as pd

from netCDF4 import Dataset
from dict_stations_arg_emas import arg_emas

# ~ var = 'precip'
# ~ nc_var = 'precip'
# ~ unit_var = 'mm'
# ~ name_var = 'Hourly total of precipitation'
# ~ std_var = 'precipitation'

var = 'temp'
nc_var = 'temp'
unit_var = 'C'
name_var = 'Hourly mean of temperature'
std_var = 'temperature'

dt_i = pd.date_range('2018-01-01','2021-12-31', freq='MS').strftime('%Y-%m-%d').tolist()
dt_f = pd.date_range('2018-01-31','2021-12-31', freq='M').strftime('%Y-%m-%d').tolist()
dt_h = pd.date_range('2018-01-01 00:00:00','2021-12-31 23:00:00', freq='H').strftime('%Y-%m-%d %H:%M:%S').tolist()

dict_dt = {1: ['2018-01-01 00:00:00', '2018-01-31 23:00:00'],
2: ['2018-02-01 00:00:00', '2018-02-28 23:00:00'],
3: ['2018-03-01 00:00:00', '2018-03-31 23:00:00'],
4: ['2018-04-01 00:00:00', '2018-04-30 23:00:00'],
5: ['2018-05-01 00:00:00', '2018-05-31 23:00:00'],
6: ['2018-06-01 00:00:00', '2018-06-30 23:00:00'],
7: ['2018-07-01 00:00:00', '2018-07-31 23:00:00'],
8: ['2018-08-01 00:00:00', '2018-08-31 23:00:00'],
9: ['2018-09-01 00:00:00', '2018-09-30 23:00:00'],
10: ['2018-10-01 00:00:00', '2018-10-31 23:00:00'],
11: ['2018-11-01 00:00:00', '2018-11-30 23:00:00'],
12: ['2018-12-01 00:00:00', '2018-12-31 23:00:00'],
13: ['2019-01-01 00:00:00', '2019-01-31 23:00:00'],
14: ['2019-02-01 00:00:00', '2019-02-28 23:00:00'],
15: ['2019-03-01 00:00:00', '2019-03-31 23:00:00'],
16: ['2019-04-01 00:00:00', '2019-04-30 23:00:00'],
17: ['2019-05-01 00:00:00', '2019-05-31 23:00:00'],
18: ['2019-06-01 00:00:00', '2019-06-30 23:00:00'],
19: ['2019-07-01 00:00:00', '2019-07-31 23:00:00'],
20: ['2019-08-01 00:00:00', '2019-08-31 23:00:00'],
21: ['2019-09-01 00:00:00', '2019-09-30 23:00:00'],
22: ['2019-10-01 00:00:00', '2019-10-31 23:00:00'],
23: ['2019-11-01 00:00:00', '2019-11-30 23:00:00'],
24: ['2019-12-01 00:00:00', '2019-12-31 23:00:00'],
25: ['2020-01-01 00:00:00', '2020-01-31 23:00:00'],
26: ['2020-02-01 00:00:00', '2020-02-29 23:00:00'],
27: ['2020-03-01 00:00:00', '2020-03-31 23:00:00'],
28: ['2020-04-01 00:00:00', '2020-04-30 23:00:00'],
29: ['2020-05-01 00:00:00', '2020-05-31 23:00:00'],
30: ['2020-06-01 00:00:00', '2020-06-30 23:00:00'],
31: ['2020-07-01 00:00:00', '2020-07-31 23:00:00'],
32: ['2020-08-01 00:00:00', '2020-08-31 23:00:00'],
33: ['2020-09-01 00:00:00', '2020-09-30 23:00:00'],
34: ['2020-10-01 00:00:00', '2020-10-31 23:00:00'],
35: ['2020-11-01 00:00:00', '2020-11-30 23:00:00'],
36: ['2020-12-01 00:00:00', '2020-12-31 23:00:00'],
37: ['2021-01-01 00:00:00', '2021-01-31 23:00:00'],
38: ['2021-02-01 00:00:00', '2021-02-28 23:00:00'],
39: ['2021-03-01 00:00:00', '2021-03-31 23:00:00'],
40: ['2021-04-01 00:00:00', '2021-04-30 23:00:00'],
41: ['2021-05-01 00:00:00', '2021-05-31 23:00:00'],
42: ['2021-06-01 00:00:00', '2021-06-30 23:00:00'],
43: ['2021-07-01 00:00:00', '2021-07-31 23:00:00'],
44: ['2021-08-01 00:00:00', '2021-08-31 23:00:00'],
45: ['2021-09-01 00:00:00', '2021-09-30 23:00:00'],
46: ['2021-10-01 00:00:00', '2021-10-31 23:00:00'],
47: ['2021-11-01 00:00:00', '2021-11-30 23:00:00'],
48: ['2021-12-01 00:00:00', '2021-12-31 23:00:00']}

for i in range(1, 88):
	
	station_value = []
	idx = arg_emas[i][0]
	print('Read weather station:', i, idx)

	for dt in range(0, 48):
		
		# ~ df = pd.read_csv(os.path.join('/home/nice/Documentos/FPS_SESA/arg_emas/BCER/', '{0}_{1}_EMAs.csv'.format(dt_i[dt], dt_f[dt])))
		# ~ df = df.loc[df['idEstacion']==int(idx)]
		# ~ df['fechaHora'] = pd.to_datetime(df['fechaHora'])
		# ~ df = df.set_index('fechaHora').resample('H').mean()	# sum.() to precip and .mean() to temp
		# ~ df.iloc[:,4].to_csv('/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas/{0}_{1}_H_{2}_{3}_emas.csv'.format(var, idx, dt_i[dt], dt_f[dt]))
		# ~ print('Write weather station: {0}_{1}_H_{2}_{3}_emas.csv'.format(var, idx, dt_i[dt], dt_f[dt]))

		# ~ df = pd.read_csv(os.path.join('/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas/', '{0}_{1}_H_{2}_{3}_emas.csv'.format(var, idx, dt_i[dt], dt_f[dt])), index_col='fechaHora')
		# ~ df.head()
		# ~ df.index = pd.DatetimeIndex(df.index)
		# ~ df = df.reindex(pd.date_range(dict_dt[dt+1][0], dict_dt[dt+1][1], freq='H'), fill_value=-999.)
		# ~ df.iloc[:,:].to_csv('/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas/{0}_{1}_{2}_{3}_emas.csv'.format(var, idx, dt_i[dt], dt_f[dt]))
		# ~ print('Write weather station: {0}_{1}_{2}_{3}_emas.csv'.format(var, idx, dt_i[dt], dt_f[dt]))

		data = np.loadtxt(os.path.join('/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas/', '{0}_{1}_{2}_{3}_emas.csv'.format(var, idx, dt_i[dt], dt_f[dt])), dtype='str', delimiter=',', unpack=True)
		data_i = data[:,1:]
		data_var = np.where(data_i[1,:] == 'null', -999., data_i[1,:])	
		data_var_i = np.array(data_var, dtype=float)
		station_value.append(data_var_i)

	data_values = np.concatenate(station_value)
	
	print(data_values)
	exit()
			
	data_dates  = []
	for i in range(len(dt_h)):
			
		data_dates.append('{0}'.format(dt_h[i]))
		
	nc_output = '/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas_nc/precip_{0}_H_2018-01-01_2021-12-31.nc'.format(idx)

	print('Create netcdf')
	# Create netcdf
	ds = Dataset(nc_output, mode='w', format='NETCDF4_CLASSIC')

	ds.Conventions 	= 'CF-1.6'
	ds.title 		= 'Argentina weather station.'
	ds.institution 	= 'EMAS.'
	ds.source 		= 'Automatic weather station.'
	ds.history 		= 'Rewrote via python script.'
	ds.references 	= 'https://centrales.bolsacer.org.ar/accounts/login/?next=/.'
	ds.comment 		= 'This script convert .csv to .nc from each Argentina weather station'
		
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


















