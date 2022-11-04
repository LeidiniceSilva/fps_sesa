# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot time series to each inmet automatic station"

import os
import csv
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_stations_urug_smn import urug_smn

csv_header = []

for idx in range(0, 4):

	csv_header.append(urug_smn[idx+1][0])
	
print(csv_header)
	
# Reading uruguai weather stations
df = pd.read_csv(os.path.join('/home/nice/Documentos/FPS_SESA/uru/', 'PP_URUG_SMN_2018_2021.csv'))
df['Date'] = pd.to_datetime(df['Date'])
df = df.set_index('Date').resample('H').sum()
var = df.values
	
with open('example.csv', 'w') as file:
	writer = csv.writer(file)
	writer.writerow(var)

exit()
		
# ~ for idx in range(0, 72):
	
	# ~ # Reading uruguai weather stations
	# ~ df = pd.read_csv(os.path.join('/home/nice/Documentos/FPS_SESA/uru/', 'PP_URUG_SMN_2018_2021.csv'))
	# ~ df['Date'] = pd.to_datetime(df['Date'])
	# ~ df = df.set_index('Date').resample('H').sum()
	# ~ var = df_i.iloc[:,idx]
	
	# ~ bins = [ 1,2,3,4,5 ]
	# ~ freq = [ 9,8,7,6,5 ]

	# ~ f = open("test.csv", "w")

	# ~ for i in xrange(len(bins)):
		# ~ f.write("{} {}\n".format(bins[i], freq[i]))
	# ~ f.close()
	
	# ~ with open('example.csv', 'w') as file:
		# ~ writer = csv.writer(file)
		# ~ writer.writerow(df_h)

	# ~ exit()



df['Data Medicao'] = pd.to_datetime(df['Data Medicao'], format='%Y-%m-%d %H:%M:%S', errors='ignore')
df_i = df.groupby(pd.Grouper(key='Data Medicao', freq='M')).mean()


idx=1
dt = 'H_2018-01-01_2021-12-31'

dict_var = {1: ['pre', 'Precipitation (mm d$\mathregular{^{-1}}$)', 'tp'],
            5: ['tmp', 'Temperature (°C)', 't2m'],
            11: ['uv', 'Wind speed (m s$\mathregular{^{-1}}$)', 'uv10']}	
            
# Getting the data 
for i in range(1, 289):
	
	if i == 4:
		continue
	if i == 12:
		continue
	if i == 45:
		continue
	if i == 55:
		continue
	if i == 77:
		continue
	if i == 98:
		continue
	if i == 99:
		continue
	if i == 118:
		continue
	if i == 122:
		continue
	if i == 130:
		continue
	if i == 135:
		continue
	if i == 151:
		continue
	if i == 155:
		continue
	if i == 159:
		continue
	if i == 160:
		continue
	if i == 163:
		continue
	if i == 164:
		continue
	if i == 181:
		continue
	if i == 183:
		continue
	if i == 186:
		continue
	if i == 187:
		continue
	if i == 188:
		continue
	if i == 209:
		continue
	if i == 216:
		continue
	if i == 228:
		continue
	if i == 236:
		continue
	if i == 246:
		continue
	if i == 268:
		continue
	if i == 287:
		continue
	
	yy = coord[i][0]
	xx = coord[i][1]

	print('Reading inmet weather station:', i, codes[i], names[i][1])
	# Reading inmet weather station	
	df = pd.read_csv(os.path.join('/home/nice/Documentos/FPS_SESA/inmet/inmet_used/', 'dados_{0}_{1}.csv'.format(codes[i], dt)), sep='[:,|_]', engine='python')
	df['Data Medicao'] = pd.to_datetime(df['Data Medicao'], format='%Y-%m-%d %H:%M:%S', errors='ignore')
	df_i = df.groupby(pd.Grouper(key='Data Medicao', freq='M')).mean()

	if idx == 1:
		var_x = df_i.iloc[:,idx]
		clim = []
		for mon in range(0, 12):
			mon_x = np.nanmean(var_x[mon::12], axis=0)
			mon_x = mon_x*24
			clim.append(mon_x)
	else:
		var_x = df_i.iloc[:,idx]
		clim = []
		for mon in range(0, 12):
			mon_x = np.nanmean(var_x[mon::12], axis=0)
			clim.append(mon_x)

	# reading era5 reanalisis
	ds = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + '{0}_sesa_era5_2018-2021.nc'.format(dict_var[idx][2]))

	if idx == 1:
		ds = ds.tp.sel(time=slice('2018-01-01','2021-12-31'))
		ds = ds.sel(latitude=yy,longitude=xx, method='nearest')
		var = ds.groupby('time.month').mean('time')
		var_i = var.values
		clim_i = var_i*24
	elif idx == 5:
		ds = ds.t2m.sel(time=slice('2018-01-01','2021-12-31'))
		ds = ds.sel(latitude=yy,longitude=xx, method='nearest')
		var = ds.groupby('time.month').mean('time')
		clim_i = var.values
	else:
		ds = ds.u10.sel(time=slice('2018-01-01','2021-12-31'))
		ds = ds.sel(latitude=yy,longitude=xx, method='nearest')
		var = ds.groupby('time.month').mean('time')
		clim_i = var.values
			
	print('Plot figure')
	# Plot figure
	fig = plt.figure()
	time = np.arange(0.5, 12 + 0.5)
	plt.plot(time, clim, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='black', label = 'INMET')
	plt.plot(time, clim_i, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label = 'ERA5')
	plt.title('{0} - {1} (2018-2021)'.format(codes[i], names[i][1]), fontsize=8, fontweight='bold')
	
	if idx == 1:
		plt.yticks(np.arange(0, 13, 1))
	elif idx == 5:
		plt.yticks(np.arange(0, 33, 3))
	else:
		plt.yticks(np.arange(0, 9, 1))
	
	plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
	plt.xlabel('Annual Cycle', fontsize=8, fontweight='bold')
	plt.ylabel('{0}'.format(dict_var[idx][1]), fontsize=8, fontweight='bold')
	plt.legend(fontsize=8)
	plt.grid()
	
	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos/FPS_SESA/figs'
	name_out = 'pyplt_annual_cycle_{0}_{1}_{2}.png'.format(dict_var[idx][0], codes[i], names[i][0])
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

