# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot time series to each inmet automatic station"

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations_code import codes
from dict_inmet_stations_latlon import coord
from dict_inmet_stations_name import names

idx=12
dt = 'H_2018-01-01_2021-12-31'

dict_var = {2: ['pre', 'Precipitation (mm d⁻¹)', 'tp'],
            6: ['tmp', 'Temperature (°C)', 't2m'],
            12: ['uv', 'Wind speed (m s⁻¹)', 'uv10']}	
                        
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
	if i == 268:
		continue
	if i == 287:
		continue
	
	yy = coord[i][0]
	xx = coord[i][1]

	print('Reading inmet weather station:', i, codes[i], names[i][1])
	# Reading inmet weather station	
	df = pd.read_csv(os.path.join('/home/nice/Downloads/FPS_SESA/inmet/inmet_used/', 'dados_{0}_{1}.csv'.format(codes[i], dt)), sep='[:,|_]', engine='python')

	if idx == 2:
		var_x = df.iloc[:,idx]
		clim = []
		for hr in range(0, 24):
			hr_x = np.nanmean(var_x[hr::24], axis=0)
			hr_x = hr_x*24
			clim.append(hr_x)
	else:
		var_x = df.iloc[:,idx]
		clim = []
		for hr in range(0, 24):
			hr_x = np.nanmean(var_x[hr::24], axis=0)
			clim.append(hr_x)
			
	a = clim[3:]
	b = clim[:3]
	clim_inmet = a+b

	# reading era5 reanalisis
	ds = xr.open_mfdataset('/home/nice/Downloads/FPS_SESA/era5/' + '{0}_sesa_era5_2018-2021.nc'.format(dict_var[idx][2]), combine='by_coords')

	if idx == 2:
		ds = ds.tp.sel(time=slice('2018-01-01','2021-12-31'))
		ds = ds.sel(latitude=yy,longitude=xx, method='nearest')
		var = ds.groupby('time.hour').mean('time')
		var_i = var.values
		clim_i = var_i*24
	elif idx == 6:
		ds = ds.t2m.sel(time=slice('2018-01-01','2021-12-31'))
		ds = ds.sel(latitude=yy,longitude=xx, method='nearest')
		var = ds.groupby('time.hour').mean('time')
		clim_i = var.values
	else:
		ds = ds.u10.sel(time=slice('2018-01-01','2021-12-31'))
		ds = ds.sel(latitude=yy,longitude=xx, method='nearest')
		var = ds.groupby('time.hour').mean('time')
		clim_i = var.values

	clim_i = clim_i.tolist()
	a_i = clim_i[3:]
	b_i = clim_i[:3]
	clim_era5 = a_i+b_i
	
	print('Plot figure')
	# Plot figure
	fig = plt.figure()
	time = np.arange(0.5, 24 + 0.5)
	plt.plot(time, clim_inmet, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='black', label = 'INMET')
	plt.plot(time, clim_era5, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label = 'ERA5')
	plt.title('{0} - {1} (2018-2021)'.format(codes[i], names[i][1]), fontsize=8, fontweight='bold')

	if idx == 2:
		plt.yticks(np.arange(0, 13, 1))
	elif idx == 6:
		plt.yticks(np.arange(0, 33, 3))
	else:
		plt.yticks(np.arange(0, 9, 1))
		
	plt.xticks(time, ('00Z', '', '02Z', '', '04Z', '', '06Z', '', '08Z', '', '10Z', '', '12Z', '', '14Z', '', '16Z', '', '18Z', '', '20Z', '', '22Z', '', '24Z'))
	plt.xlabel('Diurnal Cycle', fontsize=8, fontweight='bold')
	plt.ylabel('{0}'.format(dict_var[idx][1]), fontsize=8, fontweight='bold')
	plt.grid(linestyle='--')
	plt.legend()

	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Downloads/FPS_SESA/figs'
	name_out = 'pyplt_diurnal_cycle_{0}_{1}_{2}.png'.format(dict_var[idx][0], codes[i], names[i][0])
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

