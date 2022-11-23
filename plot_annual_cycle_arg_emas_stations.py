# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot time series to each Argentina weather station"

import os
import csv
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_stations_arg_emas import arg_emas

dt = 'H_2018-01-01_2021-12-31'
       
# Getting the data 
for j in range(1, 88):
	
	print('Reading Argentina weather station:', j, arg_emas[j][0])	

	# Reading Argentina weather stations
	ds_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas_nc/' + 'precip_{0}_{1}.nc'.format(arg_emas[j][0], dt))
	ds_i = ds_i.precip.sel(time=slice('2018-01-01','2021-12-31'))
	ds_i = ds_i.groupby('time.month').mean('time')
	values_ds_i = ds_i.values
	clim_ds_i = values_ds_i*24
	
	# Reading ERA5 reanalisis
	ds_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'tp_sesa_era5_2018-2021.nc')
	ds_ii = ds_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
	ds_ii = ds_ii.sel(latitude=arg_emas[j][2], longitude=arg_emas[j][1], method='nearest')
	ds_ii = ds_ii.groupby('time.month').mean('time')
	values_ds_ii = ds_ii.values
	clim_ds_ii = values_ds_ii*24
			
	print('Plot figure')
	# Plot figure
	fig = plt.figure()
	time = np.arange(0.5, 12 + 0.5)
	plt.plot(time, clim_ds_i, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='black', label='AMES')
	plt.plot(time, clim_ds_ii, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label='ERA5')
	plt.title('{0} (2018-2021)'.format(arg_emas[j][0]), fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 13, 1))
	plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
	plt.xlabel('Annual Cycle', fontsize=8, fontweight='bold')
	plt.ylabel('Precipitation (mm d$\mathregular{^{-1}}$)', fontsize=8, fontweight='bold')
	plt.legend(fontsize=8)
	plt.grid()
	
	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos/FPS_SESA/figs'
	name_out = 'pyplt_annual_cycle_{0}.png'.format(arg_emas[j][0])
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

