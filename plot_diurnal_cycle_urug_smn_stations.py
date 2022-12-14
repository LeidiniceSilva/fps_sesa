# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot time series to each Uruguai automatic station"

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_stations_urug_smn import urug_smn

dt = 'H_2018-01-01_2021-12-31'

# Getting the data 
for j in range(1, 72):
	
	print('Reading Uruguai weather station:', j, urug_smn[j][0])

	# Reading uruguai weather stations
	ds_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/urug_smn/urug_smn_nc/' + 'pre_{0}_{1}.nc'.format(urug_smn[j][0], dt))
	ds_i = ds_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
	ds_i = ds_i.groupby('time.hour').mean('time')
	values_ds_i = ds_i.values
	clim_ds_i = values_ds_i*24
	
	clim_ds_i = clim_ds_i.tolist()	
	clim_ds_ia = clim_ds_i[3:]
	clim_ds_ib = clim_ds_i[:3]
	clim_inmet = clim_ds_ia+clim_ds_ib
	
	# Reading ERA5 reanalisis
	ds_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'tp_era5_sesa_hr_2018-2021.nc')
	ds_ii = ds_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
	ds_ii = ds_ii.sel(latitude=urug_smn[j][1], longitude=urug_smn[j][2], method='nearest')
	ds_ii = ds_ii.groupby('time.hour').mean('time')
	values_ds_ii = ds_ii.values
	clim_ds_ii = values_ds_ii*24

	clim_ds_ii = clim_ds_ii.tolist()	
	clim_ds_iia = clim_ds_ii[3:]
	clim_ds_iib = clim_ds_ii[:3]
	clim_era5 = clim_ds_iia+clim_ds_iib
	
	print('Plot figure')
	# Plot figure
	fig = plt.figure()
	time = np.arange(0.5, 24 + 0.5)
	plt.plot(time, clim_inmet, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='black', label = 'INMET')
	plt.plot(time, clim_era5, linewidth=1.5, linestyle='--', markersize=5, marker='.', markerfacecolor='white', color='blue', label = 'ERA5')
	plt.title('{0} (2018-2021)'.format(urug_smn[j][0]), fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 13, 1))
	plt.xticks(time, ('00Z', '', '02Z', '', '04Z', '', '06Z', '', '08Z', '', '10Z', '', '12Z', '', '14Z', '', '16Z', '', '18Z', '', '20Z', '', '22Z', ''))
	plt.xlabel('Diurnal Cycle', fontsize=8, fontweight='bold')
	plt.ylabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
	plt.grid()
	plt.legend()

	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos/FPS_SESA/figs'
	name_out = 'pyplt_diurnal_cycle_{0}.png'.format(urug_smn[j][0])
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

