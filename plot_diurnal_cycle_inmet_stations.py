# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot time series to each INMET automatic station"

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_stations_inmet import inmet

dt = 'H_2018-01-01_2021-12-31'

# Getting the data 
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
	if j == 246:
		continue
	if j == 268:
		continue
	if j == 287:
		continue

	print('Reading INMET weather station:', j, inmet[j][0], inmet[j][1])

	# Reading inmet weather station	
	ds_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/' + 'uv_{0}_{1}.nc'.format(inmet[j][0], dt))
	ds_i = ds_i.uv.sel(time=slice('2018-01-01','2021-12-31'))
	ds_i = ds_i.groupby('time.hour').mean('time')
	values_ds_i = ds_i.values
	clim_ds_i = values_ds_i
	
	clim_ds_i = clim_ds_i.tolist()	
	clim_ds_ia = clim_ds_i[3:]
	clim_ds_ib = clim_ds_i[:3]
	clim_inmet = clim_ds_ia+clim_ds_ib
	
	# reading era5 reanalisis
	ds_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'uv_era5_sesa_hr_2018-2021.nc')
	ds_ii = ds_ii.u10.sel(time=slice('2018-01-01','2021-12-31'))
	ds_ii = ds_ii.sel(latitude=inmet[j][2], longitude=inmet[j][3], method='nearest')
	ds_ii = ds_ii.groupby('time.hour').mean('time')
	values_ds_ii = ds_ii.values
	clim_ds_ii = values_ds_ii

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
	plt.title('{0} - {1} (2018-2021)'.format(inmet[j][0], inmet[j][1]), fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 10, 1))
	plt.xticks(time, ('00Z', '', '02Z', '', '04Z', '', '06Z', '', '08Z', '', '10Z', '', '12Z', '', '14Z', '', '16Z', '', '18Z', '', '20Z', '', '22Z', ''))
	plt.xlabel('Diurnal Cycle', fontsize=8, fontweight='bold')
	plt.ylabel('Wind speed (m s⁻¹)', fontsize=8, fontweight='bold')
	plt.grid()
	plt.legend()

	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos'
	name_out = 'pyplt_diurnal_cycle_uv_{0}_{1}.png'.format(inmet[j][0], inmet[j][1])
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

