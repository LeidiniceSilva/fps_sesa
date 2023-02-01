# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot annual cycle for subregions on sesa"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import itertools

from dict_stations_inmet import inmet
from dict_stations_arg_emas import arg_emas
from dict_stations_urug_smn import urug_smn

subregion_i = []
subregion_ii = []
subregion_iii = []
subregion_iv = []
subregion_v = []
subregion_vi = []
subregion_vii = []
subregion_viii = []
subregion_ix = []

era5_i = []
era5_ii = []
era5_iii = []
era5_iv = []
era5_v = []
era5_vi = []
era5_vii = []
era5_viii = []
era5_ix = []

# Getting the data subregion I
for j in range(1, 40):
	
	if j == 4:
		continue
	if j == 12:
		continue

	# Reading inmet weather station	
	ds_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[j][0]))
	ds_i = ds_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
	ds_i = ds_i.groupby('time.month').mean('time')
	subregion_i.append(ds_i.values*24)

	# reading era5 reanalisis
	ds_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'tp_sesa_era5_2018-2021.nc')
	ds_ii = ds_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
	ds_ii = ds_ii.sel(latitude=inmet[j][2], longitude=inmet[j][3], method='nearest')
	ds_ii = ds_ii.groupby('time.month').mean('time')
	era5_i.append(ds_ii.values*24)

# Getting the data subregion II
for j in range(40, 80):
	
	if j == 45:
		continue
	if j == 55:
		continue
	if j == 77:
		continue

	# Reading inmet weather station	
	ds_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[j][0]))
	ds_i = ds_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
	ds_i = ds_i.groupby('time.month').mean('time')
	subregion_ii.append(ds_i.values*24)

	# reading era5 reanalisis
	ds_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'tp_sesa_era5_2018-2021.nc')
	ds_ii = ds_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
	ds_ii = ds_ii.sel(latitude=inmet[j][2], longitude=inmet[j][3], method='nearest')
	ds_ii = ds_ii.groupby('time.month').mean('time')
	era5_ii.append(ds_ii.values*24)

# Getting the data subregion III
for j in range(80, 120):
	
	if j == 98:
		continue
	if j == 99:
		continue
	if j == 118:
		continue

	# Reading inmet weather station	
	ds_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/inmet/inmet_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[j][0]))
	ds_i = ds_i.pre.sel(time=slice('2018-01-01','2021-12-31'))
	ds_i = ds_i.groupby('time.month').mean('time')
	subregion_iii.append(ds_i.values*24)

	# reading era5 reanalisis
	ds_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + 'tp_sesa_era5_2018-2021.nc')
	ds_ii = ds_ii.tp.sel(time=slice('2018-01-01','2021-12-31'))
	ds_ii = ds_ii.sel(latitude=inmet[j][2], longitude=inmet[j][3], method='nearest')
	ds_ii = ds_ii.groupby('time.month').mean('time')
	era5_iii.append(ds_ii.values*24)
	
era5_subregion_i = np.nanmean(era5_i, axis=0)
era5_subregion_ii = np.nanmean(era5_ii, axis=0)
era5_subregion_iii = np.nanmean(era5_iii, axis=0)
era5_subregion_iv = np.nanmean(era5_i, axis=0)
era5_subregion_v = np.nanmean(era5_i, axis=0)
era5_subregion_vi = np.nanmean(era5_i, axis=0)
era5_subregion_vii = np.nanmean(era5_i, axis=0)
era5_subregion_viii = np.nanmean(era5_i, axis=0)
era5_subregion_ix = np.nanmean(era5_i, axis=0)

print('Plot figure')
# Plot figure
fig = plt.figure(figsize=(11, 9))
time = np.arange(0.5, 12 + 0.5)

ax = fig.add_subplot(3, 3, 1)  
for i in range(0, len(subregion_i)):
	plt.plot(time, subregion_i[i])
	plt.plot(time, era5_subregion_i, color='k', marker='.')
	plt.title(u'(a) Subregion I', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Precipitation (mm d$\mathregular{^{-1}}$)', fontsize=8, fontweight='bold')
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)

ax = fig.add_subplot(3, 3, 2)  
for i in range(0, len(subregion_ii)):
	plt.plot(time, subregion_ii[i])
	plt.plot(time, era5_subregion_ii, color='k', marker='.')
	plt.title(u'(b) Subregion II', loc='left', fontsize=8, fontweight='bold')
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)

ax = fig.add_subplot(3, 3, 3)  
for i in range(0, len(subregion_iii)):
	plt.plot(time, subregion_iii[i])
	plt.plot(time, era5_subregion_iii, color='k', marker='.')
	plt.title(u'(a) Subregion III', loc='left', fontsize=8, fontweight='bold')
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)

ax = fig.add_subplot(3, 3, 4)  
for i in range(0, len(subregion_i)):
	plt.plot(time, subregion_i[i])
	plt.plot(time, era5_subregion_iv, color='k', marker='.')
	plt.title(u'(c) Subregion VI', loc='left', fontsize=8, fontweight='bold')
	plt.ylabel(u'Precipitation (mm d$\mathregular{^{-1}}$)', fontsize=8, fontweight='bold')
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)

ax = fig.add_subplot(3, 3, 5)  
for i in range(0, len(subregion_i)):
	plt.plot(time, subregion_i[i])
	plt.plot(time, era5_subregion_v, color='k', marker='.')
	plt.title(u'(d) Subregion V', loc='left', fontsize=8, fontweight='bold')
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)

ax = fig.add_subplot(3, 3, 6)  
for i in range(0, len(subregion_i)):
	plt.plot(time, subregion_i[i])
	plt.plot(time, era5_subregion_vi, color='k', marker='.')
	plt.title(u'(e) Subregion VI', loc='left', fontsize=8, fontweight='bold')
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)

ax = fig.add_subplot(3, 3, 7)  
for i in range(0, len(subregion_i)):
	plt.plot(time, subregion_i[i])
	plt.plot(time, era5_subregion_vii, color='k', marker='.')
	plt.title(u'(f) Subregion VII', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Months', fontsize=8, fontweight='bold')
	plt.ylabel(u'Precipitation (mm d$\mathregular{^{-1}}$)', fontsize=8, fontweight='bold')
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)

ax = fig.add_subplot(3, 3, 8)  
for i in range(0, len(subregion_i)):
	plt.plot(time, subregion_i[i])
	plt.plot(time, era5_subregion_viii, color='k', marker='.')
	plt.title(u'(g) Subregion VIII', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Months', fontsize=8, fontweight='bold')
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)

ax = fig.add_subplot(3, 3, 9)  
for i in range(0, len(subregion_i)):
	plt.plot(time, subregion_i[i])
	plt.plot(time, era5_subregion_ix, color='k', marker='.')
	plt.title(u'(h) Subregion IX', loc='left', fontsize=8, fontweight='bold')
	plt.xlabel(u'Months', fontsize=8, fontweight='bold')
	plt.xticks(time, ('J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D'), fontsize=8)
	plt.yticks(np.arange(0, 16, 2), fontsize=8)
	plt.ylim(0, 14)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_annual_cycle_all_subregions_2018-2021.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()
