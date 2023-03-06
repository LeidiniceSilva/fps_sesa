# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot annual cycle of the INMET weather station"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_sesa_inmet_stations import inmet


def import_inmet(var, dt):
	
	mean_i = []
	mean_ii = []
	mean_iii = []

	# Select lat and lon 
	for i in range(1, 155):
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])

		if var == 't2m':
			# reading regcm 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.tas.sel(time=slice('2019-01-01','2021-12-31'))
			d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_i = d_i.groupby('time.month').mean('time')
			values_i = d_i.values
			mean_i.append(values_i-273.15)
			# Reading inmet 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'tmp_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_ii = d_ii.tmp.sel(time=slice('2019-01-01','2021-12-31'))
			d_ii = d_ii.groupby('time.month').mean('time')
			values_ii = d_ii.values
			mean_ii.append(values_ii)	
			# reading era5 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_iii = d_iii.t2m.sel(time=slice('2019-01-01','2021-12-31'))
			d_iii = d_iii.groupby('time.month').mean('time')
			d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			values_iii = d_iii.values
			mean_iii.append(values_iii-273.15)
		else:
			# reading regcm 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.sfcWind.sel(time=slice('2019-01-01','2021-12-31'))
			d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_i = d_i.groupby('time.month').mean('time')
			values_i = d_i.values
			mean_i.append(values_i)
			# Reading inmet 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'uv_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_ii = d_ii.uv.sel(time=slice('2019-01-01','2021-12-31'))
			d_ii = d_ii.groupby('time.month').mean('time')
			values_ii = d_ii.values
			mean_ii.append(values_ii)	
			# reading era5 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_iii = d_iii.u10.sel(time=slice('2019-01-01','2021-12-31'))
			d_iii = d_iii.groupby('time.month').mean('time')
			d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			values_iii = d_iii.values
			mean_iii.append(values_iii)
				
	return mean_i, mean_ii, mean_iii


var = 't2m'
dt = 'H_2018-01-01_2021-12-31'

print('Import dataset')
# Import dataset
mean_i, mean_ii, mean_iii = import_inmet(var, dt)			
		
clim_regcm = mean_i
clim_inmet = mean_ii
clim_era5 = mean_iii

list_hc = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 0,
3, 1, 3, 1, 1, 1, 3, 1, 3, 4, 3, 1, 3, 1, 3, 3, 3, 3, 1, 3, 0, 3, 3, 3, 3,
3, 4, 3, 4, 0, 0, 0, 4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 0, 4, 4, 4, 4, 4, 3, 4,
3, 0, 4, 0, 4, 4, 4, 4, 4, 0, 4, 0, 0, 4, 0, 4, 4, 4, 4, 4, 4, 0, 4, 0, 2,
4, 4, 2, 4, 4, 2, 2, 2, 2, 4, 2, 2, 4, 2, 2, 2, 2, 0, 0, 0, 2, 0, 2, 2, 2, 
0, 2, 2, 2, 2, 2, 4, 0]

count_i = []
count_ii = []
count_iii = []
count_iv = []
count_v = []

for count, idx in enumerate(list_hc):
	if idx == 0:
		count_i.append(count)
	if idx == 1:
		count_ii.append(count)
	if idx == 2:
		count_iii.append(count)
	if idx == 3:
		count_iv.append(count)
	if idx == 4:
		count_v.append(count)

regcm_i = []
regcm_ii = []
regcm_iii = []
regcm_iv = []
regcm_v = []

inmet_i = []
inmet_ii = []
inmet_iii = []
inmet_iv = []
inmet_v = []

era5_i = []
era5_ii = []
era5_iii = []
era5_iv = []
era5_v = []

for c_i in count_i:
	regcm_i.append(clim_regcm[c_i])
	inmet_i.append(clim_inmet[c_i])
	era5_i.append(clim_era5[c_i])

for c_ii in count_ii:
	regcm_ii.append(clim_regcm[c_ii])
	inmet_ii.append(clim_inmet[c_ii])
	era5_ii.append(clim_era5[c_ii])
	
for c_iii in count_iii:
	regcm_iii.append(clim_regcm[c_iii])
	inmet_iii.append(clim_inmet[c_iii])
	era5_iii.append(clim_era5[c_iii])
	
for c_iv in count_iv:
	regcm_iv.append(clim_regcm[c_iv])
	inmet_iv.append(clim_inmet[c_iv])
	era5_iv.append(clim_era5[c_iv])
	
for c_v in count_v:
	regcm_v.append(clim_regcm[c_v])
	inmet_v.append(clim_inmet[c_v])
	era5_v.append(clim_era5[c_v])

regcm_cluster_i = np.nanmean(regcm_i, axis=0)
inmet_cluster_i = np.nanmean(inmet_i, axis=0)
era5_cluster_i = np.nanmean(era5_i, axis=0)

regcm_cluster_ii = np.nanmean(regcm_ii, axis=0)
inmet_cluster_ii = np.nanmean(inmet_ii, axis=0)
era5_cluster_ii = np.nanmean(era5_ii, axis=0)

regcm_cluster_iii = np.nanmean(regcm_iii, axis=0)
inmet_cluster_iii = np.nanmean(inmet_iii, axis=0)
era5_cluster_iii = np.nanmean(era5_iii, axis=0)

regcm_cluster_iv = np.nanmean(regcm_iv, axis=0)
inmet_cluster_iv = np.nanmean(inmet_iv, axis=0)
era5_cluster_iv = np.nanmean(era5_iv, axis=0)

regcm_cluster_v = np.nanmean(regcm_v, axis=0)
inmet_cluster_v = np.nanmean(inmet_v, axis=0)
era5_cluster_v = np.nanmean(era5_v, axis=0)

print('Plot figure')
# Plot figure
fig = plt.figure(figsize=(10, 8))
time = np.arange(0.5, 12 + 0.5)

if var == 't2m':
	legend = 'Temperature (°C)'
else:
	legend = 'Wind 10m (m s⁻¹)'
	
ax = fig.add_subplot(3, 2, 1)
plt.plot(time, regcm_cluster_i, linewidth=1.5, linestyle='-', markersize=5, marker='o', markerfacecolor='white', color='black', label = 'RegCM')
plt.plot(time, inmet_cluster_i, linewidth=1.5, linestyle='-', markersize=5, marker='s', markerfacecolor='white', color='blue', label = 'INMET')
plt.plot(time, era5_cluster_i, linewidth=1.5, linestyle='-', markersize=5, marker='^', markerfacecolor='white', color='red', label = 'ERA5')
plt.title('(a) Cluster I', loc='left', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
if var == 't2m':
	plt.ylim(12, 30)
	plt.yticks(np.arange(12, 32, 2))
else:
	plt.ylim(0, 5)
	plt.yticks(np.arange(0, 5.5, 0.5))
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
plt.setp(ax.get_xticklabels(), visible=False)
plt.legend(loc=9, ncol=3, frameon=False)

ax = fig.add_subplot(3, 2, 2)
plt.plot(time, regcm_cluster_ii, linewidth=1.5, linestyle='-', markersize=5, marker='o', markerfacecolor='white', color='black', label = 'RegCM')
plt.plot(time, inmet_cluster_ii, linewidth=1.5, linestyle='-', markersize=5, marker='s', markerfacecolor='white', color='blue', label = 'INMET')
plt.plot(time, era5_cluster_ii, linewidth=1.5, linestyle='-', markersize=5, marker='^', markerfacecolor='white', color='red', label = 'ERA5')
plt.title('(b) Cluster II', loc='left', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
if var == 't2m':
	plt.ylim(12, 30)
	plt.yticks(np.arange(12, 32, 2))
else:
	plt.ylim(0, 5)
	plt.yticks(np.arange(0, 5.5, 0.5))
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(3, 2, 3)
plt.plot(time, regcm_cluster_iii, linewidth=1.5, linestyle='-', markersize=5, marker='o', markerfacecolor='white', color='black', label = 'RegCM')
plt.plot(time, inmet_cluster_iii, linewidth=1.5, linestyle='-', markersize=5, marker='s', markerfacecolor='white', color='blue', label = 'INMET')
plt.plot(time, era5_cluster_iii, linewidth=1.5, linestyle='-', markersize=5, marker='^', markerfacecolor='white', color='red', label = 'ERA5')
plt.title('(c) Cluster III', loc='left', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
if var == 't2m':
	plt.ylim(12, 30)
	plt.yticks(np.arange(12, 32, 2))
else:
	plt.ylim(0, 5)
	plt.yticks(np.arange(0, 5.5, 0.5))
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(3, 2, 4)
plt.plot(time, regcm_cluster_iv, linewidth=1.5, linestyle='-', markersize=5, marker='o', markerfacecolor='white', color='black', label = 'RegCM')
plt.plot(time, inmet_cluster_iv, linewidth=1.5, linestyle='-', markersize=5, marker='s', markerfacecolor='white', color='blue', label = 'INMET')
plt.plot(time, era5_cluster_iv, linewidth=1.5, linestyle='-', markersize=5, marker='^', markerfacecolor='white', color='red', label = 'ERA5')
plt.title('(d) Cluster IV', loc='left', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
plt.xlabel('Months', fontweight='bold')
if var == 't2m':
	plt.ylim(12, 30)
	plt.yticks(np.arange(12, 32, 2))
else:
	plt.ylim(0, 5)
	plt.yticks(np.arange(0, 5.5, 0.5))
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))

ax = fig.add_subplot(3, 2, 5)
plt.plot(time, regcm_cluster_v, linewidth=1.5, linestyle='-', markersize=5, marker='o', markerfacecolor='white', color='black', label = 'RegCM')
plt.plot(time, inmet_cluster_v, linewidth=1.5, linestyle='-', markersize=5, marker='s', markerfacecolor='white', color='blue', label = 'INMET')
plt.plot(time, era5_cluster_v, linewidth=1.5, linestyle='-', markersize=5, marker='^', markerfacecolor='white', color='red', label = 'ERA5')
plt.title('(e) Cluster V', loc='left', fontweight='bold')
plt.ylabel('{0}'.format(legend), fontweight='bold')
plt.xlabel('Months', fontweight='bold')
if var == 't2m':
	plt.ylim(12, 30)
	plt.yticks(np.arange(12, 32, 2))
else:
	plt.ylim(0, 5)
	plt.yticks(np.arange(0, 5.5, 0.5))
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_stations_cluster_annual_cycle_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()
