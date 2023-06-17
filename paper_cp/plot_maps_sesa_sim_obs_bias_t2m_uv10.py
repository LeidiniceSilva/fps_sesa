# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "02/09/2023"
__description__ = "This script plot bias maps from regcm and database"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_sesa_inmet_stations import inmet


def import_inmet(var, dt):
	
	ix = []		  
	iy = []
	bias_i = []
	bias_ii = []

	# Select lat and lon 
	for i in range(1, 101):

		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		if var == 't2m':
			# reading regcm 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.tas.sel(time=slice('2018-06-01','2021-12-31'))
			d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_i = d_i.groupby('time.year').mean('time')
			values_i = d_i.values
			list_i = np.nanmean(values_i-273.15)
			# Reading inmet 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'tmp_{0}_{1}.nc'.format(inmet[i][0], dt))
			d_ii = d_ii.tmp.sel(time=slice('2019-01-01','2021-12-31'))
			d_ii = d_ii.groupby('time.year').mean('time')
			values_ii = d_ii.values
			list_ii = np.nanmean(values_ii)
			# reading era5 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_iii = d_iii.t2m.sel(time=slice('2019-01-01','2021-12-31'))
			d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = d_iii.values
			list_iii = np.nanmean(values_iii-273.15)
		else:
			# reading regcm 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/reg4/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.sfcWind.sel(time=slice('2018-06-01','2021-12-31'))
			d_i = d_i.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_i = d_i.groupby('time.year').mean('time')
			values_i = d_i.values
			list_i = np.nanmean(values_i)
			# Reading inmet 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/inmet/inmet_nc/' + 'uv_{0}_{1}.nc'.format(inmet[i][0], dt))
			d_ii = d_ii.uv.sel(time=slice('2019-01-01','2021-12-31'))
			d_ii = d_ii.groupby('time.year').mean('time')
			values_ii = d_ii.values
			list_ii = np.nanmean(values_ii)
			# reading era5 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/era5/' + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_iii = d_iii.u10.sel(time=slice('2019-01-01','2021-12-31'))
			d_iii = d_iii.sel(lat=inmet[i][2], lon=inmet[i][3], method='nearest')
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = d_iii.values
			list_iii = np.nanmean(values_iii)

		# calculate bias
		mean_i = list_i - list_ii
		bias_i.append(mean_i)

		mean_ii = list_i - list_iii
		bias_ii.append(mean_ii)
				
	return iy, ix, bias_i, bias_ii
		

def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-60., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
	my_map.drawmeridians(np.arange(-60.,-48.,4.), labels=[0,0,0,1], size=6, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-35.,-17.,4.), labels=[1,0,0,0], size=6, linewidth=0.5, color='black') 
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)

	return my_map


var = 'uv10'
dt = 'H_2018-01-01_2021-12-31'

print('Import dataset')
# Import dataset
iy, ix, bias_i, bias_ii = import_inmet(var, dt)			

lon_xx = ix
lat_yy = iy

bias_tot_i = bias_i
bias_tot_ii = bias_ii

print('Plot figure')
# Plot figure   
fig = plt.figure(figsize=(4, 5))

if var == 't2m':
	color='bwr'
	v_min = -3
	v_max = 3
	legend = 'Bias of temperature (°C)'
else:
	color='PiYG'
	v_min = -3
	v_max = 3
	legend = 'Bias of wind 10m (m s⁻¹)'
	
ax = fig.add_subplot(2, 1, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, bias_tot_i, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM - INMET', loc='left', fontsize=6, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=6)  

ax = fig.add_subplot(2, 1, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 4, bias_tot_ii, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) RegCM - ERA5', loc='left', fontsize=6, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=6, fontweight='bold')
cbar.ax.tick_params(labelsize=6)  
	
print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/figs_sesa'
name_out = 'pyplt_maps_bias_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
plt.show()
exit()


