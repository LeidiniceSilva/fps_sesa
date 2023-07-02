# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot maps of correlation"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_stations import urug_smn
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap


def import_inmet(var):
	
	ix = []		  
	iy = []
	corr_i = []
	corr_ii = []
	corr_iii = []
	corr_iv = []
	corr_v = []
	corr_vi = []
	
	# Select lat and lon 
	for i in range(1, 101):
		yy=inmet[i][2]
		xx=inmet[i][3]
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])

		if var == 't2m':
			# reading regcm usp 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
			d_i = d_i.groupby('time.year').mean('time')
			values_i = d_i.values
			values_i = values_i-273.15

			# reading wrf ncar 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'tas_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_ii = d_ii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
			d_ii = d_ii.groupby('time.year').mean('time')
			values_ii = d_ii.values
			values_ii = values_ii-273.15

			# reading wrf ucan 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_iii = d_iii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = d_iii.values
			values_iii = values_iii-273.15
					
			# Reading inmet 
			d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc_sesa/tmp/' + 'tmp_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_iv = d_iv.tmp.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.groupby('time.year').mean('time')
			values_iv = d_iv.values
			
			# reading era5 
			d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_v = d_v.t2m.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
			d_v = d_v.groupby('time.year').mean('time')
			values_v = d_v.values
			values_v = values_v-273.15
			
		else:
			# reading regcm usp 
			d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
			d_i = d_i.groupby('time.year').mean('time')
			values_i = d_i.values

			# reading wrf ncar 
			d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'sfcWind_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_ii = d_ii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
			d_ii = d_ii.groupby('time.year').mean('time')
			values_ii = d_ii.values

			# reading wrf ucan 
			d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_iii = d_iii.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = d_iii.values
					
			# Reading inmet 
			d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc_sesa/uv/' + 'uv_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_iv = d_iv.uv.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.groupby('time.year').mean('time')
			values_iv = d_iv.values

			# reading era5 
			d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_v = d_v.u10.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
			d_v = d_v.groupby('time.year').mean('time')
			values_v = d_v.values

		# calculate correlation regcm usp
		corr_i.append(np.corrcoef(values_i, values_iv)[0][1])
		corr_ii.append(np.corrcoef(values_i, values_v)[0][1])

		# calculate correlation wrf ncar
		corr_iii.append(np.corrcoef(values_ii, values_iv)[0][1])
		corr_iv.append(np.corrcoef(values_ii, values_v)[0][1])

		# calculate correlation wrf ucan
		corr_v.append(np.corrcoef(values_iii, values_iv)[0][1])
		corr_vi.append(np.corrcoef(values_iii, values_v)[0][1])	
						
	return iy, ix, corr_i, corr_ii, corr_iii, corr_iv, corr_v, corr_vi
		

def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-60., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
	my_map.drawmeridians(np.arange(-60.,-48.,4.), labels=[0,0,0,1], size=6, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-35.,-17.,4.), labels=[1,0,0,0], size=6, linewidth=0.5, color='black') 
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=.5)

	return my_map


var = 'uv10'

# Import dataset
lat_x, lon_x, corr_i_x, corr_ii_x, corr_iii_x, corr_iv_x, corr_v_x, corr_vi_x = import_inmet(var)			

lon_xx = lon_x 
lat_yy = lat_x

reg_usp_inmet_smn = corr_i_x 
reg_usp_reanalise = corr_ii_x

wrf_ncar_inmet_smn = corr_iii_x 
wrf_ncar_reanalise = corr_iv_x 

wrf_ucan_inmet_smn = corr_v_x 
wrf_ucan_reanalise = corr_vi_x 

# Plot figure   
fig = plt.figure(figsize=(7, 6))

if var == 't2m':
	color='PRGn'
	v_min = -0.9
	v_max = 0.9
	legend = 'Correlation of temperature'
else:
	color='PRGn'
	v_min = -0.9
	v_max = 0.9
	legend = 'Correlation of wind 10m'
font_size = 8

ax = fig.add_subplot(2, 3, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM USP (INMET)', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) WRF NCAR (INMET)', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) WRF UCAN (INMET)', loc='left', fontsize=font_size, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size) 

ax = fig.add_subplot(2, 3, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) RegCM USP (ERA5)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) WRF NCAR (ERA5)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(f) WRF UCAN (ERA5)', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)
	
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_maps_corr_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
