# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot maps of climatology"

import os
import conda
import cmocean
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_stations import urug_smn
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap


def import_inmet():
	
	ix = []		  
	iy = []
	clim_i = []
	clim_ii = []
	clim_iii = []
	clim_iv = []
	clim_v = []

	# Select lat and lon 
	for i in range(1, 101):
		yy=inmet[i][2]
		xx=inmet[i][3]
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.groupby('time.year').mean('time')
		values_i = np.nanmean(d_i.values)
		clim_i.append(values_i*86400)
		
		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		values_ii = np.nanmean(d_ii.values)
		clim_ii.append(values_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		values_iii = np.nanmean(d_iii.values)
		clim_iii.append(values_iii*86400)
				
		# Reading inmet 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc_sesa/pre/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.groupby('time.year').mean('time')
		values_iv = np.nanmean(d_iv.values)
		clim_iv.append(values_iv*24)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		values_v = np.nanmean(d_v.values)
		clim_v.append(values_v)
						
	return iy, ix, clim_i, clim_ii, clim_iii, clim_iv, clim_v
		

def import_smn():

	ix = []		  
	iy = []
	clim_i = []
	clim_ii = []
	clim_iii = []
	clim_iv = []
	clim_v = []

	# Select lat and lon 
	for i in range(1, 72):
		yy=urug_smn[i][1]
		xx=urug_smn[i][2]
		iy.append(urug_smn[i][1])
		ix.append(urug_smn[i][2])
		
		print('Reading weather station:', i, urug_smn[i][0])
		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.groupby('time.year').mean('time')
		values_i = np.nanmean(d_i.values)
		clim_i.append(values_i*86400)
		
		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		values_ii = np.nanmean(d_ii.values)
		clim_ii.append(values_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		values_iii = np.nanmean(d_iii.values)
		clim_iii.append(values_iii*86400)
				
		# Reading smn 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn/urug_smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(urug_smn[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.groupby('time.year').mean('time')
		values_iv = np.nanmean(d_iv.values)
		clim_iv.append(values_iv*24)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		values_v = np.nanmean(d_v.values)
		clim_v.append(values_v)
		
	return iy, ix, clim_i, clim_ii, clim_iii, clim_iv, clim_v


def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-60., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
	my_map.drawmeridians(np.arange(-60.,-48.,4.), labels=[0,0,0,1], size=6, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-35.,-17.,4.), labels=[1,0,0,0], size=6, linewidth=0.5, color='black') 
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=.5)

	return my_map
	

var = 'pr'

# Import dataset
lat_x, lon_x, clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x = import_inmet()			
lat_y, lon_y, clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y = import_smn()			

lon_xx = lon_x + lon_y
lat_yy = lat_x + lat_y

reg_usp = clim_i_x + clim_i_y
wrf_ncar = clim_ii_x + clim_ii_y
wrf_ucan = clim_iii_x + clim_iii_y
inmet_smn = clim_iv_x + clim_iv_y
era5 = clim_v_x + clim_v_y

# Plot figure   
fig = plt.figure(figsize=(7, 6))

color = cmocean.cm.dense
v_min = 0
v_max = 8
legend = 'Precipitation (mm d⁻¹)'
font_size = 8
	
ax = fig.add_subplot(2, 3, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) WRF NCAR', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) WRF UCAN', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='max')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)   

ax = fig.add_subplot(2, 3, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, era5, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='max')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)   

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_maps_clim_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

