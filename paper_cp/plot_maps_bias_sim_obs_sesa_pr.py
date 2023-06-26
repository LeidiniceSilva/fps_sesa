# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot maps of bias"

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


def import_inmet():
	
	ix = []		  
	iy = []
	bias_i = []
	bias_ii = []
	bias_iii = []
	bias_iv = []
	bias_v = []
	bias_vi = []

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
		values_i = d_i.values
		values_i = np.nanmean(values_i*86400)

		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		values_ii = d_ii.values
		values_ii = np.nanmean(values_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		values_iii = d_iii.values
		values_iii = np.nanmean(values_iii*86400)
				
		# Reading inmet 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc_sesa/pre/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.groupby('time.year').mean('time')
		values_iv = d_iv.values
		values_iv = np.nanmean(values_iv*24)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		values_v = d_v.values
		values_v = np.nanmean(values_v)

		# calculate bias regcm usp
		bias_i.append(values_i - values_iv)	
		bias_ii.append(values_i - values_v)

		# calculate bias wrf ncar
		bias_iii.append(values_ii - values_iv)	
		bias_iv.append(values_ii - values_v)

		# calculate bias wrf ucan
		bias_v.append(values_iii - values_iv)	
		bias_vi.append(values_iii - values_v)
						
	return iy, ix, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi
		

def import_smn():

	ix = []		  
	iy = []
	bias_i = []
	bias_ii = []
	bias_iii = []
	bias_iv = []
	bias_v = []
	bias_vi = []

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
		values_i = d_i.values
		values_i = np.nanmean(values_i*86400)

		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		values_ii = d_ii.values
		values_ii = np.nanmean(values_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		values_iii = d_iii.values
		values_iii = np.nanmean(values_iii*86400)
				
		# Reading smn 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn/urug_smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(urug_smn[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.groupby('time.year').mean('time')
		values_iv = d_iv.values
		values_iv = np.nanmean(values_iv*24)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		values_v = d_v.values
		values_v = np.nanmean(values_v)

		# calculate bias regcm usp
		bias_i.append(values_i - values_iv)	
		bias_ii.append(values_i - values_v)

		# calculate bias wrf ncar
		bias_iii.append(values_ii - values_iv)	
		bias_iv.append(values_ii - values_v)

		# calculate bias wrf ucan
		bias_v.append(values_iii - values_iv)	
		bias_vi.append(values_iii - values_v)
		
	return iy, ix, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi
	
	
def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-60., llcrnrlat=-35., urcrnrlon=-48.,urcrnrlat=-17., resolution='c')
	my_map.drawmeridians(np.arange(-60.,-48.,4.), labels=[0,0,0,1], size=6, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-35.,-17.,4.), labels=[1,0,0,0], size=6, linewidth=0.5, color='black') 
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=.5)

	return my_map


var = 'pr'

# Import dataset
lat_x, lon_x, bias_i_x, bias_ii_x, bias_iii_x, bias_iv_x, bias_v_x, bias_vi_x = import_inmet()			
lat_y, lon_y, bias_i_y, bias_ii_y, bias_iii_y, bias_iv_y, bias_v_y, bias_vi_y = import_smn()			

lon_xx = lon_x + lon_y
lat_yy = lat_x + lat_y

reg_usp_inmet_smn = bias_i_x + bias_i_y
reg_usp_reanalise = bias_ii_x + bias_ii_y

wrf_ncar_inmet_smn = bias_iii_x + bias_iii_y
wrf_ncar_reanalise = bias_iv_x + bias_iv_y

wrf_ucan_inmet_smn = bias_v_x + bias_v_y
wrf_ucan_reanalise = bias_vi_x + bias_vi_y

# Plot figure   
fig = plt.figure(figsize=(7, 6))

color = 'BrBG'
v_min = -3
v_max = 3
legend = 'Bias of precipitation (mm d⁻¹)'
font_size = 8

ax = fig.add_subplot(2, 3, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM USP - INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) WRF NCAR - INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) WRF UCAN - INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size) 

ax = fig.add_subplot(2, 3, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) RegCM USP - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) WRF NCAR - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 3, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(f) WRF UCAN - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
cbar = my_map.colorbar(shrink=0.8, pad=0.1, extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_maps_bias_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

