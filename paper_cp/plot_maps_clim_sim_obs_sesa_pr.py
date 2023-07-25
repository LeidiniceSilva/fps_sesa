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
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap


def import_inmet():

	iy, ix = [], []
	mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii = [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 100):
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
		d_i = np.nanmean(d_i.values)
		mean_i.append(d_i*86400)
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		d_ii = np.nanmean(d_ii.values)
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		d_iii = np.nanmean(d_iii.values)
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
		d_iv = d_iv.groupby('time.year').mean('time')
		d_iv = np.nanmean(d_iv.values)
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		d_v = np.nanmean(d_v.values)
		mean_v.append(d_v*86400)
		
		# Reading inmet 
		d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc/pre/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		d_vi = np.nanmean(d_vi.values)
		mean_vi.append(d_vi*24)
		
		# reading era5 
		d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
		d_vii = d_vii.groupby('time.year').mean('time')
		d_vii = np.nanmean(d_vii.values)
		mean_vii.append(d_vii)
				
	return iy, ix, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def import_smn_i():

	iy, ix = [], []
	mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii = [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 72):
		yy=smn_i[i][1]
		xx=smn_i[i][2]
		iy.append(smn_i[i][1])
		ix.append(smn_i[i][2])
		
		print('Reading weather station:', i, smn_i[i][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.groupby('time.year').mean('time')
		d_i = np.nanmean(d_i.values)
		mean_i.append(d_i*86400)
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		d_ii = np.nanmean(d_ii.values)
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		d_iii = np.nanmean(d_iii.values)
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
		d_iv = d_iv.groupby('time.year').mean('time')
		d_iv = np.nanmean(d_iv.values)
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		d_v = np.nanmean(d_v.values)
		mean_v.append(d_v*86400)
		
		# Reading smn 
		d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn_i/smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		d_vi = np.nanmean(d_vi.values)
		mean_vi.append(d_vi*24)
		
		# reading era5 
		d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
		d_vii = d_vii.groupby('time.year').mean('time')
		d_vii = np.nanmean(d_vii.values)
		mean_vii.append(d_vii)
				
	return iy, ix, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def import_smn_ii():
	
	iy, ix = [], []
	mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii = [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 86):
		yy=smn_ii[i][1]
		xx=smn_ii[i][2]
		iy.append(smn_ii[i][1])
		ix.append(smn_ii[i][2])
		
		print('Reading weather station:', i, smn_ii[i][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.groupby('time.year').mean('time')
		d_i = np.nanmean(d_i.values)
		mean_i.append(d_i*86400)
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		d_ii = np.nanmean(d_ii.values)
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		d_iii = np.nanmean(d_iii.values)
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
		d_iv = d_iv.groupby('time.year').mean('time')
		d_iv = np.nanmean(d_iv.values)
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		d_v = np.nanmean(d_v.values)
		mean_v.append(d_v*86400)
		
		# Reading smn 
		d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn_ii/smn_nc/' + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		d_vi = np.nanmean(d_vi.values)
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
		d_vii = d_vii.groupby('time.year').mean('time')
		d_vii = np.nanmean(d_vii.values)
		mean_vii.append(d_vii)
				
	return iy, ix, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-70., llcrnrlat=-40., urcrnrlon=-45.,urcrnrlat=-15., resolution='c')
	my_map.drawmeridians(np.arange(-70,-45,5), labels=[0,0,0,1], size=8, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-40,-15,5), labels=[1,0,0,0], size=8, linewidth=0.5, color='black')
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=.5)

	return my_map
	

var = 'pr'

# Import dataset
lat_x, lon_x, clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x, clim_vi_x, clim_vii_x = import_inmet()			
lat_y, lon_y, clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y, clim_vi_y, clim_vii_y = import_smn_i()			
lat_z, lon_z, clim_i_z, clim_ii_z, clim_iii_z, clim_iv_z, clim_v_z, clim_vi_z, clim_vii_z = import_smn_ii()			

lat_yy = lat_x + lat_y + lat_z
lon_xx = lon_x + lon_y + lon_z

reg_usp = clim_i_x + clim_i_y + clim_i_z
reg_ictp_i = clim_ii_x + clim_ii_y + clim_ii_z
reg_ictp_ii = clim_iii_x + clim_iii_y + clim_iii_z
wrf_ncar = clim_iv_x + clim_iv_y + clim_iv_z
wrf_ucan = clim_v_x + clim_v_y + clim_v_z
inmet_smn = clim_vi_x + clim_vi_y + clim_vi_z
era5 = clim_vii_x + clim_vi_y + clim_vii_z

# Plot figure   
fig = plt.figure(figsize=(7, 7))

color = cmocean.cm.dense
v_min = 0
v_max = 8
legend = 'Precipitation (mm d⁻¹)'
font_size = 8
	
ax = fig.add_subplot(3, 3, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM4', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_i, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) RegCM5 Holtslag', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_ii, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) RegCM5 UW-PBL', loc='left', fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.25, 0.019, 0.50]), extend='max')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(3, 3, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) WRF NCAR', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) WRF UCAN', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(f) INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 3, 7)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, era5, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(g) ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold') 

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_maps_clim_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

