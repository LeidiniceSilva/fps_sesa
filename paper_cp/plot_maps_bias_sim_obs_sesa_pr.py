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
	bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x = [], [], [], [], [], [], [], [], [], []

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
		values_i = np.nanmean(d_i.values)
		values_i = values_i*86400
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		values_ii = np.nanmean(d_ii.values)
		values_ii = values_ii*86400
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		values_iii = np.nanmean(d_iii.values)
		values_iii = values_iii*86400
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
		d_iv = d_iv.groupby('time.year').mean('time')
		values_iv = np.nanmean(d_iv.values)
		values_iv = values_iv*86400
			
		# reading wrf ucan 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		values_v = np.nanmean(d_v.values)
		values_v = values_v*86400
			
		# Reading inmet 
		d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc/pre/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		values_vi = np.nanmean(d_vi.values)
		values_vi = values_vi*24
		
		# reading era5 
		d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
		d_vii = d_vii.groupby('time.year').mean('time')
		values_vii = np.nanmean(d_vii.values)
		values_vii = values_vii

		# calculate bias regcm usp
		bias_i.append(values_i - values_vi)	
		bias_ii.append(values_i - values_vii)

		# calculate bias regcm ictp 1
		bias_iii.append(values_ii - values_vi)	
		bias_iv.append(values_ii - values_vii)

		# calculate bias regcm ictp 2
		bias_v.append(values_iii - values_vi)	
		bias_vi.append(values_iii - values_vii)

		# calculate bias wrf ncar
		bias_vii.append(values_iv - values_vi)	
		bias_viii.append(values_iv - values_vii)

		# calculate bias wrf ucan
		bias_ix.append(values_v - values_vi)	
		bias_x.append(values_v - values_vii)

	return iy, ix, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x


def import_smn_i():

	iy, ix = [], []
	bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x = [], [], [], [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 72):
		yy=smn_i[i][1]
		xx=smn_i[i][2]
		iy.append(smn_i[i][1])
		ix.append(smn_i[i][2])
		
		print('Reading weather station:', i, smn_ii[i][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.groupby('time.year').mean('time')
		values_i = np.nanmean(d_i.values)
		values_i = values_i*86400
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		values_ii = np.nanmean(d_ii.values)
		values_ii = values_ii*86400
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		values_iii = np.nanmean(d_iii.values)
		values_iii = values_iii*86400
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
		d_iv = d_iv.groupby('time.year').mean('time')
		values_iv = np.nanmean(d_iv.values)
		values_iv = values_iv*86400
			
		# reading wrf ucan 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		values_v = np.nanmean(d_v.values)
		values_v = values_v*86400
			
		# Reading smn 
		d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn_i/smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		values_vi = np.nanmean(d_vi.values)
		values_vi = values_vi*24
		
		# reading era5 
		d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
		d_vii = d_vii.groupby('time.year').mean('time')
		values_vii = np.nanmean(d_vii.values)

		# calculate bias regcm usp
		bias_i.append(values_i - values_vi)	
		bias_ii.append(values_i - values_vii)

		# calculate bias regcm ictp 1
		bias_iii.append(values_ii - values_vi)	
		bias_iv.append(values_ii - values_vii)

		# calculate bias regcm ictp 2
		bias_v.append(values_iii - values_vi)	
		bias_vi.append(values_iii - values_vii)

		# calculate bias wrf ncar
		bias_vii.append(values_iv - values_vi)	
		bias_viii.append(values_iv - values_vii)

		# calculate bias wrf ucan
		bias_ix.append(values_v - values_vi)	
		bias_x.append(values_v - values_vii)
			
	return iy, ix, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x
	

def import_smn_ii():
	
	iy, ix = [], []
	bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x = [], [], [], [], [], [], [], [], [], []

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
		values_i = np.nanmean(d_i.values)
		values_i = values_i*86400
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.groupby('time.year').mean('time')
		values_ii = np.nanmean(d_ii.values)
		values_ii = values_ii*86400
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.groupby('time.year').mean('time')
		values_iii = np.nanmean(d_iii.values)
		values_iii = values_iii*86400
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
		d_iv = d_iv.groupby('time.year').mean('time')
		values_iv = np.nanmean(d_iv.values)
		values_iv = values_iv*86400
			
		# reading wrf ucan 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.groupby('time.year').mean('time')
		values_v = np.nanmean(d_v.values)
		values_v = values_v*86400
			
		# Reading smn 
		d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn_ii/smn_nc/' + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.year').mean('time')
		values_vi = np.nanmean(d_vi.values)
		values_vi = values_vi
		
		# reading era5 
		d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
		d_vii = d_vii.groupby('time.year').mean('time')
		values_vii = np.nanmean(d_vii.values)

		# calculate bias regcm usp
		bias_i.append(values_i - values_vi)	
		bias_ii.append(values_i - values_vii)

		# calculate bias regcm ictp 1
		bias_iii.append(values_ii - values_vi)	
		bias_iv.append(values_ii - values_vii)

		# calculate bias regcm ictp 2
		bias_v.append(values_iii - values_vi)	
		bias_vi.append(values_iii - values_vii)

		# calculate bias wrf ncar
		bias_vii.append(values_iv - values_vi)	
		bias_viii.append(values_iv - values_vii)

		# calculate bias wrf ucan
		bias_ix.append(values_v - values_vi)	
		bias_x.append(values_v - values_vii)
			
	return iy, ix, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x


def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-70., llcrnrlat=-40., urcrnrlon=-45.,urcrnrlat=-15., resolution='c')
	my_map.drawmeridians(np.arange(-70,-45,5), labels=[0,0,0,1], size=8, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-40,-15,5), labels=[1,0,0,0], size=8, linewidth=0.5, color='black')
	my_map.readshapefile('/home/nice/Documentos/github_projects/shp/shp_america_sul/america_sul', 'america_sul', drawbounds=True, color='black', linewidth=.5)

	return my_map
	

var = 'pr'

# Import dataset
lat_x, lon_x, bias_i_x, bias_ii_x, bias_iii_x, bias_iv_x, bias_v_x, bias_vi_x, bias_vii_x, bias_viii_x, bias_ix_x, bias_x_x = import_inmet()			
lat_y, lon_y, bias_i_y, bias_ii_y, bias_iii_y, bias_iv_y, bias_v_y, bias_vi_y, bias_vii_y, bias_viii_y, bias_ix_y, bias_x_y = import_smn_i()			
lat_z, lon_z, bias_i_z, bias_ii_z, bias_iii_z, bias_iv_z, bias_v_z, bias_vi_z, bias_vii_z, bias_viii_z, bias_ix_z, bias_x_z = import_smn_ii()			

lat_yy = lat_x + lat_y + lat_z
lon_xx = lon_x + lon_y + lon_z

reg_usp_inmet_smn = bias_i_x + bias_i_y + bias_i_z
reg_usp_reanalise = bias_ii_x + bias_ii_y + bias_ii_z

reg_ictp_i_inmet_smn = bias_iii_x + bias_iii_y + bias_iii_z
reg_ictp_i_reanalise = bias_iv_x + bias_iv_y + bias_iv_z

reg_ictp_ii_inmet_smn = bias_v_x + bias_v_y + bias_v_z
reg_ictp_ii_reanalise = bias_vi_x + bias_vi_y + bias_vi_z

wrf_ncar_inmet_smn = bias_vii_x + bias_vii_y + bias_vii_z
wrf_ncar_reanalise = bias_viii_x + bias_viii_y + bias_viii_z

wrf_ucan_inmet_smn = bias_ix_x + bias_ix_y + bias_ix_z
wrf_ucan_reanalise = bias_x_x + bias_x_y + bias_x_z

# Plot figure   
fig = plt.figure(figsize=(11, 5))

color = 'BrBG'
v_min = -3
v_max = 3
legend = 'Bias of precipitation (mm d⁻¹)'
font_size = 8
	
ax = fig.add_subplot(2, 5, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) RegCM4 - INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.25, 0.010, 0.50]), extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(2, 5, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_i_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) RegCM5 Holtslag - \nINMET+SMN', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_ii_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) RegCM5 UW-PBL - \nINMET+SMN', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) WRF415 - INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) WRF433 - INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(f) RegCM4 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 7)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_i_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(g) RegCM5 Holtslag - \nERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 8)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_ii_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(h) RegCM5 UW-PBL - \nERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 9)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(i) WRF415 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 10)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(j) WRF433 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_maps_bias_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
