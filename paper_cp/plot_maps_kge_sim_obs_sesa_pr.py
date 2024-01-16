# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot maps of kge function"

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap

path = '/afs/ictp.it/home/m/mda_silv/Documents'


def compute_kge(model, obs):

	"""
	The input arrays must have the same dimensions
	Param model: Numpy array with model data
	Param obs: Numpy array with obs data
	Return: Kling-Gupta Efficiency
	"""

	p1 = np.corrcoef(obs, model)[0][1]
	p2 = np.nanmean(obs)
	p3 = np.nanmean(model)
	p4 = np.nanstd(obs, ddof=0)
	p5 = np.nanstd(model, ddof=0)
	p6 = p3/p2
	p7 = p5/p4
	p8 = np.sqrt((p1 -1)**2 + (p6 -1)**2 + (p7 -1)**2)
	kge = 1 - p8

	return kge
	

def basemap():
	
	my_map = Basemap(projection='cyl', llcrnrlon=-70., llcrnrlat=-40., urcrnrlon=-45.,urcrnrlat=-15., resolution='c')
	my_map.drawmeridians(np.arange(-70,-45,5), labels=[0,0,0,1], size=font_size, linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-40,-15,5), labels=[1,0,0,0], size=font_size, linewidth=0.5, color='black')
	my_map.readshapefile('{0}/github_projects/shp/shp_america_sul/america_sul'.format(path), 'america_sul', drawbounds=True, color='black', linewidth=.5)

	return my_map
	
	
def import_inmet():

	iy, ix = [], []
	kge_i, kge_ii, kge_iii, kge_iv, kge_v, kge_vi, kge_vii, kge_viii, kge_ix, kge_x = [], [], [], [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 100):
		yy=inmet[i][2]
		xx=inmet[i][3]
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		values_i = values_i*86400
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.month').mean('time')
		values_ii = d_ii.values
		values_ii = values_ii*86400
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.month').mean('time')
		values_iii = d_iii.values
		values_iii = values_iii*86400
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.month').mean('time')
		values_iv = d_iv.values
		values_iv = values_iv*86400
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.month').mean('time')
		values_v = d_v.values
		values_v = values_v*86400
			
		# Reading inmet 
		d_vi = xr.open_dataset('{0}/FPS_SESA/database/obs/inmet/inmet_nc_sesa/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.month').mean('time')
		values_vi = d_vi.values
		values_vi = values_vi*24
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/FPS_SESA/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.groupby('time.month').mean('time')
		values_vii = d_vii.values
		values_vii = values_vii

		# calculate kge regcm usp
		kge_i.append(compute_kge(values_i, values_vi))	
		kge_ii.append(compute_kge(values_i, values_vii))

		# calculate kge regcm ictp 1
		kge_iii.append(compute_kge(values_ii, values_vi))	
		kge_iv.append(compute_kge(values_ii, values_vii))

		# calculate kge regcm ictp 2
		kge_v.append(compute_kge(values_iii, values_vi))
		kge_vi.append(compute_kge(values_iii, values_vii))

		# calculate kge wrf ncar
		kge_vii.append(compute_kge(values_iv, values_vi))	
		kge_viii.append(compute_kge(values_iv, values_vii))

		# calculate kge wrf ucan
		kge_ix.append(compute_kge(values_v, values_vi))
		kge_x.append(compute_kge(values_v, values_vii))

	return iy, ix, kge_i, kge_ii, kge_iii, kge_iv, kge_v, kge_vi, kge_vii, kge_viii, kge_ix, kge_x


def import_smn_i():

	iy, ix = [], []
	kge_i, kge_ii, kge_iii, kge_iv, kge_v, kge_vi, kge_vii, kge_viii, kge_ix, kge_x = [], [], [], [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 72):
		yy=smn_i[i][1]
		xx=smn_i[i][2]
		iy.append(smn_i[i][1])
		ix.append(smn_i[i][2])
		
		print('Reading weather station:', i, smn_i[i][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		values_i = values_i*86400
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.month').mean('time')
		values_ii = d_ii.values
		values_ii = values_ii*86400
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.month').mean('time')
		values_iii = d_iii.values
		values_iii = values_iii*86400
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.month').mean('time')
		values_iv = d_iv.values
		values_iv = values_iv*86400
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.month').mean('time')
		values_v = d_v.values
		values_v = values_v*86400
			
		# Reading smn 
		d_vi = xr.open_dataset('{0}/FPS_SESA/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.month').mean('time')
		values_vi = d_vi.values
		values_vi = values_vi*24
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/FPS_SESA/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.groupby('time.month').mean('time')
		values_vii = d_vii.values

		# calculate kge regcm usp
		kge_i.append(compute_kge(values_i, values_vi))	
		kge_ii.append(compute_kge(values_i, values_vii))

		# calculate kge regcm ictp 1
		kge_iii.append(compute_kge(values_ii, values_vi))	
		kge_iv.append(compute_kge(values_ii, values_vii))

		# calculate kge regcm ictp 2
		kge_v.append(compute_kge(values_iii, values_vi))
		kge_vi.append(compute_kge(values_iii, values_vii))

		# calculate kge wrf ncar
		kge_vii.append(compute_kge(values_iv, values_vi))	
		kge_viii.append(compute_kge(values_iv, values_vii))

		# calculate kge wrf ucan
		kge_ix.append(compute_kge(values_v, values_vi))
		kge_x.append(compute_kge(values_v, values_vii))

	return iy, ix, kge_i, kge_ii, kge_iii, kge_iv, kge_v, kge_vi, kge_vii, kge_viii, kge_ix, kge_x
	

def import_smn_ii():
	
	iy, ix = [], []
	kge_i, kge_ii, kge_iii, kge_iv, kge_v, kge_vi, kge_vii, kge_viii, kge_ix, kge_x = [], [], [], [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 86):
		yy=smn_ii[i][1]
		xx=smn_ii[i][2]
		iy.append(smn_ii[i][1])
		ix.append(smn_ii[i][2])
		
		print('Reading weather station:', i, smn_ii[i][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		values_i = values_i*86400
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.month').mean('time')
		values_ii = d_ii.values
		values_ii = values_ii*86400
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.month').mean('time')
		values_iii = d_iii.values
		values_iii = values_iii*86400
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.month').mean('time')
		values_iv = d_iv.values
		values_iv = values_iv*86400
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.month').mean('time')
		values_v = d_v.values
		values_v = values_v*86400
			
		# Reading smn 
		d_vi = xr.open_dataset('{0}/FPS_SESA/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.month').mean('time')
		values_vi = d_vi.values
		values_vi = values_vi
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/FPS_SESA/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.groupby('time.month').mean('time')
		values_vii = d_vii.values

		# calculate kge regcm usp
		kge_i.append(compute_kge(values_i, values_vi))	
		kge_ii.append(compute_kge(values_i, values_vii))

		# calculate kge regcm ictp 1
		kge_iii.append(compute_kge(values_ii, values_vi))	
		kge_iv.append(compute_kge(values_ii, values_vii))

		# calculate kge regcm ictp 2
		kge_v.append(compute_kge(values_iii, values_vi))
		kge_vi.append(compute_kge(values_iii, values_vii))

		# calculate kge wrf ncar
		kge_vii.append(compute_kge(values_iv, values_vi))	
		kge_viii.append(compute_kge(values_iv, values_vii))

		# calculate kge wrf ucan
		kge_ix.append(compute_kge(values_v, values_vi))
		kge_x.append(compute_kge(values_v, values_vii))

	return iy, ix, kge_i, kge_ii, kge_iii, kge_iv, kge_v, kge_vi, kge_vii, kge_viii, kge_ix, kge_x
	
	
var = 'pr'

# Import dataset
lat_x, lon_x, kge_i_x, kge_ii_x, kge_iii_x, kge_iv_x, kge_v_x, kge_vi_x, kge_vii_x, kge_viii_x, kge_ix_x, kge_x_x = import_inmet()			
lat_y, lon_y, kge_i_y, kge_ii_y, kge_iii_y, kge_iv_y, kge_v_y, kge_vi_y, kge_vii_y, kge_viii_y, kge_ix_y, kge_x_y = import_smn_i()			
lat_z, lon_z, kge_i_z, kge_ii_z, kge_iii_z, kge_iv_z, kge_v_z, kge_vi_z, kge_vii_z, kge_viii_z, kge_ix_z, kge_x_z = import_smn_ii()			

lat_yy = lat_x + lat_y + lat_z
lon_xx = lon_x + lon_y + lon_z

reg_usp_inmet_smn = kge_i_x + kge_i_y + kge_i_z
reg_usp_reanalise = kge_ii_x + kge_ii_y + kge_ii_z

reg_ictp_i_inmet_smn = kge_iii_x + kge_iii_y + kge_iii_z
reg_ictp_i_reanalise = kge_iv_x + kge_iv_y + kge_iv_z

reg_ictp_ii_inmet_smn = kge_v_x + kge_v_y + kge_v_z
reg_ictp_ii_reanalise = kge_vi_x + kge_vi_y + kge_vi_z

wrf_ncar_inmet_smn = kge_vii_x + kge_vii_y + kge_vii_z
wrf_ncar_reanalise = kge_viii_x + kge_viii_y + kge_viii_z

wrf_ucan_inmet_smn = kge_ix_x + kge_ix_y + kge_ix_z
wrf_ucan_reanalise = kge_x_x + kge_x_y + kge_x_z

# Plot figure   
fig = plt.figure(figsize=(11, 4.5))

color = 'PRGn'
v_min = -0.9
v_max = 0.9
legend = 'KGE of precipitation'
font_size = 7
	
ax = fig.add_subplot(2, 5, 1)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(a) Reg4 vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.25, 0.010, 0.50]), extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(2, 5, 2)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_i_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(b) Reg5-Holt vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 3)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_ii_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(c) Reg5-UW vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 4)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(d) WRF-NCAR vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 5)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_inmet_smn, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(e) WRF-UCAN vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 6)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_usp_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(f) Reg4 vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')
plt.ylabel(u'Latitude', labelpad=20, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 7)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_i_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(g) Reg5-Holt vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 8)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, reg_ictp_ii_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(h) Reg5-UW vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 9)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ncar_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(i) WRF-NCAR vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(2, 5, 10)
my_map = basemap()
pltfig = my_map.scatter(lon_xx, lat_yy, 5, wrf_ucan_reanalise, cmap=color, marker='o', vmin=v_min, vmax=v_max)
plt.title('(j) WRF-UCAN vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel(u'Longitude', labelpad=15, fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '{0}/FPS_SESA/figs/paper_cp'.format(path)
name_out = 'pyplt_maps_kge_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
