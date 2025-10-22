# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sept 22, 2025"
__description__ = "This script plot maps of bias"

import os
import sys
import numpy as np
import pandas as pd
import xarray as xr
import cartopy.crs as ccrs
import cartopy.feature as cfeat
import matplotlib.pyplot as plt
import matplotlib.cm as cm

from matplotlib.path import Path
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'pr'
font_size = 7
path = '/home/mda_silv/users/FPS_SESA'


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
	
	
def import_inmet():

	iy, ix = [], []
	kge_i, kge_ii, kge_iii, kge_iv, kge_v, kge_vi, kge_vii, kge_viii, kge_ix, kge_x = [], [], [], [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 99):
		yy=inmet[i][2]
		xx=inmet[i][3]
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		values_i = values_i*86400
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.month').mean('time')
		values_ii = d_ii.values
		values_ii = values_ii*86400
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.month').mean('time')
		values_iii = d_iii.values
		values_iii = values_iii*86400
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.month').mean('time')
		values_iv = d_iv.values
		values_iv = values_iv*86400
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.month').mean('time')
		values_v = d_v.values
		values_v = values_v*86400
			
		# Reading inmet 
		d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.month').mean('time')
		values_vi = d_vi.values
		values_vi = values_vi*24
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
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
		d_i = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		values_i = values_i*86400
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.month').mean('time')
		values_ii = d_ii.values
		values_ii = values_ii*86400
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.month').mean('time')
		values_iii = d_iii.values
		values_iii = values_iii*86400
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.month').mean('time')
		values_iv = d_iv.values
		values_iv = values_iv*86400
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.month').mean('time')
		values_v = d_v.values
		values_v = values_v*86400
			
		# Reading smn 
		d_vi = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.month').mean('time')
		values_vi = d_vi.values
		values_vi = values_vi*24
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
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
		d_i = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.month').mean('time')
		values_i = d_i.values
		values_i = values_i*86400
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.month').mean('time')
		values_ii = d_ii.values
		values_ii = values_ii*86400
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.month').mean('time')
		values_iii = d_iii.values
		values_iii = values_iii*86400
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.month').mean('time')
		values_iv = d_iv.values
		values_iv = values_iv*86400
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.month').mean('time')
		values_v = d_v.values
		values_v = values_v*86400
			
		# Reading smn 
		d_vi = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.month').mean('time')
		values_vi = d_vi.values
		values_vi = values_vi
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
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
		

def configure_subplot(ax):

	lon_bounds = [-70, -46]
	lat_bounds = [-36, -18]

	states_provinces = cfeat.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lines', scale='50m', facecolor='none')

	ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]], crs=ccrs.PlateCarree())
	ax.set_xticks(np.arange(lon_bounds[0], lon_bounds[1], 4), crs=ccrs.PlateCarree())
	ax.set_yticks(np.arange(lat_bounds[0], lat_bounds[1], 4), crs=ccrs.PlateCarree())
	ax.xaxis.set_major_formatter(LongitudeFormatter())
	ax.yaxis.set_major_formatter(LatitudeFormatter())
	ax.grid(c='k', ls='--', alpha=0.5)  


	for label in ax.get_xticklabels() + ax.get_yticklabels():
		label.set_fontsize(font_size)
		
	ax.add_feature(states_provinces, edgecolor='0.05')
	ax.add_feature(cfeat.BORDERS, linewidth=0.75)
	ax.coastlines(linewidth=0.75)
	

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
fig, axes = plt.subplots(2,5, figsize=(13, 4), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3, ax4, ax5), (ax6, ax7, ax8, ax9, ax10) = axes

color = 'PRGn'
v_min = -0.9
v_max = 0.9
legend = 'KGE of precipitation'

st1 = ax1.scatter(lon_xx, lat_yy, 20, reg_usp_inmet_smn, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax1.set_title('(a) Reg4 vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

st2 = ax2.scatter(lon_xx, lat_yy, 20, reg_ictp_i_inmet_smn, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax2.set_title('(b) Reg5-Holt vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

st3 = ax3.scatter(lon_xx, lat_yy, 20, reg_ictp_ii_inmet_smn, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax3.set_title('(c) Reg5-UW vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

st4 = ax4.scatter(lon_xx, lat_yy, 20, wrf_ncar_inmet_smn, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax4.set_title('(d) WRF-NCAR vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

st5 = ax5.scatter(lon_xx, lat_yy, 20, wrf_ucan_inmet_smn, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax5.set_title('(e) WRF-UCAN vs. INMET+SMN', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

st6 = ax6.scatter(lon_xx, lat_yy, 20, reg_usp_reanalise, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax6.set_title('(f) Reg4 vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax6.set_xlabel(u'Longitude',fontsize=font_size, fontweight='bold')
ax6.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

st7 = ax7.scatter(lon_xx, lat_yy, 20, reg_ictp_i_reanalise, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax7.set_title('(g) Reg5-Holt vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax7.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

st8 = ax8.scatter(lon_xx, lat_yy, 20, reg_ictp_ii_reanalise, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax8.set_title('(h) Reg5-UW vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax8.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

st9 = ax9.scatter(lon_xx, lat_yy, 20, wrf_ncar_reanalise, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax9.set_title('(i) WRF-NCAR vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax9.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax9)

st10 = ax10.scatter(lon_xx, lat_yy, 20, wrf_ucan_reanalise, cmap=color, marker='o', edgecolor='black', linewidth=0.5, vmin=v_min, vmax=v_max)
ax10.set_title('(j) WRF-UCAN vs. ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax10.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax10)

cbar = plt.colorbar(st1, cax=fig.add_axes([0.25, -0.01, 0.5, 0.03]), orientation='horizontal', extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_maps_kge_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()




