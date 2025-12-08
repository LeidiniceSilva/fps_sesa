# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot maps of climatology"

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
from matplotlib.patches import Polygon
from matplotlib.colors import BoundaryNorm
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 'uv10' # t2m or uv10
font_size = 8
path = '/home/mda_silv/users/FPS_SESA'


def import_inmet():

	iy, ix = [], []
	mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii = [], [], [], [], [], [], [], []

	# Select lat and lon 
	for i in range(1, 99):
		yy=inmet[i][2]
		xx=inmet[i][3]
		iy.append(inmet[i][2])
		ix.append(inmet[i][3])
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])		
		if var == 't2m':
			# reading regcm usp 
			d_0 = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_0 = d_0.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_0 = d_0.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_0 = d_0.groupby('time.year').mean('time')
			d_0 = np.nanmean(d_0.values)
			mean_.append(d_0-273.15)
				
			# reading regcm ictp pbl 1 3 km 
			d_i = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_20180601-20211231.nc')
			d_i = d_i.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.groupby('time.year').mean('time')
			d_i = np.nanmean(d_i.values)
			mean_i.append(d_i)
					
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.groupby('time.year').mean('time')
			d_ii = np.nanmean(d_ii.values)
			mean_ii.append(d_ii-273.15)
			
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.groupby('time.year').mean('time')
			d_iii = np.nanmean(d_iii.values)
			mean_iii.append(d_iii-273.15)
							
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'tas_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.groupby('time.year').mean('time')
			d_iv = np.nanmean(d_iv.values)
			mean_iv.append(d_iv-273.15)
				
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.groupby('time.year').mean('time')
			d_v = np.nanmean(d_v.values)
			mean_v.append(d_v-273.15)
			
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/tmp/'.format(path) + 'tmp_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.tmp.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.year').mean('time')
			d_vi = np.nanmean(d_vi.values)
			mean_vi.append(d_vi)
			
			# reading era5 
			d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.t2m.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.groupby('time.year').mean('time')
			d_vii = np.nanmean(d_vii.values)
			mean_vii.append(d_vii-273.15)
			
		else:
			# reading regcm usp 
			d_0 = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_0 = d_0.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_0 = d_0.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_0 = d_0.groupby('time.year').mean('time')
			d_0 = np.nanmean(d_0.values)
			mean_.append(d_0)
				
			# reading regcm ictp pbl 1 3 km 
			d_i = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_20180601-20211231.nc')
			d_i = d_i.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.groupby('time.year').mean('time')
			d_i = np.nanmean(d_i.values)
			mean_i.append(d_i)
					
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.groupby('time.year').mean('time')
			d_ii = np.nanmean(d_ii.values)
			mean_ii.append(d_ii)
			
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.groupby('time.year').mean('time')
			d_iii = np.nanmean(d_iii.values)
			mean_iii.append(d_iii)
							
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'sfcWind_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.groupby('time.year').mean('time')
			d_iv = np.nanmean(d_iv.values)
			mean_iv.append(d_iv)
				
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.groupby('time.year').mean('time')
			d_v = np.nanmean(d_v.values)
			mean_v.append(d_v)
			
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/uv/'.format(path) + 'uv_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.uv.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.year').mean('time')
			d_vi = np.nanmean(d_vi.values)
			mean_vi.append(d_vi)
			
			# reading era5 
			d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.u10.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.groupby('time.year').mean('time')
			d_vii = np.nanmean(d_vii.values)
			mean_vii.append(d_vii)
	
	return iy, ix, mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def configure_subplot(ax):

	lon_bounds = [-62, -46]
	lat_bounds = [-36, -16]

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
lat_yy, lon_xx, reg_usp, reg_ictp, reg_ictp_i, reg_ictp_ii, wrf_ncar, wrf_ucan, inmet_smn, era5 = import_inmet()			

# Plot figure   
fig, axes = plt.subplots(3,3, figsize=(7, 8), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3), (ax4, ax5, ax6), (ax7, ax8, ax9) = axes
fig.delaxes(ax9)

if var == 't2m':
	vmin = 4
	vmax = 28
	n_classes = 16
	bins = np.linspace(vmin, vmax, n_classes + 1)
	cmap = cm.get_cmap('jet', n_classes)
	norm = BoundaryNorm(bins, cmap.N)
	legend = 'Temperature (°C)'
else:
	vmin = 0
	vmax = 6
	n_classes = 16
	bins = np.linspace(vmin, vmax, n_classes + 1)
	cmap = cm.get_cmap('viridis_r', n_classes)
	norm = BoundaryNorm(bins, cmap.N)
	legend = 'Wind 10m (m s⁻¹)'   
	    
st1 = ax1.scatter(lon_xx, lat_yy, 20, inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax1.set_title('(a) INMET', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

st2 = ax2.scatter(lon_xx, lat_yy, 20, era5, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax2.set_title('(b) ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

st3 = ax3.scatter(lon_xx, lat_yy, 20, reg_usp, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax3.set_title('(c) Reg4', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

st4 = ax4.scatter(lon_xx, lat_yy, 20, reg_ictp, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax4.set_title('(d) Reg5-Holt3', loc='left', fontsize=font_size, fontweight='bold')
ax4.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

st5 = ax5.scatter(lon_xx, lat_yy, 20, reg_ictp_i, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax5.set_title('(e) Reg5-Holt', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

st6 = ax6.scatter(lon_xx, lat_yy, 20, reg_ictp_ii, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax6.set_title('(f) Reg5-UW', loc='left', fontsize=font_size, fontweight='bold')
ax6.set_xlabel(u'Longitude',fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

st7 = ax7.scatter(lon_xx, lat_yy, 20, wrf_ncar, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax7.set_title('(g) WRF-NCAR', loc='left', fontsize=font_size, fontweight='bold')
ax7.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
ax7.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

st8 = ax8.scatter(lon_xx, lat_yy, 20, wrf_ucan, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax8.set_title('(h) WRF-UCAN', loc='left', fontsize=font_size, fontweight='bold')
ax8.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

cbar = plt.colorbar(st8, cax=fig.add_axes([0.91, 0.25, 0.019, 0.50]), extend='max')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_maps_clim_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()

