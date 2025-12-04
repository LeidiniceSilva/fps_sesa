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
from matplotlib.patches import Polygon
from matplotlib.colors import BoundaryNorm
from cartopy import config
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

var = 't2m' # t2m or uv10
font_size = 7
path = '/home/mda_silv/users/FPS_SESA'


def import_inmet():

	iy, ix = [], []
	bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x, bias_xi, bias_xii = [], [], [], [], [], [], [], [], [], [], [], []

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
			values_0 = np.nanmean(d_0.values)
			values_0 = values_0-273.15
				
			# reading regcm ictp pbl 1 3 km 
			d_i = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_20180601-20211231.nc')
			d_i = d_i.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.groupby('time.year').mean('time')
			values_i = np.nanmean(d_i.values)
		
			d_i = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.groupby('time.year').mean('time')
			values_i = np.nanmean(d_i.values)
			values_i = values_i-273.15
					
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.groupby('time.year').mean('time')
			values_ii = np.nanmean(d_ii.values)
			values_ii = values_ii-273.15
			
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = np.nanmean(d_iii.values)
			values_iii = values_iii-273.15
							
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'tas_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.groupby('time.year').mean('time')
			values_iv = np.nanmean(d_iv.values)
			values_iv = values_iv-273.15
				
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.groupby('time.year').mean('time')
			values_v = np.nanmean(d_v.values)
			values_v = values_v-273.15
				
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/tmp/'.format(path) + 'tmp_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.tmp.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.year').mean('time')
			values_vi = np.nanmean(d_vi.values)
			values_vi = values_vi
			
			# reading era5 
			d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.t2m.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.groupby('time.year').mean('time')
			values_vii = np.nanmean(d_vii.values)
			values_vii = values_vii-273.15
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
			values_ii = np.nanmean(d_ii.values)
			
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.groupby('time.year').mean('time')
			values_iii = np.nanmean(d_iii.values)
							
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'sfcWind_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.groupby('time.year').mean('time')
			values_iv = np.nanmean(d_iv.values)
				
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.groupby('time.year').mean('time')
			values_v = np.nanmean(d_v.values)
				
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/uv/'.format(path) + 'uv_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.uv.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.year').mean('time')
			values_vi = np.nanmean(d_vi.values)
			
			# reading era5 
			d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.u10.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.groupby('time.year').mean('time')
			values_vii = np.nanmean(d_vii.values)

		# calculate bias regcm usp
		bias_i.append(values_0 - values_vi)	
		bias_ii.append(values_0 - values_vii)

		# calculate bias regcm ictp
		bias_iii.append(values_i - values_vi)	
		bias_iv.append(values_i - values_vii)

		# calculate bias regcm ictp 1
		bias_v.append(values_ii - values_vi)	
		bias_vi.append(values_ii - values_vii)

		# calculate bias regcm ictp 2
		bias_vii.append(values_iii - values_vi)	
		bias_viii.append(values_iii - values_vii)

		# calculate bias wrf ncar
		bias_ix.append(values_iv - values_vi)	
		bias_x.append(values_iv - values_vii)

		# calculate bias wrf ucan
		bias_xi.append(values_v - values_vi)	
		bias_xii.append(values_v - values_vii)

	return iy, ix, bias_i, bias_ii, bias_iii, bias_iv, bias_v, bias_vi, bias_vii, bias_viii, bias_ix, bias_x, bias_xi, bias_xii


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
lat_yy, lon_xx, reg_usp_inmet_smn, reg_usp_reanalise, reg_ictp_inmet_smn, reg_ictp_reanalise, reg_ictp_i_inmet_smn, reg_ictp_i_reanalise, reg_ictp_ii_inmet_smn, reg_ictp_ii_reanalise, wrf_ncar_inmet_smn, wrf_ncar_reanalise, wrf_ucan_inmet_smn, wrf_ucan_reanalise = import_inmet()			

# Plot figure   
fig, axes = plt.subplots(3,4, figsize=(9, 7), subplot_kw={"projection": ccrs.PlateCarree()})
(ax1, ax2, ax3, ax4), (ax5, ax6, ax7, ax8), (ax9, ax10, ax11, ax12) = axes

vmin = -3
vmax = 3
n_classes = 16
bins = np.linspace(vmin, vmax, n_classes + 1)

if var == 't2m':
	cmap = cm.get_cmap('bwr', n_classes)
	norm = BoundaryNorm(bins, cmap.N)
	legend = 'Bias of temperature (°C)'
else:
	cmap = cm.get_cmap('PuOr', n_classes)
	norm = BoundaryNorm(bins, cmap.N)
	legend = 'Bias of wind 10m (m s⁻¹)'

st1 = ax1.scatter(lon_xx, lat_yy, 20, reg_usp_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax1.set_title('(a) Reg4 - INMET', loc='left', fontsize=font_size, fontweight='bold')
ax1.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax1)

st2 = ax2.scatter(lon_xx, lat_yy, 20, reg_usp_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax2.set_title('(b) Reg4 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax2)

st3 = ax3.scatter(lon_xx, lat_yy, 20, reg_ictp_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax3.set_title('(c) Reg5-holt3 - INMET', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax3)

st4 = ax4.scatter(lon_xx, lat_yy, 20, reg_ictp_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax4.set_title('(d) Reg5-holt3 - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax4)

st5 = ax5.scatter(lon_xx, lat_yy, 20, reg_ictp_i_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax5.set_title('(e) Reg5-holt - INMET', loc='left', fontsize=font_size, fontweight='bold')
ax5.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax5)

st6 = ax6.scatter(lon_xx, lat_yy, 20, reg_ictp_i_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax6.set_title('(f) Reg5-holt - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax6)

st7 = ax7.scatter(lon_xx, lat_yy, 20, reg_ictp_i_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax7.set_title('(g) Reg5-UW - INMET', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax7)

st8 = ax8.scatter(lon_xx, lat_yy, 20, reg_ictp_ii_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax8.set_title('(h) Reg5-UW - ERA5', loc='left', fontsize=font_size, fontweight='bold')
configure_subplot(ax8)

st9 = ax9.scatter(lon_xx, lat_yy, 20, wrf_ncar_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax9.set_title('(i) WRF-NCAR - INMET', loc='left', fontsize=font_size, fontweight='bold')
ax9.set_ylabel(u'Latitude', fontsize=font_size, fontweight='bold')
ax9.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax9)

st10 = ax10.scatter(lon_xx, lat_yy, 20, wrf_ncar_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax10.set_title('(j) WRF-NCAR - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax10.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax10)

st11 = ax11.scatter(lon_xx, lat_yy, 20, wrf_ucan_inmet_smn, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax11.set_title('(k) WRF-UCAN - INMET', loc='left', fontsize=font_size, fontweight='bold')
ax11.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax11)

st12 = ax12.scatter(lon_xx, lat_yy, 20, wrf_ucan_reanalise, cmap=cmap, norm=norm, marker='o', edgecolor='black', linewidth=0.5)
ax12.set_title('(l) WRF-UCAN - ERA5', loc='left', fontsize=font_size, fontweight='bold')
ax12.set_xlabel(u'Longitude', fontsize=font_size, fontweight='bold')
configure_subplot(ax12)

cbar = plt.colorbar(st1, cax=fig.add_axes([0.91, 0.25, 0.015, 0.50]), extend='both')
cbar.set_label('{0}'.format(legend), fontsize=font_size, fontweight='bold')
cbar.ax.tick_params(labelsize=font_size)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_maps_bias_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()






	


