# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot annual cycle"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii
from matplotlib.patches import Polygon

var = 'pr'
path = '/home/mda_silv/users/FPS_SESA'

skip_list_inmet_i = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] 
	
skip_list_inmet_ii = [2, 3, 4, 14, 19, 20, 21, 24, 25, 26, 27, 28, 32, 33, 34, 35, 38, 40, 41, 44, 45, 48, 52, 54, 55, 56, 59, 60, 62, 64, 68, 
70, 77, 79, 80, 82, 83, 92, 93, 96, 100, 106, 107, 111, 113, 120, 127, 130, 133, 135, 136, 140, 141, 144, 152, 154, 155, 160, 161, 163, 167, 168, 
173, 177, 180, 181, 182, 184, 186, 187, 188, 193, 197, 199, 204, 206, 207, 210, 212, 215, 216, 219, 220, 224, 225, 226, 229, 233, 237, 239, 240, 
241, 243, 248, 249, 251, 253, 254, 256, 261, 262, 264, 266, 269, 275, 276, 277, 280, 281, 282, 293, 295, 296, 298, 300, 303, 306, 308, 314, 315, 
316, 317, 319, 322, 325, 330, 331, 334, 337, 341, 344, 347, 348, 350, 353, 354, 357, 358, 360, 361, 362, 364, 370, 383, 384, 385, 389, 390, 392, 
393, 395, 396, 400, 401, 402, 404, 405, 408, 415, 416, 418, 423, 424, 427, 434, 440, 441, 443, 446, 448, 450, 451, 454, 455, 459, 465, 467, 471, 
474, 477, 481, 483, 488, 489, 492, 496, 504, 509, 513, 514, 516, 518, 519, 520, 523, 526, 528, 534, 538, 541, 544, 546, 552, 553, 557, 559]

skip_list_smn_ii = [39, 51, 55, 58, 64, 65, 66, 72, 75, 83, 86, 90, 91, 92]


def import_inmet():
	
	mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii  = [], [], [], [], [], [], [], []
	for i in range(1, 567):
		if i in skip_list_inmet_i:
			continue
		if i in skip_list_inmet_ii:
			continue
		if inmet[i][3] <= -48 and inmet[i][2] <= -16.5:
			yy=inmet[i][2]
			xx=inmet[i][3]
		
			print('Reading weather station:', i, inmet[i][0])		
			# reading regcm usp 
			d_ = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_ = d_.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_ = d_.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ = d_.groupby('time.month').mean('time')
			d_ = d_.values
			mean_.append(d_*86400)

			# reading regcm ictp pbl 1 3 km 
			d_i = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_20180601-20211231.nc')
			d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.groupby('time.month').mean('time')
			d_i = d_i.values
			mean_i.append(d_i)
						
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.groupby('time.month').mean('time')
			d_ii = d_ii.values
			mean_ii.append(d_ii*86400)
		
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.groupby('time.month').mean('time')
			d_iii = d_iii.values
			mean_iii.append(d_iii*86400)
						
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.groupby('time.month').mean('time')
			d_iv = d_iv.values
			mean_iv.append(d_iv*86400)
			
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.groupby('time.month').mean('time')
			d_v = d_v.values
			mean_v.append(d_v*86400)
		
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/pre/'.format(path) + 'pre_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.groupby('time.month').mean('time')
			d_vi = d_vi.values
			mean_vi.append(d_vi)
		
			# reading era5 
			d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.groupby('time.month').mean('time')
			d_vii = d_vii.values
			mean_vii.append(d_vii)
				
	return mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def import_smn_i():
	
	mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii  = [], [], [], [], [], [], [], []
	for i in range(1, 73):
		yy=smn_i[i][1]
		xx=smn_i[i][2]

		print('Reading weather station:', i, smn_i[i][0])	
		# reading regcm usp 
		d_ = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_ = d_.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ = d_.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ = d_.groupby('time.month').mean('time')
		d_ = d_.values
		mean_.append(d_*86400)

		# reading regcm ictp pbl 1 3 km 
		d_i = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_20180601-20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.month').mean('time')
		d_i = d_i.values
		mean_i.append(d_i)
				
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.month').mean('time')
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.month').mean('time')
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.month').mean('time')
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.month').mean('time')
		d_v = d_v.values
		mean_v.append(d_v*86400)
						
		# Reading smn 
		d_vi = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.month').mean('time')
		d_vi = d_vi.values
		mean_vi.append(d_vi*24)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.groupby('time.month').mean('time')
		d_vii = d_vii.values
		mean_vii.append(d_vii)
				
	return mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def import_smn_ii():
	
	mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii  = [], [], [], [], [], [], [], []
	for i in range(1, 110):
		if i in skip_list_smn_ii:
			continue
		yy=smn_ii[i][1]
		xx=smn_ii[i][2]
				
		print('Reading weather station:', i, smn_ii[i][0])	
		# reading regcm usp 
		d_ = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_ = d_.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ = d_.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ = d_.groupby('time.month').mean('time')
		d_ = d_.values
		mean_.append(d_*86400)

		# reading regcm ictp pbl 1 3 km 
		d_i = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v1_mon_20180601-20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.groupby('time.month').mean('time')
		d_i = d_i.values
		mean_i.append(d_i)
		
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.groupby('time.month').mean('time')
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.groupby('time.month').mean('time')
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.groupby('time.month').mean('time')
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.groupby('time.month').mean('time')
		d_v = d_v.values
		mean_v.append(d_v*86400)
						
		# Reading smn 
		d_vi = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.groupby('time.month').mean('time')
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')	
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.groupby('time.month').mean('time')
		d_vii = d_vii.values
		mean_vii.append(d_vii)
				
	return mean_, mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii
	

def portrait_diagram(ax, columns, data, line_label='Bias (Correlation)'):

	n_rows = 1  
	n_cols = len(columns)
	cell_width = 1
	cell_height = 1

	for j in range(n_cols):
		x = j * cell_width
		y = 0  
		era5_text, inmet_text = data[j]

		# ERA5
		era5_triangle = Polygon([[x, y+cell_height], [x, y], [x+cell_width, y+cell_height]], color='lightgray')
		ax.add_patch(era5_triangle)
		ax.text(x + 0.3*cell_width, y + 0.7*cell_height, era5_text, ha='center', va='center', fontsize=8)
	
		# INMET
		inmet_triangle = Polygon([[x, y], [x+cell_width, y], [x+cell_width, y+cell_height]], color='white')
		ax.add_patch(inmet_triangle)
		ax.text(x + 0.7*cell_width, y + 0.3*cell_height, inmet_text, ha='center', va='center', fontsize=8)

		ax.plot([x, x+cell_width, x+cell_width, x, x], [y, y, y+cell_height, y+cell_height, y], color='black')

	ax.set_xlim(0, n_cols*cell_width)
	ax.set_ylim(0, n_rows*cell_height)
	ax.set_xticks([x + cell_width/2 for x in range(n_cols)])
	ax.set_xticklabels(columns, fontsize=8)
	ax.set_yticks([0.5])
	ax.set_yticklabels([line_label], rotation=90, va='center', fontsize=8, fontweight='bold')
	ax.invert_yaxis()
	
	
# Import dataset
clim_0_x, clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x, clim_vi_x, clim_vii_x = import_inmet()			
clim_0_y, clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y, clim_vi_y, clim_vii_y = import_smn_i()			
clim_0_z, clim_i_z, clim_ii_z, clim_iii_z, clim_iv_z, clim_v_z, clim_vi_z, clim_vii_z = import_smn_ii()			

reg_usp      = clim_0_x + clim_0_y + clim_0_z
reg_ictp     = clim_i_x + clim_i_y + clim_i_z
reg_ictp_i_  = clim_ii_x + clim_ii_y + clim_ii_z
reg_ictp_ii_ = clim_iii_x + clim_iii_y + clim_iii_z
wrf_ncar     = clim_iv_x + clim_iv_y + clim_iv_z
wrf_ucan     = clim_v_x + clim_v_y + clim_v_z
inmet_smn    = clim_vi_x + clim_vi_y + clim_vi_z
era5         = clim_vii_x + clim_vii_y + clim_vii_z

list_hc = [1, 2, 3, 2, 0, 1, 1, 0, 2, 2, 0, 3, 0, 2, 3, 0, 1, 2, 0, 3, 0, 4, 2, 4, 3, 1, 4, 2, 4, 2, 2, 2, 1, 2, 4, 2, 2, 3, 2, 4, 4, 4, 0, 2, 4, 3, 2, 0, 0, 0, 3, 2, 2, 2, 1, 2, 4, 1, 4, 3, 4, 3, 0, 2, 0, 3, 2, 3, 2, 4, 0, 1, 4, 2, 4, 4, 0, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 0, 3, 2, 0, 0, 0, 4, 2, 3, 2, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 4, 0, 0, 4, 0, 4, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 2, 4, 3, 1, 4, 1, 2, 1, 1, 1, 4, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 2, 2, 2, 2, 2]

print(len(reg_usp))
print(len(reg_ictp))
print(len(reg_ictp_i_))
print(len(reg_ictp_ii_))
print(len(wrf_ncar))
print(len(wrf_ucan))
print(len(inmet_smn))
print(len(era5))
print(len(list_hc))

count_i, count_ii, count_iii, count_iv, count_v = [], [], [], [], []
for count, idx in enumerate(list_hc):
	if idx == 0:
		count_i.append(count)
	if idx == 1:
		count_ii.append(count)
	if idx == 2:
		count_iii.append(count)
	if idx == 3:
		count_iv.append(count)
	if idx == 4:
		count_v.append(count)

reg_usp_i,     reg_usp_ii,     reg_usp_iii,     reg_usp_iv,     reg_usp_v     = [], [], [], [], []
reg_ictp_i,    reg_ictp_ii,    reg_ictp_iii,    reg_ictp_iv,    reg_ictp_v    = [], [], [], [], []
reg_ictp_i_i,  reg_ictp_i_ii,  reg_ictp_i_iii,  reg_ictp_i_iv,  reg_ictp_i_v  = [], [], [], [], []
reg_ictp_ii_i, reg_ictp_ii_ii, reg_ictp_ii_iii, reg_ictp_ii_iv, reg_ictp_ii_v = [], [], [], [], []
wrf_ncar_i,    wrf_ncar_ii,    wrf_ncar_iii,    wrf_ncar_iv,    wrf_ncar_v    = [], [], [], [], []
wrf_ucan_i,    wrf_ucan_ii,    wrf_ucan_iii,    wrf_ucan_iv,    wrf_ucan_v    = [], [], [], [], []
inmet_smn_i,   inmet_smn_ii,   inmet_smn_iii,   inmet_smn_iv,   inmet_smn_v   = [], [], [], [], []
era5_i,        era5_ii,        era5_iii,        era5_iv,        era5_v        = [], [], [], [], []

for c_i in count_i:
	reg_usp_i.append(reg_usp[c_i])
	reg_ictp_i.append(reg_ictp[c_i])
	reg_ictp_i_i.append(reg_ictp_i_[c_i])
	reg_ictp_ii_i.append(reg_ictp_ii_[c_i])
	wrf_ncar_i.append(wrf_ncar[c_i])
	wrf_ucan_i.append(wrf_ucan[c_i])
	inmet_smn_i.append(inmet_smn[c_i])
	era5_i.append(era5[c_i])

for c_ii in count_ii:
	reg_usp_ii.append(reg_usp[c_ii])
	reg_ictp_ii.append(reg_ictp[c_ii])
	reg_ictp_i_ii.append(reg_ictp_i_[c_ii])
	reg_ictp_ii_ii.append(reg_ictp_ii_[c_ii])
	wrf_ncar_ii.append(wrf_ncar[c_ii])
	wrf_ucan_ii.append(wrf_ucan[c_ii])
	inmet_smn_ii.append(inmet_smn[c_ii])
	era5_ii.append(era5[c_ii])
	
for c_iii in count_iii:
	reg_usp_iii.append(reg_usp[c_iii])
	reg_ictp_iii.append(reg_ictp[c_iii])
	reg_ictp_i_iii.append(reg_ictp_i_[c_iii])
	reg_ictp_ii_iii.append(reg_ictp_ii_[c_iii])
	wrf_ncar_iii.append(wrf_ncar[c_iii])
	wrf_ucan_iii.append(wrf_ucan[c_iii])
	inmet_smn_iii.append(inmet_smn[c_iii])
	era5_iii.append(era5[c_iii])
	
for c_iv in count_iv:
	reg_usp_iv.append(reg_usp[c_iv])
	reg_ictp_iv.append(reg_ictp[c_iv])
	reg_ictp_i_iv.append(reg_ictp_i_[c_iv])
	reg_ictp_ii_iv.append(reg_ictp_ii_[c_iv])
	wrf_ncar_iv.append(wrf_ncar[c_iv])
	wrf_ucan_iv.append(wrf_ucan[c_iv])
	inmet_smn_iv.append(inmet_smn[c_iv])
	era5_iv.append(era5[c_iv])
	
for c_v in count_v:
	reg_usp_v.append(reg_usp[c_v])
	reg_ictp_v.append(reg_ictp[c_v])
	reg_ictp_i_v.append(reg_ictp_i_[c_v])
	reg_ictp_ii_v.append(reg_ictp_ii_[c_v])
	wrf_ncar_v.append(wrf_ncar[c_v])
	wrf_ucan_v.append(wrf_ucan[c_v])
	inmet_smn_v.append(inmet_smn[c_v])
	era5_v.append(era5[c_v])

# Group I
# Average
reg_usp_c_i     = np.nanmean(reg_usp_i, axis=0)
reg_ictp_c_i    = np.nanmean(reg_ictp_i, axis=0)
reg_ictp_i_c_i  = np.nanmean(reg_ictp_i_i, axis=0)
reg_ictp_ii_c_i = np.nanmean(reg_ictp_ii_i, axis=0)
wrf_ncar_c_i    = np.nanmean(wrf_ncar_i, axis=0)
wrf_ucan_c_i    = np.nanmean(wrf_ucan_i, axis=0)
inmet_smn_c_i   = np.nanmean(inmet_smn_i, axis=0)
era5_c_i        = np.nanmean(era5_i, axis=0)

# Bias
b_reg_usp_inmet_ci     = np.nanmean(reg_usp_c_i - inmet_smn_c_i)
b_reg_ictp_inmet_ci    = np.nanmean(reg_ictp_c_i - inmet_smn_c_i)
b_reg_ictp_i_inmet_ci  = np.nanmean(reg_ictp_i_c_i - inmet_smn_c_i)
b_reg_ictp_ii_inmet_ci = np.nanmean(reg_ictp_ii_c_i - inmet_smn_c_i)
b_wrf_ncar_inmet_ci    = np.nanmean(wrf_ncar_c_i - inmet_smn_c_i)
b_wrf_ucan_inmet_ci    = np.nanmean(wrf_ucan_c_i - inmet_smn_c_i)

b_reg_usp_era5_ci     = np.nanmean(reg_usp_c_i - era5_c_i)
b_reg_ictp_era5_ci    = np.nanmean(reg_ictp_c_i - era5_c_i)
b_reg_ictp_i_era5_ci  = np.nanmean(reg_ictp_i_c_i - era5_c_i)
b_reg_ictp_ii_era5_ci = np.nanmean(reg_ictp_ii_c_i - era5_c_i)
b_wrf_ncar_era5_ci    = np.nanmean(wrf_ncar_c_i - era5_c_i)
b_wrf_ucan_era5_ci    = np.nanmean(wrf_ucan_c_i - era5_c_i)

# Correlation
r_reg_usp_inmet_ci     = np.corrcoef(reg_usp_c_i, inmet_smn_c_i)[0, 1]
r_reg_ictp_inmet_ci    = np.corrcoef(reg_ictp_c_i, inmet_smn_c_i)[0, 1]
r_reg_ictp_i_inmet_ci  = np.corrcoef(reg_ictp_i_c_i, inmet_smn_c_i)[0, 1]
r_reg_ictp_ii_inmet_ci = np.corrcoef(reg_ictp_ii_c_i, inmet_smn_c_i)[0, 1]
r_wrf_ncar_inmet_ci    = np.corrcoef(wrf_ncar_c_i, inmet_smn_c_i)[0, 1]
r_wrf_ucan_inmet_ci    = np.corrcoef(wrf_ucan_c_i, inmet_smn_c_i)[0, 1]

r_reg_usp_era5_ci     = np.corrcoef(reg_usp_c_i, era5_c_i)[0, 1]
r_reg_ictp_era5_ci    = np.corrcoef(reg_ictp_c_i, era5_c_i)[0, 1]
r_reg_ictp_i_era5_ci  = np.corrcoef(reg_ictp_i_c_i, era5_c_i)[0, 1]
r_reg_ictp_ii_era5_ci = np.corrcoef(reg_ictp_ii_c_i, era5_c_i)[0, 1]
r_wrf_ncar_era5_ci    = np.corrcoef(wrf_ncar_c_i, era5_c_i)[0, 1]
r_wrf_ucan_era5_ci    = np.corrcoef(wrf_ucan_c_i, era5_c_i)[0, 1]

# Group II
# Average
reg_usp_c_ii     = np.nanmean(reg_usp_ii, axis=0)
reg_ictp_c_ii    = np.nanmean(reg_ictp_ii, axis=0)
reg_ictp_i_c_ii  = np.nanmean(reg_ictp_i_ii, axis=0)
reg_ictp_ii_c_ii = np.nanmean(reg_ictp_ii_ii, axis=0)
wrf_ncar_c_ii    = np.nanmean(wrf_ncar_ii, axis=0)
wrf_ucan_c_ii    = np.nanmean(wrf_ucan_ii, axis=0)
inmet_smn_c_ii   = np.nanmean(inmet_smn_ii, axis=0)
era5_c_ii        = np.nanmean(era5_ii, axis=0)

# Bias
b_reg_usp_inmet_cii     = np.nanmean(reg_usp_c_ii - inmet_smn_c_ii)
b_reg_ictp_inmet_cii    = np.nanmean(reg_ictp_c_ii - inmet_smn_c_ii)
b_reg_ictp_i_inmet_cii  = np.nanmean(reg_ictp_i_c_ii - inmet_smn_c_ii)
b_reg_ictp_ii_inmet_cii = np.nanmean(reg_ictp_ii_c_ii - inmet_smn_c_ii)
b_wrf_ncar_inmet_cii    = np.nanmean(wrf_ncar_c_ii - inmet_smn_c_ii)
b_wrf_ucan_inmet_cii    = np.nanmean(wrf_ucan_c_ii - inmet_smn_c_ii)

b_reg_usp_era5_cii     = np.nanmean(reg_usp_c_ii - era5_c_ii)
b_reg_ictp_era5_cii    = np.nanmean(reg_ictp_c_ii - era5_c_ii)
b_reg_ictp_i_era5_cii  = np.nanmean(reg_ictp_i_c_ii - era5_c_ii)
b_reg_ictp_ii_era5_cii = np.nanmean(reg_ictp_ii_c_ii - era5_c_ii)
b_wrf_ncar_era5_cii    = np.nanmean(wrf_ncar_c_ii - era5_c_ii)
b_wrf_ucan_era5_cii    = np.nanmean(wrf_ucan_c_ii - era5_c_ii)

# Correlation
r_reg_usp_inmet_cii     = np.corrcoef(reg_usp_c_ii, inmet_smn_c_ii)[0, 1]
r_reg_ictp_inmet_cii    = np.corrcoef(reg_ictp_c_ii, inmet_smn_c_ii)[0, 1]
r_reg_ictp_i_inmet_cii  = np.corrcoef(reg_ictp_i_c_ii, inmet_smn_c_ii)[0, 1]
r_reg_ictp_ii_inmet_cii = np.corrcoef(reg_ictp_ii_c_ii, inmet_smn_c_ii)[0, 1]
r_wrf_ncar_inmet_cii    = np.corrcoef(wrf_ncar_c_ii, inmet_smn_c_ii)[0, 1]
r_wrf_ucan_inmet_cii    = np.corrcoef(wrf_ucan_c_ii, inmet_smn_c_ii)[0, 1]

r_reg_usp_era5_cii     = np.corrcoef(reg_usp_c_ii, era5_c_ii)[0, 1]
r_reg_ictp_era5_cii    = np.corrcoef(reg_ictp_c_ii, era5_c_ii)[0, 1]
r_reg_ictp_i_era5_cii  = np.corrcoef(reg_ictp_i_c_ii, era5_c_ii)[0, 1]
r_reg_ictp_ii_era5_cii = np.corrcoef(reg_ictp_ii_c_ii, era5_c_ii)[0, 1]
r_wrf_ncar_era5_cii    = np.corrcoef(wrf_ncar_c_ii, era5_c_ii)[0, 1]
r_wrf_ucan_era5_cii    = np.corrcoef(wrf_ucan_c_ii, era5_c_ii)[0, 1]

# Group III
# Average
reg_usp_c_iii     = np.nanmean(reg_usp_iii, axis=0)
reg_ictp_c_iii    = np.nanmean(reg_ictp_iii, axis=0)
reg_ictp_i_c_iii  = np.nanmean(reg_ictp_i_iii, axis=0)
reg_ictp_ii_c_iii = np.nanmean(reg_ictp_ii_iii, axis=0)
wrf_ncar_c_iii    = np.nanmean(wrf_ncar_iii, axis=0)
wrf_ucan_c_iii    = np.nanmean(wrf_ucan_iii, axis=0)
inmet_smn_c_iii   = np.nanmean(inmet_smn_iii, axis=0)
era5_c_iii        = np.nanmean(era5_iii, axis=0)

# Bias
b_reg_usp_inmet_ciii     = np.nanmean(reg_usp_c_iii - inmet_smn_c_iii)
b_reg_ictp_inmet_ciii    = np.nanmean(reg_ictp_c_iii - inmet_smn_c_iii)
b_reg_ictp_i_inmet_ciii  = np.nanmean(reg_ictp_i_c_iii - inmet_smn_c_iii)
b_reg_ictp_ii_inmet_ciii = np.nanmean(reg_ictp_ii_c_iii - inmet_smn_c_iii)
b_wrf_ncar_inmet_ciii    = np.nanmean(wrf_ncar_c_iii - inmet_smn_c_iii)
b_wrf_ucan_inmet_ciii    = np.nanmean(wrf_ucan_c_iii - inmet_smn_c_iii)

b_reg_usp_era5_ciii     = np.nanmean(reg_usp_c_iii - era5_c_iii)
b_reg_ictp_era5_ciii    = np.nanmean(reg_ictp_c_iii - era5_c_iii)
b_reg_ictp_i_era5_ciii  = np.nanmean(reg_ictp_i_c_iii - era5_c_iii)
b_reg_ictp_ii_era5_ciii = np.nanmean(reg_ictp_ii_c_iii - era5_c_iii)
b_wrf_ncar_era5_ciii    = np.nanmean(wrf_ncar_c_iii - era5_c_iii)
b_wrf_ucan_era5_ciii    = np.nanmean(wrf_ucan_c_iii - era5_c_iii)

# Correlation
r_reg_usp_inmet_ciii     = np.corrcoef(reg_usp_c_iii, inmet_smn_c_iii)[0, 1]
r_reg_ictp_inmet_ciii    = np.corrcoef(reg_ictp_c_iii, inmet_smn_c_iii)[0, 1]
r_reg_ictp_i_inmet_ciii  = np.corrcoef(reg_ictp_i_c_iii, inmet_smn_c_iii)[0, 1]
r_reg_ictp_ii_inmet_ciii = np.corrcoef(reg_ictp_ii_c_iii, inmet_smn_c_iii)[0, 1]
r_wrf_ncar_inmet_ciii    = np.corrcoef(wrf_ncar_c_iii, inmet_smn_c_iii)[0, 1]
r_wrf_ucan_inmet_ciii    = np.corrcoef(wrf_ucan_c_iii, inmet_smn_c_iii)[0, 1]

r_reg_usp_era5_ciii     = np.corrcoef(reg_usp_c_iii, era5_c_iii)[0, 1]
r_reg_ictp_era5_ciii    = np.corrcoef(reg_ictp_c_iii, era5_c_iii)[0, 1]
r_reg_ictp_i_era5_ciii  = np.corrcoef(reg_ictp_i_c_iii, era5_c_iii)[0, 1]
r_reg_ictp_ii_era5_ciii = np.corrcoef(reg_ictp_ii_c_iii, era5_c_iii)[0, 1]
r_wrf_ncar_era5_ciii    = np.corrcoef(wrf_ncar_c_iii, era5_c_iii)[0, 1]
r_wrf_ucan_era5_ciii    = np.corrcoef(wrf_ucan_c_iii, era5_c_iii)[0, 1]

# Group IV
# Average
reg_usp_c_iv     = np.nanmean(reg_usp_iv, axis=0)
reg_ictp_c_iv    = np.nanmean(reg_ictp_iv, axis=0)
reg_ictp_i_c_iv  = np.nanmean(reg_ictp_i_iv, axis=0)
reg_ictp_ii_c_iv = np.nanmean(reg_ictp_ii_iv, axis=0)
wrf_ncar_c_iv    = np.nanmean(wrf_ncar_iv, axis=0)
wrf_ucan_c_iv    = np.nanmean(wrf_ucan_iv, axis=0)
inmet_smn_c_iv   = np.nanmean(inmet_smn_iv, axis=0)
era5_c_iv        = np.nanmean(era5_iv, axis=0)

# Bias
b_reg_usp_inmet_civ     = np.nanmean(reg_usp_c_iv - inmet_smn_c_iv)
b_reg_ictp_inmet_civ    = np.nanmean(reg_ictp_c_iv - inmet_smn_c_iv)
b_reg_ictp_i_inmet_civ  = np.nanmean(reg_ictp_i_c_iv - inmet_smn_c_iv)
b_reg_ictp_ii_inmet_civ = np.nanmean(reg_ictp_ii_c_iv - inmet_smn_c_iv)
b_wrf_ncar_inmet_civ    = np.nanmean(wrf_ncar_c_iv - inmet_smn_c_iv)
b_wrf_ucan_inmet_civ    = np.nanmean(wrf_ucan_c_iv - inmet_smn_c_iv)

b_reg_usp_era5_civ     = np.nanmean(reg_usp_c_iv - era5_c_iv)
b_reg_ictp_era5_civ    = np.nanmean(reg_ictp_c_iv - era5_c_iv)
b_reg_ictp_i_era5_civ  = np.nanmean(reg_ictp_i_c_iv - era5_c_iv)
b_reg_ictp_ii_era5_civ = np.nanmean(reg_ictp_ii_c_iv - era5_c_iv)
b_wrf_ncar_era5_civ    = np.nanmean(wrf_ncar_c_iv - era5_c_iv)
b_wrf_ucan_era5_civ    = np.nanmean(wrf_ucan_c_iv - era5_c_iv)

# Correlation
r_reg_usp_inmet_civ     = np.corrcoef(reg_usp_c_iv, inmet_smn_c_iv)[0, 1]
r_reg_ictp_inmet_civ    = np.corrcoef(reg_ictp_c_iv, inmet_smn_c_iv)[0, 1]
r_reg_ictp_i_inmet_civ  = np.corrcoef(reg_ictp_i_c_iv, inmet_smn_c_iv)[0, 1]
r_reg_ictp_ii_inmet_civ = np.corrcoef(reg_ictp_ii_c_iv, inmet_smn_c_iv)[0, 1]
r_wrf_ncar_inmet_civ    = np.corrcoef(wrf_ncar_c_iv, inmet_smn_c_iv)[0, 1]
r_wrf_ucan_inmet_civ    = np.corrcoef(wrf_ucan_c_iv, inmet_smn_c_iv)[0, 1]

r_reg_usp_era5_civ     = np.corrcoef(reg_usp_c_iv, era5_c_iv)[0, 1]
r_reg_ictp_era5_civ    = np.corrcoef(reg_ictp_c_iv, era5_c_iv)[0, 1]
r_reg_ictp_i_era5_civ  = np.corrcoef(reg_ictp_i_c_iv, era5_c_iv)[0, 1]
r_reg_ictp_ii_era5_civ = np.corrcoef(reg_ictp_ii_c_iv, era5_c_iv)[0, 1]
r_wrf_ncar_era5_civ    = np.corrcoef(wrf_ncar_c_iv, era5_c_iv)[0, 1]
r_wrf_ucan_era5_civ    = np.corrcoef(wrf_ucan_c_iv, era5_c_iv)[0, 1]

# Group V
# Average
reg_usp_c_v     = np.nanmean(reg_usp_v, axis=0)
reg_ictp_c_v    = np.nanmean(reg_ictp_v, axis=0)
reg_ictp_i_c_v  = np.nanmean(reg_ictp_i_v, axis=0)
reg_ictp_ii_c_v = np.nanmean(reg_ictp_ii_v, axis=0)
wrf_ncar_c_v    = np.nanmean(wrf_ncar_v, axis=0)
wrf_ucan_c_v    = np.nanmean(wrf_ucan_v, axis=0)
inmet_smn_c_v   = np.nanmean(inmet_smn_v, axis=0)
era5_c_v        = np.nanmean(era5_v, axis=0)

# Bias
b_reg_usp_inmet_cv     = np.nanmean(reg_usp_c_v - inmet_smn_c_v)
b_reg_ictp_inmet_cv    = np.nanmean(reg_ictp_c_v - inmet_smn_c_v)
b_reg_ictp_i_inmet_cv  = np.nanmean(reg_ictp_i_c_v - inmet_smn_c_v)
b_reg_ictp_ii_inmet_cv = np.nanmean(reg_ictp_ii_c_v - inmet_smn_c_v)
b_wrf_ncar_inmet_cv    = np.nanmean(wrf_ncar_c_v - inmet_smn_c_v)
b_wrf_ucan_inmet_cv    = np.nanmean(wrf_ucan_c_v - inmet_smn_c_v)

b_reg_usp_era5_cv     = np.nanmean(reg_usp_c_v - era5_c_v)
b_reg_ictp_era5_cv    = np.nanmean(reg_ictp_c_v - era5_c_v)
b_reg_ictp_i_era5_cv  = np.nanmean(reg_ictp_i_c_v - era5_c_v)
b_reg_ictp_ii_era5_cv = np.nanmean(reg_ictp_ii_c_v - era5_c_v)
b_wrf_ncar_era5_cv    = np.nanmean(wrf_ncar_c_v - era5_c_v)
b_wrf_ucan_era5_cv    = np.nanmean(wrf_ucan_c_v - era5_c_v)

# Correlation
r_reg_usp_inmet_cv     = np.corrcoef(reg_usp_c_v, inmet_smn_c_v)[0, 1]
r_reg_ictp_inmet_cv    = np.corrcoef(reg_ictp_c_v, inmet_smn_c_v)[0, 1]
r_reg_ictp_i_inmet_cv  = np.corrcoef(reg_ictp_i_c_v, inmet_smn_c_v)[0, 1]
r_reg_ictp_ii_inmet_cv = np.corrcoef(reg_ictp_ii_c_v, inmet_smn_c_v)[0, 1]
r_wrf_ncar_inmet_cv    = np.corrcoef(wrf_ncar_c_v, inmet_smn_c_v)[0, 1]
r_wrf_ucan_inmet_cv    = np.corrcoef(wrf_ucan_c_v, inmet_smn_c_v)[0, 1]

r_reg_usp_era5_cv     = np.corrcoef(reg_usp_c_v, era5_c_v)[0, 1]
r_reg_ictp_era5_cv    = np.corrcoef(reg_ictp_c_v, era5_c_v)[0, 1]
r_reg_ictp_i_era5_cv  = np.corrcoef(reg_ictp_i_c_v, era5_c_v)[0, 1]
r_reg_ictp_ii_era5_cv = np.corrcoef(reg_ictp_ii_c_v, era5_c_v)[0, 1]
r_wrf_ncar_era5_cv    = np.corrcoef(wrf_ncar_c_v, era5_c_v)[0, 1]
r_wrf_ucan_era5_cv    = np.corrcoef(wrf_ucan_c_v, era5_c_v)[0, 1]

# All groups 
# Average
reg_usp_c     = np.nanmean(reg_usp, axis=0)
reg_ictp_c    = np.nanmean(reg_ictp, axis=0)
reg_ictp_i_c  = np.nanmean(reg_ictp_i_, axis=0)
reg_ictp_ii_c = np.nanmean(reg_ictp_ii_, axis=0)
wrf_ncar_c    = np.nanmean(wrf_ncar, axis=0)
wrf_ucan_c    = np.nanmean(wrf_ucan, axis=0)
inmet_smn_c   = np.nanmean(inmet_smn, axis=0)
era5_c        = np.nanmean(era5, axis=0)

# Bias
b_reg_usp_inmet_c     = np.nanmean(reg_usp_c - inmet_smn_c)
b_reg_ictp_inmet_c    = np.nanmean(reg_ictp_c - inmet_smn_c)
b_reg_ictp_i_inmet_c  = np.nanmean(reg_ictp_i_c - inmet_smn_c)
b_reg_ictp_ii_inmet_c = np.nanmean(reg_ictp_ii_c - inmet_smn_c)
b_wrf_ncar_inmet_c    = np.nanmean(wrf_ncar_c - inmet_smn_c)
b_wrf_ucan_inmet_c    = np.nanmean(wrf_ucan_c - inmet_smn_c)

b_reg_usp_era5_c     = np.nanmean(reg_usp_c - era5_c)
b_reg_ictp_era5_c    = np.nanmean(reg_ictp_c - era5_c)
b_reg_ictp_i_era5_c  = np.nanmean(reg_ictp_i_c - era5_c)
b_reg_ictp_ii_era5_c = np.nanmean(reg_ictp_ii_c - era5_c)
b_wrf_ncar_era5_c    = np.nanmean(wrf_ncar_c - era5_c)
b_wrf_ucan_era5_c    = np.nanmean(wrf_ucan_c - era5_c)

# Correlation
r_reg_usp_inmet_c     = np.corrcoef(reg_usp_c, inmet_smn_c)[0, 1]
r_reg_ictp_inmet_c    = np.corrcoef(reg_ictp_c, inmet_smn_c)[0, 1]
r_reg_ictp_i_inmet_c  = np.corrcoef(reg_ictp_i_c, inmet_smn_c)[0, 1]
r_reg_ictp_ii_inmet_c = np.corrcoef(reg_ictp_ii_c, inmet_smn_c)[0, 1]
r_wrf_ncar_inmet_c    = np.corrcoef(wrf_ncar_c, inmet_smn_c)[0, 1]
r_wrf_ucan_inmet_c    = np.corrcoef(wrf_ucan_c, inmet_smn_c)[0, 1]

r_reg_usp_era5_c     = np.corrcoef(reg_usp_c, era5_c)[0, 1]
r_reg_ictp_era5_c    = np.corrcoef(reg_ictp_c, era5_c)[0, 1]
r_reg_ictp_i_era5_c  = np.corrcoef(reg_ictp_i_c, era5_c)[0, 1]
r_reg_ictp_ii_era5_c = np.corrcoef(reg_ictp_ii_c, era5_c)[0, 1]
r_wrf_ncar_era5_c    = np.corrcoef(wrf_ncar_c, era5_c)[0, 1]
r_wrf_ucan_era5_c    = np.corrcoef(wrf_ucan_c, era5_c)[0, 1]

data_c_i = [[f"{b_e_ci:.2f}\n({r_e_ci:.2f})", f"{b_i_ci:.2f}\n({r_i_ci:.2f})"]
    for b_e_ci, b_i_ci, r_e_ci, r_i_ci in zip([b_reg_usp_era5_ci,  b_reg_ictp_era5_ci,  b_reg_ictp_i_era5_ci,  b_reg_ictp_ii_era5_ci,  b_wrf_ncar_era5_ci,  b_wrf_ucan_era5_ci],
        [b_reg_usp_inmet_ci, b_reg_ictp_inmet_ci, b_reg_ictp_i_inmet_ci, b_reg_ictp_ii_inmet_ci, b_wrf_ncar_inmet_ci, b_wrf_ucan_inmet_ci],
        [r_reg_usp_era5_ci,  r_reg_ictp_era5_ci,  r_reg_ictp_i_era5_ci,  r_reg_ictp_ii_era5_ci,  r_wrf_ncar_era5_ci,  r_wrf_ucan_era5_ci],
        [r_reg_usp_inmet_ci, r_reg_ictp_inmet_ci, r_reg_ictp_i_inmet_ci, r_reg_ictp_ii_inmet_ci, r_wrf_ncar_inmet_ci, r_wrf_ucan_inmet_ci])]

data_c_ii = [[f"{b_e_cii:.2f}\n({r_e_cii:.2f})", f"{b_i_cii:.2f}\n({r_i_cii:.2f})"]
    for b_e_cii, b_i_cii, r_e_cii, r_i_cii in zip([b_reg_usp_era5_cii,  b_reg_ictp_era5_cii,  b_reg_ictp_i_era5_cii,  b_reg_ictp_ii_era5_cii,  b_wrf_ncar_era5_cii,  b_wrf_ucan_era5_cii],
        [b_reg_usp_inmet_cii, b_reg_ictp_inmet_cii, b_reg_ictp_i_inmet_cii, b_reg_ictp_ii_inmet_cii, b_wrf_ncar_inmet_cii, b_wrf_ucan_inmet_cii],
        [r_reg_usp_era5_cii,  r_reg_ictp_era5_cii,  r_reg_ictp_i_era5_cii,  r_reg_ictp_ii_era5_cii,  r_wrf_ncar_era5_cii,  r_wrf_ucan_era5_cii],
        [r_reg_usp_inmet_cii, r_reg_ictp_inmet_cii, r_reg_ictp_i_inmet_cii, r_reg_ictp_ii_inmet_cii, r_wrf_ncar_inmet_cii, r_wrf_ucan_inmet_cii])]


data_c_iii = [[f"{b_e_ciii:.2f}\n({r_e_ciii:.2f})", f"{b_i_ciii:.2f}\n({r_i_ciii:.2f})"]
    for b_e_ciii, b_i_ciii, r_e_ciii, r_i_ciii in zip(	[b_reg_usp_era5_ciii,  b_reg_ictp_era5_ciii,  b_reg_ictp_i_era5_ciii,  b_reg_ictp_ii_era5_ciii,  b_wrf_ncar_era5_ciii,  b_wrf_ucan_era5_ciii],
        [b_reg_usp_inmet_ciii, b_reg_ictp_inmet_ciii, b_reg_ictp_i_inmet_ciii, b_reg_ictp_ii_inmet_ciii, b_wrf_ncar_inmet_ciii, b_wrf_ucan_inmet_ciii],
        [r_reg_usp_era5_ciii,  r_reg_ictp_era5_ciii,  r_reg_ictp_i_era5_ciii,  r_reg_ictp_ii_era5_ciii,  r_wrf_ncar_era5_ciii,  r_wrf_ucan_era5_ciii],
        [r_reg_usp_inmet_ciii, r_reg_ictp_inmet_ciii, r_reg_ictp_i_inmet_ciii, r_reg_ictp_ii_inmet_ciii, r_wrf_ncar_inmet_ciii, r_wrf_ucan_inmet_ciii])]


data_c_iv = [[f"{b_e_civ:.2f}\n({r_e_civ:.2f})", f"{b_i_civ:.2f}\n({r_i_civ:.2f})"]
    for b_e_civ, b_i_civ, r_e_civ, r_i_civ in zip([b_reg_usp_era5_civ,  b_reg_ictp_era5_civ,  b_reg_ictp_i_era5_civ,  b_reg_ictp_ii_era5_civ,  b_wrf_ncar_era5_civ,  b_wrf_ucan_era5_civ],
        [b_reg_usp_inmet_civ, b_reg_ictp_inmet_civ, b_reg_ictp_i_inmet_civ, b_reg_ictp_ii_inmet_civ, b_wrf_ncar_inmet_civ, b_wrf_ucan_inmet_civ],
        [r_reg_usp_era5_civ,  r_reg_ictp_era5_civ,  r_reg_ictp_i_era5_civ,  r_reg_ictp_ii_era5_civ,  r_wrf_ncar_era5_civ,  r_wrf_ucan_era5_civ],
        [r_reg_usp_inmet_civ, r_reg_ictp_inmet_civ, r_reg_ictp_i_inmet_civ, r_reg_ictp_ii_inmet_civ, r_wrf_ncar_inmet_civ, r_wrf_ucan_inmet_civ])]

data_c_v = [[f"{b_e_cv:.2f}\n({r_e_cv:.2f})", f"{b_i_cv:.2f}\n({r_i_cv:.2f})"]
    for b_e_cv, b_i_cv, r_e_cv, r_i_cv in zip([b_reg_usp_era5_cv,  b_reg_ictp_era5_cv,  b_reg_ictp_i_era5_cv,  b_reg_ictp_ii_era5_cv,  b_wrf_ncar_era5_cv,  b_wrf_ucan_era5_cv],
        [b_reg_usp_inmet_cv, b_reg_ictp_inmet_cv, b_reg_ictp_i_inmet_cv, b_reg_ictp_ii_inmet_cv, b_wrf_ncar_inmet_cv, b_wrf_ucan_inmet_cv],
        [r_reg_usp_era5_cv,  r_reg_ictp_era5_cv,  r_reg_ictp_i_era5_cv,  r_reg_ictp_ii_era5_cv,  r_wrf_ncar_era5_cv,  r_wrf_ucan_era5_cv],
        [r_reg_usp_inmet_cv, r_reg_ictp_inmet_cv, r_reg_ictp_i_inmet_cv, r_reg_ictp_ii_inmet_cv, r_wrf_ncar_inmet_cv, r_wrf_ucan_inmet_cv])]

data_c = [[f"{b_e_c:.2f}\n({r_e_c:.2f})", f"{b_i_c:.2f}\n({r_i_c:.2f})"]
    for b_e_c, b_i_c, r_e_c, r_i_c in zip([b_reg_usp_era5_c,  b_reg_ictp_era5_c,  b_reg_ictp_i_era5_c,  b_reg_ictp_ii_era5_c,  b_wrf_ncar_era5_c,  b_wrf_ucan_era5_c],
        [b_reg_usp_inmet_c, b_reg_ictp_inmet_c, b_reg_ictp_i_inmet_c, b_reg_ictp_ii_inmet_c, b_wrf_ncar_inmet_c, b_wrf_ucan_inmet_c],
        [r_reg_usp_era5_c,  r_reg_ictp_era5_c,  r_reg_ictp_i_era5_c,  r_reg_ictp_ii_era5_c,  r_wrf_ncar_era5_c,  r_wrf_ucan_era5_c],
        [r_reg_usp_inmet_c, r_reg_ictp_inmet_c, r_reg_ictp_i_inmet_c, r_reg_ictp_ii_inmet_c, r_wrf_ncar_inmet_c, r_wrf_ucan_inmet_c])]

# Plot figure
fig = plt.figure(figsize=(10, 12))
time = np.arange(0.5, 12 + 0.5)
font_size = 8

legend = 'Precipitation (mm d$^-$$^1$)'
vmin = 0
vmax = 10
vmax_ = 11
int_ = 1
	
columns = ['Reg4', 'Reg5-Holt3', 'Reg5-Holt', 'Reg5-UW', 'WRF-NCAR', 'WRF-UCAN']

ax = fig.add_subplot(6, 2, 1)
plt.plot(time, inmet_smn_c_i,   linewidth=1.5, color='black', label='INMET+SMN')
plt.plot(time, era5_c_i,        linewidth=1.5, color='red',   label='ERA5')
plt.plot(time, reg_usp_c_i,     linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='blue',    label='Reg4')
plt.plot(time, reg_ictp_c_i,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='magenta', label='Reg5-Holt3')
plt.plot(time, reg_ictp_i_c_i,  linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='gray',    label='Reg5-Holt')
plt.plot(time, reg_ictp_ii_c_i, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='brown',   label='Reg5-UW')
plt.plot(time, wrf_ncar_c_i,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='green',   label='WRF-NCAR')
plt.plot(time, wrf_ucan_c_i,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='orange',  label='WRF-UCAN')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(6, 2, 2)
plt.plot(time, inmet_smn_c_ii,   linewidth=1.5, color='black', label='INMET+SMN')
plt.plot(time, era5_c_ii,        linewidth=1.5, color='red',   label='ERA5')
plt.plot(time, reg_usp_c_ii,     linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='blue',    label='Reg4')
plt.plot(time, reg_ictp_c_ii,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='magenta', label='Reg5-Holt3')
plt.plot(time, reg_ictp_i_c_ii,  linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='gray',    label='Reg5-Holt')
plt.plot(time, reg_ictp_ii_c_ii, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='brown',   label='Reg5-UW')
plt.plot(time, wrf_ncar_c_ii,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='green',   label='WRF-NCAR')
plt.plot(time, wrf_ucan_c_ii,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='orange',  label='WRF-UCAN')
plt.title('(b) Cluster II', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(6, 2, 3)
portrait_diagram(ax, columns, data_c_i)
ax.xaxis.set_visible(False)

ax = fig.add_subplot(6, 2, 4)
portrait_diagram(ax, columns, data_c_ii)
ax.xaxis.set_visible(False)

ax = fig.add_subplot(6, 2, 5)
plt.plot(time, inmet_smn_c_iii,   linewidth=1.5, color='black', label='INMET+SMN')
plt.plot(time, era5_c_iii,        linewidth=1.5, color='red',   label='ERA5')
plt.plot(time, reg_usp_c_iii,     linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='blue',    label='Reg4')
plt.plot(time, reg_ictp_c_iii,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='magenta', label='Reg5-Holt3')
plt.plot(time, reg_ictp_i_c_iii,  linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='gray',    label='Reg5-Holt')
plt.plot(time, reg_ictp_ii_c_iii, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='brown',   label='Reg5-UW')
plt.plot(time, wrf_ncar_c_iii,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='green',   label='WRF-NCAR')
plt.plot(time, wrf_ucan_c_iii,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='orange',  label='WRF-UCAN')
plt.title('(c) Cluster III', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(6, 2, 6)
plt.plot(time, inmet_smn_c_iv,   linewidth=1.5, color='black', label='INMET+SMN')
plt.plot(time, era5_c_iv,        linewidth=1.5, color='red',   label='ERA5')
plt.plot(time, reg_usp_c_iv,     linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='blue',    label='Reg4')
plt.plot(time, reg_ictp_c_iv,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='magenta', label='Reg5-Holt3')
plt.plot(time, reg_ictp_i_c_iv,  linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='gray',    label='Reg5-Holt')
plt.plot(time, reg_ictp_ii_c_iv, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='brown',   label='Reg5-UW')
plt.plot(time, wrf_ncar_c_iv,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='green',   label='WRF-NCAR')
plt.plot(time, wrf_ucan_c_iv,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='orange',  label='WRF-UCAN')
plt.title('(d) Cluster IV', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(6, 2, 7)
portrait_diagram(ax, columns, data_c_iii)
ax.xaxis.set_visible(False)

ax = fig.add_subplot(6, 2, 8)
portrait_diagram(ax, columns, data_c_iv)
ax.xaxis.set_visible(False)

ax = fig.add_subplot(6, 2, 9)
plt.plot(time, inmet_smn_c_v,   linewidth=1.5, color='black', label='INMET+SMN')
plt.plot(time, era5_c_v,        linewidth=1.5, color='red',   label='ERA5')
plt.plot(time, reg_usp_c_v,     linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='blue',    label='Reg4')
plt.plot(time, reg_ictp_c_v,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='magenta', label='Reg5-Holt3')
plt.plot(time, reg_ictp_i_c_v,  linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='gray',    label='Reg5-Holt')
plt.plot(time, reg_ictp_ii_c_v, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='brown',   label='Reg5-UW')
plt.plot(time, wrf_ncar_c_v,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='green',   label='WRF-NCAR')
plt.plot(time, wrf_ucan_c_v,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='orange',  label='WRF-UCAN')
plt.title('(e) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.grid(linestyle='--')

ax = fig.add_subplot(6, 2, 10)
plt.plot(time, inmet_smn_c,   linewidth=1.5, color='black', label='INMET+SMN')
plt.plot(time, era5_c,        linewidth=1.5, color='red',   label='ERA5')
plt.plot(time, reg_usp_c,     linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='blue',    label='Reg4')
plt.plot(time, reg_ictp_c,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='magenta', label='Reg5-Holt3')
plt.plot(time, reg_ictp_i_c,  linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='gray',    label='Reg5-Holt')
plt.plot(time, reg_ictp_ii_c, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='brown',   label='Reg5-UW')
plt.plot(time, wrf_ncar_c,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='green',   label='WRF-NCAR')
plt.plot(time, wrf_ucan_c,    linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='orange',  label='WRF-UCAN')
plt.title('(f) All clusters', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'), fontsize=font_size)
plt.legend(ncol=8, fontsize=font_size, loc=(-1.30, -1.55))
plt.grid(linestyle='--')

ax = fig.add_subplot(6, 2, 11)
portrait_diagram(ax, columns, data_c_v)

ax = fig.add_subplot(6, 2, 12)
portrait_diagram(ax, columns, data_c)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_annual_cycle_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
