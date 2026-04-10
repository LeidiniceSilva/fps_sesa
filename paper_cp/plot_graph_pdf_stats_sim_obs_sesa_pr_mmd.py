# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot pdf"

import os
import argparse
import numpy as np
import xarray as xr
import seaborn as sns
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

parser = argparse.ArgumentParser(description='Process variable')
parser.add_argument('--var', required=True, choices=['pr'], help='Variable name')
args = parser.parse_args()
var = args.var

dict_var = {'pr': ['pre', 'tp']}

vmin, vmax, step_ = 0.1, 500, 0.5

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

	mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii = [], [], [], [], [], [], [], []

	for i in range(1, 567):
		if i in skip_list_inmet_i:
			continue
		if i in skip_list_inmet_ii:
			continue 
		if inmet[i][3] <= -48 and inmet[i][2] <= -16.5:
			station_code = inmet[i][0]
			station_name = inmet[i][1]
			print(i, station_code, station_name)
		
			# Reading inmet 
			d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/{1}/'.format(path, dict_var[var][0]) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(dict_var[var][0], station_code))
			d_i = d_i[dict_var[var][0]].sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.resample(time='1D').mean()
			mean_i.append(d_i.values*24)

			# Reading era5 
			d_ii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/{0}/'.format(dict_var[var][1]) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(dict_var[var][1], station_code, station_name))
			d_ii = d_ii[dict_var[var][1]].sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.resample(time='1D').mean()
			mean_ii.append(d_ii.values*24)

			# Reading regcm usp
			d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iii = d_iii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.resample(time='1D').mean()
			mean_iii.append(d_iii.values)
			
			# Reading regcm ictp 
			d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_iv = d_iv[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.resample(time='1D').mean()
			mean_iv.append(d_iv.values*24)
	
			# Reading regcm ictp pbl 1
			d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_v = d_v[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.resample(time='1D').mean()
			mean_v.append(d_v.values)

			# Reading regcm ictp pbl 2
			d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_vi = d_vi[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.resample(time='1D').mean()
			mean_vi.append(d_vi.values)

			# Reading wrf ncar
			d_vii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ncar/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_vii = d_vii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.resample(time='1D').mean()
			mean_vii.append(d_vii.values)
		
			# Reading wrf ucan
			d_viii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ucan/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
			d_viii = d_viii[var].sel(time=slice('2018-06-01','2021-05-31'))
			d_viii = d_viii.resample(time='1D').mean()
			mean_viii.append(d_viii.values)
		
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii


def import_smn_i():

	mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii = [], [], [], [], [], [], [], []

	for i in range(1, 73):
		station_code = f'SMN{i:03d}'
		station_name = smn_i[i][0]		
		print(i, station_code, station_name)
		
		# Reading smn 
		d_i = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/smn_i/smn_nc/'.format(path) + '{0}_{1}_H_2018-01-01_2021-12-31.nc'.format(dict_var[var][0], smn_i[i][0]))
		d_i = d_i[dict_var[var][0]].sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.resample(time='1D').mean()
		mean_i.append(d_i.values*24)

		# Reading era5 
		d_ii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/{0}/'.format(dict_var[var][1]) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(dict_var[var][1], station_code, station_name))
		d_ii = d_ii[dict_var[var][1]].sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.resample(time='1D').mean()
		mean_ii.append(d_ii.values*24)

		# Reading regcm usp
		d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_iii = d_iii[var].sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.resample(time='1D').mean()
		mean_iii.append(d_iii.values)
			
		# Reading regcm ictp 
		d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_iv = d_iv[var].sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.resample(time='1D').mean()
		mean_iv.append(d_iv.values*24)
	
		# Reading regcm ictp pbl 1
		d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_v = d_v[var].sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.resample(time='1D').mean()
		mean_v.append(d_v.values)

		# Reading regcm ictp pbl 2
		d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_vi = d_vi[var].sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.resample(time='1D').mean()
		mean_vi.append(d_vi.values)

		# Reading wrf ncar
		d_vii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ncar/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_vii = d_vii[var].sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.resample(time='1D').mean()
		mean_vii.append(d_vii.values)
		
		# Reading wrf ucan
		d_viii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ucan/{0}/'.format(var) + '{0}_{1}_{2}_H_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_viii = d_viii[var].sel(time=slice('2018-06-01','2021-05-31'))
		d_viii = d_viii.resample(time='1D').mean()
		mean_viii.append(d_viii.values)
			
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii


def import_smn_ii():
	
	mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii = [], [], [], [], [], [], [], []

	for i in range(1, 110):
		if i in skip_list_smn_ii:
			continue
		station_code = f'SMN{i:03d}'
		station_name = smn_ii[i][0]
		print(i, station_code, station_name)
					
		# Reading smn 
		d_i = xr.open_dataset('/home/mda_silv/users/FPS_SESA/database/obs/smn_ii/smn_nc/{1}/'.format(path, dict_var[var][0]) + '{0}_{1}_D_1979-01-01_2021-12-31.nc'.format(dict_var[var][0], smn_ii[i][0]))
		d_i = d_i[dict_var[var][0]].sel(time=slice('2018-06-01','2021-05-31'))
		mean_i.append(d_i.values)

		# Reading era5 
		d_ii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/{0}/'.format(dict_var[var][1]) + '{0}_{1}_{2}_D_2018-06-01-2021-05-31.nc'.format(dict_var[var][1], station_code, station_name))
		d_ii = d_ii[dict_var[var][1]].sel(time=slice('2018-06-01','2021-05-31'))
		mean_ii.append(d_ii.values)

		# Reading regcm usp
		d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/{0}/'.format(var) + '{0}_{1}_{2}_D_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_iii = d_iii[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_iii.append(d_iii.values/24)
			
		# Reading regcm ictp 
		d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/{0}/'.format(var) + '{0}_{1}_{2}_D_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_iv = d_iv[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_iv.append(d_iv.values)
	
		# Reading regcm ictp pbl 1
		d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/{0}/'.format(var) + '{0}_{1}_{2}_D_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_v = d_v[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_v.append(d_v.values/24)

		# Reading regcm ictp pbl 2
		d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/{0}/'.format(var) + '{0}_{1}_{2}_D_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_vi = d_vi[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_vi.append(d_vi.values/24)

		# Reading wrf ncar
		d_vii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ncar/{0}/'.format(var) + '{0}_{1}_{2}_D_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_vii = d_vii[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_vii.append(d_vii.values/24)
		
		# Reading wrf ucan
		d_viii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ucan/{0}/'.format(var) + '{0}_{1}_{2}_D_2018-06-01-2021-05-31.nc'.format(var, station_code, station_name))
		d_viii = d_viii[var].sel(time=slice('2018-06-01','2021-05-31'))
		mean_viii.append(d_viii.values/24)

	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii


def mask_like(reference, target):

	reference = np.asarray(reference)
	masked_target = np.asarray(target, dtype=float).copy()
	masked_target[np.isnan(reference)] = np.nan
	
	return masked_target
	
    
def compute_pdf(value, min_val=vmin, max_val=vmax, step=step_):

	valid = value[~np.isnan(value)]	
	valid = valid[(valid > 0) & (valid < 500)]
	
	perc = np.nanpercentile(valid, 99)
	max_ = np.nanmax(valid)
    
	bins = np.arange(min_val, max_val + step, step)
	wet = valid[valid >= min_val]
	hist, _ = np.histogram(wet, bins=bins)
	hist = hist.astype(float)
	hist[hist < 1] = np.nan
	pdf = hist / (np.nansum(hist) * step)
	bin_centers = (bins[:-1] + bins[1:] + step) / 2

	return perc, max_, bin_centers, pdf
	

# Import dataset
mean_i_x, mean_ii_x, mean_iii_x, mean_iv_x, mean_v_x, mean_vi_x, mean_vii_x, mean_viii_x = import_inmet()			
mean_i_y, mean_ii_y, mean_iii_y, mean_iv_y, mean_v_y, mean_vi_y, mean_vii_y, mean_viii_y = import_smn_i()			
mean_i_z, mean_ii_z, mean_iii_z, mean_iv_z, mean_v_z, mean_vi_z, mean_vii_z, mean_viii_z = import_smn_ii()			

inmet_smn    = mean_i_x    + mean_i_y    + mean_i_z
era5         = mean_ii_x   + mean_ii_y   + mean_ii_z
reg_usp      = mean_iii_x  + mean_iii_y  + mean_iii_z
reg_ictp     = mean_iv_x   + mean_iv_y   + mean_iv_z
reg_ictp_i_  = mean_v_x    + mean_v_y    + mean_v_z
reg_ictp_ii_ = mean_vi_x   + mean_vi_y   + mean_vi_z
wrf_ncar     = mean_vii_x  + mean_vii_y  + mean_vii_z
wrf_ucan     = mean_viii_x + mean_viii_y + mean_viii_z
	
list_hc = [1, 2, 3, 2, 0, 1, 1, 0, 2, 2, 0, 3, 0, 2, 3, 0, 1, 2, 0, 3, 0, 4, 2, 4, 3, 1, 4, 2, 4, 2, 2, 2, 1, 2, 4, 2, 2, 3, 2, 4, 4, 4, 0, 2, 4, 3, 2, 0, 0, 0, 3, 2, 2, 2, 1, 2, 4, 1, 4, 3, 4, 3, 0, 2, 0, 3, 2, 3, 2, 4, 0, 1, 4, 2, 4, 4, 0, 0, 2, 0, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 0, 3, 2, 0, 0, 0, 4, 2, 3, 2, 2, 2, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 2, 1, 1, 4, 0, 0, 4, 0, 4, 1, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 2, 1, 2, 4, 4, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 4, 0, 2, 4, 3, 1, 4, 1, 2, 1, 1, 1, 4, 1, 1, 2, 0, 0, 0, 0, 0, 0, 0, 2, 0, 2, 0, 0, 0, 0, 0, 0, 1, 1, 1, 4, 4, 4, 4, 2, 2, 4, 4, 2, 4, 2, 2, 2, 2, 2]
list_hc = list_hc[:len(inmet_smn)]
print(len(inmet_smn))

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

inmet_smn_i,   inmet_smn_ii,   inmet_smn_iii,   inmet_smn_iv,   inmet_smn_v   = [], [], [], [], []
era5_i,        era5_ii,        era5_iii,        era5_iv,        era5_v        = [], [], [], [], []
reg_usp_i,     reg_usp_ii,     reg_usp_iii,     reg_usp_iv,     reg_usp_v     = [], [], [], [], []
reg_ictp_i,    reg_ictp_ii,    reg_ictp_iii,    reg_ictp_iv,    reg_ictp_v    = [], [], [], [], []
reg_ictp_i_i,  reg_ictp_i_ii,  reg_ictp_i_iii,  reg_ictp_i_iv,  reg_ictp_i_v  = [], [], [], [], []
reg_ictp_ii_i, reg_ictp_ii_ii, reg_ictp_ii_iii, reg_ictp_ii_iv, reg_ictp_ii_v = [], [], [], [], []
wrf_ncar_i,    wrf_ncar_ii,    wrf_ncar_iii,    wrf_ncar_iv,    wrf_ncar_v    = [], [], [], [], []
wrf_ucan_i,    wrf_ucan_ii,    wrf_ucan_iii,    wrf_ucan_iv,    wrf_ucan_v    = [], [], [], [], []

for c_i in count_i:
	inmet_smn_i.append(inmet_smn[c_i])
	era5_i.append(era5[c_i])
	reg_usp_i.append(reg_usp[c_i])
	reg_ictp_i.append(reg_ictp[c_i])
	reg_ictp_i_i.append(reg_ictp_i_[c_i])
	reg_ictp_ii_i.append(reg_ictp_ii_[c_i])
	wrf_ncar_i.append(wrf_ncar[c_i])
	wrf_ucan_i.append(wrf_ucan[c_i])

for c_ii in count_ii:
	inmet_smn_ii.append(inmet_smn[c_ii])
	era5_ii.append(era5[c_ii])
	reg_usp_ii.append(reg_usp[c_ii])
	reg_ictp_ii.append(reg_ictp[c_ii])
	reg_ictp_i_ii.append(reg_ictp_i_[c_ii])
	reg_ictp_ii_ii.append(reg_ictp_ii_[c_ii])
	wrf_ncar_ii.append(wrf_ncar[c_ii])
	wrf_ucan_ii.append(wrf_ucan[c_ii])
	
for c_iii in count_iii:
	inmet_smn_iii.append(inmet_smn[c_iii])
	era5_iii.append(era5[c_iii])
	reg_usp_iii.append(reg_usp[c_iii])
	reg_ictp_iii.append(reg_ictp[c_iii])
	reg_ictp_i_iii.append(reg_ictp_i_[c_iii])
	reg_ictp_ii_iii.append(reg_ictp_ii_[c_iii])
	wrf_ncar_iii.append(wrf_ncar[c_iii])
	wrf_ucan_iii.append(wrf_ucan[c_iii])
	
for c_iv in count_iv:
	inmet_smn_iv.append(inmet_smn[c_iv])
	era5_iv.append(era5[c_iv])
	reg_usp_iv.append(reg_usp[c_iv])
	reg_ictp_iv.append(reg_ictp[c_iv])
	reg_ictp_i_iv.append(reg_ictp_i_[c_iv])
	reg_ictp_ii_iv.append(reg_ictp_ii_[c_iv])
	wrf_ncar_iv.append(wrf_ncar[c_iv])
	wrf_ucan_iv.append(wrf_ucan[c_iv])
	
for c_v in count_v:
	inmet_smn_v.append(inmet_smn[c_v])
	era5_v.append(era5[c_v])
	reg_usp_v.append(reg_usp[c_v])
	reg_ictp_v.append(reg_ictp[c_v])
	reg_ictp_i_v.append(reg_ictp_i_[c_v])
	reg_ictp_ii_v.append(reg_ictp_ii_[c_v])
	wrf_ncar_v.append(wrf_ncar[c_v])
	wrf_ucan_v.append(wrf_ucan[c_v])
	
# Group I
inmet_smn_c_i   = np.concatenate(inmet_smn_i)
inmet_smn_c_ii  = np.concatenate(inmet_smn_ii)
inmet_smn_c_iii = np.concatenate(inmet_smn_iii)
inmet_smn_c_iv  = np.concatenate(inmet_smn_iv)
inmet_smn_c_v   = np.concatenate(inmet_smn_v)

# Group II
era5_c_i_   = np.concatenate(era5_i)
era5_c_ii_  = np.concatenate(era5_ii)
era5_c_iii_ = np.concatenate(era5_iii)
era5_c_iv_  = np.concatenate(era5_iv)
era5_c_v_   = np.concatenate(era5_v)

era5_c_i = mask_like(inmet_smn_c_i, era5_c_i_)
era5_c_ii = mask_like(inmet_smn_c_ii, era5_c_ii_)
era5_c_iii = mask_like(inmet_smn_c_iii, era5_c_iii_)
era5_c_iv = mask_like(inmet_smn_c_iv, era5_c_iv_)
era5_c_v = mask_like(inmet_smn_c_v, era5_c_v_)

# Group III
reg_usp_c_i_   = np.concatenate(reg_usp_i)
reg_usp_c_ii_  = np.concatenate(reg_usp_ii)
reg_usp_c_iii_ = np.concatenate(reg_usp_iii)
reg_usp_c_iv_  = np.concatenate(reg_usp_iv)
reg_usp_c_v_   = np.concatenate(reg_usp_v)

reg_usp_c_i = mask_like(inmet_smn_c_i, reg_usp_c_i_)
reg_usp_c_ii = mask_like(inmet_smn_c_ii, reg_usp_c_ii_)
reg_usp_c_iii = mask_like(inmet_smn_c_iii, reg_usp_c_iii_)
reg_usp_c_iv = mask_like(inmet_smn_c_iv, reg_usp_c_iv_)
reg_usp_c_v = mask_like(inmet_smn_c_v, reg_usp_c_v_)

# Group IV
reg_ictp_c_i_   = np.concatenate(reg_ictp_i)
reg_ictp_c_ii_  = np.concatenate(reg_ictp_ii)
reg_ictp_c_iii_ = np.concatenate(reg_ictp_iii)
reg_ictp_c_iv_  = np.concatenate(reg_ictp_iv)
reg_ictp_c_v_   = np.concatenate(reg_ictp_v)

reg_ictp_c_i = mask_like(inmet_smn_c_i, reg_ictp_c_i_)
reg_ictp_c_ii = mask_like(inmet_smn_c_ii, reg_ictp_c_ii_)
reg_ictp_c_iii = mask_like(inmet_smn_c_iii, reg_ictp_c_iii_)
reg_ictp_c_iv = mask_like(inmet_smn_c_iv, reg_ictp_c_iv_)
reg_ictp_c_v = mask_like(inmet_smn_c_v, reg_ictp_c_v_)

# Group V
reg_ictp_i_c_i_   = np.concatenate(reg_ictp_i_i)
reg_ictp_i_c_ii_  = np.concatenate(reg_ictp_i_ii)
reg_ictp_i_c_iii_ = np.concatenate(reg_ictp_i_iii)
reg_ictp_i_c_iv_  = np.concatenate(reg_ictp_i_iv)
reg_ictp_i_c_v_   = np.concatenate(reg_ictp_i_v)

reg_ictp_i_c_i = mask_like(inmet_smn_c_i, reg_ictp_i_c_i_)
reg_ictp_i_c_ii = mask_like(inmet_smn_c_ii, reg_ictp_i_c_ii_)
reg_ictp_i_c_iii = mask_like(inmet_smn_c_iii, reg_ictp_i_c_iii_)
reg_ictp_i_c_iv = mask_like(inmet_smn_c_iv, reg_ictp_i_c_iv_)
reg_ictp_i_c_v = mask_like(inmet_smn_c_v, reg_ictp_i_c_v_)

# Group VI
reg_ictp_ii_c_i_   = np.concatenate(reg_ictp_ii_i)
reg_ictp_ii_c_ii_  = np.concatenate(reg_ictp_ii_ii)
reg_ictp_ii_c_iii_ = np.concatenate(reg_ictp_ii_iii)
reg_ictp_ii_c_iv_  = np.concatenate(reg_ictp_ii_iv)
reg_ictp_ii_c_v_   = np.concatenate(reg_ictp_ii_v)

reg_ictp_ii_c_i = mask_like(inmet_smn_c_i, reg_ictp_ii_c_i_)
reg_ictp_ii_c_ii = mask_like(inmet_smn_c_ii, reg_ictp_ii_c_ii_)
reg_ictp_ii_c_iii = mask_like(inmet_smn_c_iii, reg_ictp_ii_c_iii_)
reg_ictp_ii_c_iv = mask_like(inmet_smn_c_iv, reg_ictp_ii_c_iv_)
reg_ictp_ii_c_v = mask_like(inmet_smn_c_v, reg_ictp_ii_c_v_)

# Group VII
wrf_ncar_c_i_   = np.concatenate(wrf_ncar_i)
wrf_ncar_c_ii_  = np.concatenate(wrf_ncar_ii)
wrf_ncar_c_iii_ = np.concatenate(wrf_ncar_iii)
wrf_ncar_c_iv_  = np.concatenate(wrf_ncar_iv)
wrf_ncar_c_v_   = np.concatenate(wrf_ncar_v)

wrf_ncar_c_i = mask_like(inmet_smn_c_i, wrf_ncar_c_i_)
wrf_ncar_c_ii = mask_like(inmet_smn_c_ii, wrf_ncar_c_ii_)
wrf_ncar_c_iii = mask_like(inmet_smn_c_iii, wrf_ncar_c_iii_)
wrf_ncar_c_iv = mask_like(inmet_smn_c_iv, wrf_ncar_c_iv_)
wrf_ncar_c_v = mask_like(inmet_smn_c_v, wrf_ncar_c_v_)

# Group VIII
wrf_ucan_c_i_   = np.concatenate(wrf_ucan_i)
wrf_ucan_c_ii_  = np.concatenate(wrf_ucan_ii)
wrf_ucan_c_iii_ = np.concatenate(wrf_ucan_iii)
wrf_ucan_c_iv_  = np.concatenate(wrf_ucan_iv)
wrf_ucan_c_v_   = np.concatenate(wrf_ucan_v)

wrf_ucan_c_i = mask_like(inmet_smn_c_i, wrf_ucan_c_i_)
wrf_ucan_c_ii = mask_like(inmet_smn_c_ii, wrf_ucan_c_ii_)
wrf_ucan_c_iii = mask_like(inmet_smn_c_iii, wrf_ucan_c_iii_)
wrf_ucan_c_iv = mask_like(inmet_smn_c_iv, wrf_ucan_c_iv_)
wrf_ucan_c_v = mask_like(inmet_smn_c_v,wrf_ucan_c_v_)

perc_inmet_smn_c_i,   max_inmet_smn_c_i,   pdf_inmet_smn_c_i,   bins_inmet_smn_c_i   = compute_pdf(inmet_smn_c_i)
perc_era5_c_i,        max_era5_c_i,        pdf_era5_c_i,        bins_era5_c_i        = compute_pdf(era5_c_i)
perc_reg_usp_c_i,     max_reg_usp_c_i,     pdf_reg_usp_c_i,     bins_reg_usp_c_i     = compute_pdf(reg_usp_c_i)
perc_reg_ictp_c_i,    max_reg_ictp_c_i,    pdf_reg_ictp_c_i,    bins_reg_ictp_c_i    = compute_pdf(reg_ictp_c_i)
perc_reg_ictp_i_c_i,  max_reg_ictp_i_c_i,  pdf_reg_ictp_i_c_i,  bins_reg_ictp_i_c_i  = compute_pdf(reg_ictp_i_c_i)
perc_reg_ictp_ii_c_i, max_reg_ictp_ii_c_i, pdf_reg_ictp_ii_c_i, bins_reg_ictp_ii_c_i = compute_pdf(reg_ictp_ii_c_i)
perc_wrf_ncar_c_i,    max_wrf_ncar_c_i,    pdf_wrf_ncar_c_i,    bins_wrf_ncar_c_i    = compute_pdf(wrf_ncar_c_i)
perc_wrf_ucan_c_i,    max_wrf_ucan_c_i,    pdf_wrf_ucan_c_i,    bins_wrf_ucan_c_i    = compute_pdf(wrf_ucan_c_i)

perc_inmet_smn_c_ii,   max_inmet_smn_c_ii,   pdf_inmet_smn_c_ii,   bins_inmet_smn_c_ii   = compute_pdf(inmet_smn_c_ii)
perc_era5_c_ii,        max_era5_c_ii,        pdf_era5_c_ii,        bins_era5_c_ii        = compute_pdf(era5_c_ii)
perc_reg_usp_c_ii,     max_reg_usp_c_ii,     pdf_reg_usp_c_ii,	   bins_reg_usp_c_ii     = compute_pdf(reg_usp_c_ii)
perc_reg_ictp_c_ii,    max_reg_ictp_c_ii,    pdf_reg_ictp_c_ii,    bins_reg_ictp_c_ii	 = compute_pdf(reg_ictp_c_ii)
perc_reg_ictp_i_c_ii,  max_reg_ictp_i_c_ii,  pdf_reg_ictp_i_c_ii,  bins_reg_ictp_i_c_ii  = compute_pdf(reg_ictp_i_c_ii)
perc_reg_ictp_ii_c_ii, max_reg_ictp_ii_c_ii, pdf_reg_ictp_ii_c_ii, bins_reg_ictp_ii_c_ii = compute_pdf(reg_ictp_ii_c_ii)
perc_wrf_ncar_c_ii,    max_wrf_ncar_c_ii,    pdf_wrf_ncar_c_ii,    bins_wrf_ncar_c_ii    = compute_pdf(wrf_ncar_c_ii)
perc_wrf_ucan_c_ii,    max_wrf_ucan_c_ii,    pdf_wrf_ucan_c_ii,    bins_wrf_ucan_c_ii    = compute_pdf(wrf_ucan_c_ii)

perc_inmet_smn_c_iii,   max_inmet_smn_c_iii,   pdf_inmet_smn_c_iii,   bins_inmet_smn_c_iii   = compute_pdf(inmet_smn_c_iii)
perc_era5_c_iii,	max_era5_c_iii,        pdf_era5_c_iii,        bins_era5_c_iii	     = compute_pdf(era5_c_iii)
perc_reg_usp_c_iii,     max_reg_usp_c_iii,     pdf_reg_usp_c_iii,     bins_reg_usp_c_iii     = compute_pdf(reg_usp_c_iii)
perc_reg_ictp_c_iii,    max_reg_ictp_c_iii,    pdf_reg_ictp_c_iii,    bins_reg_ictp_c_iii    = compute_pdf(reg_ictp_c_iii)
perc_reg_ictp_i_c_iii,  max_reg_ictp_i_c_iii,  pdf_reg_ictp_i_c_iii,  bins_reg_ictp_i_c_iii  = compute_pdf(reg_ictp_i_c_iii)
perc_reg_ictp_ii_c_iii, max_reg_ictp_ii_c_iii, pdf_reg_ictp_ii_c_iii, bins_reg_ictp_ii_c_iii = compute_pdf(reg_ictp_ii_c_iii)
perc_wrf_ncar_c_iii,    max_wrf_ncar_c_iii,    pdf_wrf_ncar_c_iii,    bins_wrf_ncar_c_iii    = compute_pdf(wrf_ncar_c_iii)
perc_wrf_ucan_c_iii,    max_wrf_ucan_c_iii,    pdf_wrf_ucan_c_iii,    bins_wrf_ucan_c_iii    = compute_pdf(wrf_ucan_c_iii)

perc_inmet_smn_c_iv,   max_inmet_smn_c_iv,   pdf_inmet_smn_c_iv,   bins_inmet_smn_c_iv   = compute_pdf(inmet_smn_c_iv)
perc_era5_c_iv,        max_era5_c_iv,        pdf_era5_c_iv,        bins_era5_c_iv        = compute_pdf(era5_c_iv)
perc_reg_usp_c_iv,     max_reg_usp_c_iv,     pdf_reg_usp_c_iv,     bins_reg_usp_c_iv	 = compute_pdf(reg_usp_c_iv)
perc_reg_ictp_c_iv,    max_reg_ictp_c_iv,    pdf_reg_ictp_c_iv,    bins_reg_ictp_c_iv    = compute_pdf(reg_ictp_c_iv)
perc_reg_ictp_i_c_iv,  max_reg_ictp_i_c_iv,  pdf_reg_ictp_i_c_iv,  bins_reg_ictp_i_c_iv  = compute_pdf(reg_ictp_i_c_iv)
perc_reg_ictp_ii_c_iv, max_reg_ictp_ii_c_iv, pdf_reg_ictp_ii_c_iv, bins_reg_ictp_ii_c_iv = compute_pdf(reg_ictp_ii_c_iv)
perc_wrf_ncar_c_iv,    max_wrf_ncar_c_iv,    pdf_wrf_ncar_c_iv,    bins_wrf_ncar_c_iv	 = compute_pdf(wrf_ncar_c_iv)
perc_wrf_ucan_c_iv,    max_wrf_ucan_c_iv,    pdf_wrf_ucan_c_iv,    bins_wrf_ucan_c_iv    = compute_pdf(wrf_ucan_c_iv)

perc_inmet_smn_c_v,   max_inmet_smn_c_v,   pdf_inmet_smn_c_v,	bins_inmet_smn_c_v   = compute_pdf(inmet_smn_c_v)
perc_era5_c_v,        max_era5_c_v,        pdf_era5_c_v,        bins_era5_c_v        = compute_pdf(era5_c_v)
perc_reg_usp_c_v,     max_reg_usp_c_v,     pdf_reg_usp_c_v,	bins_reg_usp_c_v     = compute_pdf(reg_usp_c_v)
perc_reg_ictp_c_v,    max_reg_ictp_c_v,    pdf_reg_ictp_c_v,	bins_reg_ictp_c_v    = compute_pdf(reg_ictp_c_v)
perc_reg_ictp_i_c_v,  max_reg_ictp_i_c_v,  pdf_reg_ictp_i_c_v,  bins_reg_ictp_i_c_v  = compute_pdf(reg_ictp_i_c_v)
perc_reg_ictp_ii_c_v, max_reg_ictp_ii_c_v, pdf_reg_ictp_ii_c_v, bins_reg_ictp_ii_c_v = compute_pdf(reg_ictp_ii_c_v)
perc_wrf_ncar_c_v,    max_wrf_ncar_c_v,    pdf_wrf_ncar_c_v,	bins_wrf_ncar_c_v    = compute_pdf(wrf_ncar_c_v)
perc_wrf_ucan_c_v,    max_wrf_ucan_c_v,    pdf_wrf_ucan_c_v,	bins_wrf_ucan_c_v    = compute_pdf(wrf_ucan_c_v)

# Plot figure
fig = plt.figure(figsize=(12, 8))
font_size = 8

legend = 'mm d⁻¹'
xvmin = 0
xvmax = 300
yvmin = 1e-4
yvmax = 0
text_ = 0.65
text_1 = 0.45
text_2 = 0.4
text_3 = 0.35
text_4 = 0.3
text_5 = 0.25
text_6 = 0.2
text_7 = 0.15
text_8 = 0.1
		
ax = fig.add_subplot(2, 3, 1)
plt.plot(pdf_inmet_smn_c_i,   bins_inmet_smn_c_i,   marker='.', markersize=4, markerfacecolor='None', markeredgecolor='black',   markeredgewidth=0.5, linestyle='None', label='WS')
plt.plot(pdf_era5_c_i,        bins_era5_c_i,        marker='.', markersize=4, markerfacecolor='None', markeredgecolor='red',     markeredgewidth=0.5, linestyle='None', label='ERA5')
plt.plot(pdf_reg_usp_c_i,     bins_reg_usp_c_i,     marker='.', markersize=4, markerfacecolor='None', markeredgecolor='blue',    markeredgewidth=0.5, linestyle='None', label='Reg4')
plt.plot(pdf_reg_ictp_c_i,    bins_reg_ictp_c_i,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='magenta', markeredgewidth=0.5, linestyle='None', label='Reg5-Holt3')
plt.plot(pdf_reg_ictp_i_c_i,  bins_reg_ictp_i_c_i,  marker='.', markersize=4, markerfacecolor='None', markeredgecolor='gray',    markeredgewidth=0.5, linestyle='None', label='Reg5-Holt')
plt.plot(pdf_reg_ictp_ii_c_i, bins_reg_ictp_ii_c_i, marker='.', markersize=4, markerfacecolor='None', markeredgecolor='brown',   markeredgewidth=0.5, linestyle='None', label='Reg5-UW')
plt.plot(pdf_wrf_ncar_c_i,    bins_wrf_ncar_c_i,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='green',   markeredgewidth=0.5, linestyle='None', label='WRF-NCAR')
plt.plot(pdf_wrf_ucan_c_i,    bins_wrf_ucan_c_i,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='orange',  markeredgewidth=0.5, linestyle='None', label='WRF-UCAN')
plt.axvline(perc_inmet_smn_c_i,   color='black',   linestyle='-', linewidth=0.75)
plt.axvline(perc_era5_c_i,        color='red',     linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_usp_c_i,     color='blue',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_c_i,    color='magenta', linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_i_c_i,  color='gray',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_ii_c_i, color='brown',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ncar_c_i,    color='green',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ucan_c_i,    color='orange',  linestyle='-', linewidth=0.75)
plt.text(text_, text_1, r'$Pr_{{\max}} = {0}$'.format(round(max_inmet_smn_c_i, 1)),   transform=ax.transAxes, color='black',   fontsize=font_size)
plt.text(text_, text_2, r'$Pr_{{\max}} = {0}$'.format(round(max_era5_c_i, 1)),        transform=ax.transAxes, color='red',     fontsize=font_size)
plt.text(text_, text_3, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_usp_c_i, 1)),     transform=ax.transAxes, color='blue',    fontsize=font_size)
plt.text(text_, text_4, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_c_i, 1)),    transform=ax.transAxes, color='magenta', fontsize=font_size)
plt.text(text_, text_5, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_i_c_i, 1)),  transform=ax.transAxes, color='gray',    fontsize=font_size)
plt.text(text_, text_6, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_ii_c_i, 1)), transform=ax.transAxes, color='brown',   fontsize=font_size)
plt.text(text_, text_7, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ncar_c_i, 1)),    transform=ax.transAxes, color='green',   fontsize=font_size)
plt.text(text_, text_8, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ucan_c_i, 1)),    transform=ax.transAxes, color='orange',  fontsize=font_size)
plt.title('(a)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CLUSTER I', loc='center', fontsize=font_size)
plt.ylabel('PDF (#)', fontsize=font_size)
plt.yscale('log')
plt.ylim(yvmin, yvmax)
plt.xlim(xvmin, xvmax)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
ax.grid(True, which='major', axis='x', linestyle='--', linewidth=0.8)
ax.grid(True, which='major', axis='y', linestyle='--', linewidth=0.8)
ax.grid(True, which='minor', axis='y', linestyle='--', linewidth=0.4)
plt.legend(loc=1, fontsize=font_size, frameon=True, framealpha=0.9)

ax = fig.add_subplot(2, 3, 2)
plt.plot(pdf_inmet_smn_c_ii,   bins_inmet_smn_c_ii,   marker='.', markersize=4, markerfacecolor='None', markeredgecolor='black',   markeredgewidth=0.5, linestyle='None', label='WS')
plt.plot(pdf_era5_c_ii,        bins_era5_c_ii,        marker='.', markersize=4, markerfacecolor='None', markeredgecolor='red',     markeredgewidth=0.5, linestyle='None', label='ERA5')
plt.plot(pdf_reg_usp_c_ii,     bins_reg_usp_c_ii,     marker='.', markersize=4, markerfacecolor='None', markeredgecolor='blue',    markeredgewidth=0.5, linestyle='None', label='Reg4')
plt.plot(pdf_reg_ictp_c_ii,    bins_reg_ictp_c_ii,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='magenta', markeredgewidth=0.5, linestyle='None', label='Reg5-Holt3')
plt.plot(pdf_reg_ictp_i_c_ii,  bins_reg_ictp_i_c_ii,  marker='.', markersize=4, markerfacecolor='None', markeredgecolor='gray',    markeredgewidth=0.5, linestyle='None', label='Reg5-Holt')
plt.plot(pdf_reg_ictp_ii_c_ii, bins_reg_ictp_ii_c_ii, marker='.', markersize=4, markerfacecolor='None', markeredgecolor='brown',   markeredgewidth=0.5, linestyle='None', label='Reg5-UW')
plt.plot(pdf_wrf_ncar_c_ii,    bins_wrf_ncar_c_ii,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='green',   markeredgewidth=0.5, linestyle='None', label='WRF-NCAR')
plt.plot(pdf_wrf_ucan_c_ii,    bins_wrf_ucan_c_ii,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='orange',  markeredgewidth=0.5, linestyle='None', label='WRF-UCAN')
plt.axvline(perc_inmet_smn_c_ii,   color='black',   linestyle='-', linewidth=0.75)
plt.axvline(perc_era5_c_ii,        color='red',     linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_usp_c_ii,     color='blue',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_c_ii,    color='magenta', linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_i_c_ii,  color='gray',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_ii_c_ii, color='brown',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ncar_c_ii,    color='green',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ucan_c_ii,    color='orange',  linestyle='-', linewidth=0.75)
plt.text(text_, text_1, r'$Pr_{{\max}} = {0}$'.format(round(max_inmet_smn_c_ii, 1)),   transform=ax.transAxes, color='black',   fontsize=font_size)
plt.text(text_, text_2, r'$Pr_{{\max}} = {0}$'.format(round(max_era5_c_ii, 1)),        transform=ax.transAxes, color='red',     fontsize=font_size)
plt.text(text_, text_3, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_usp_c_ii, 1)),     transform=ax.transAxes, color='blue',    fontsize=font_size)
plt.text(text_, text_4, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_c_ii, 1)),    transform=ax.transAxes, color='magenta', fontsize=font_size)
plt.text(text_, text_5, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_i_c_ii, 1)),  transform=ax.transAxes, color='gray',    fontsize=font_size)
plt.text(text_, text_6, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_ii_c_ii, 1)), transform=ax.transAxes, color='brown',   fontsize=font_size)
plt.text(text_, text_7, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ncar_c_ii, 1)),    transform=ax.transAxes, color='green',   fontsize=font_size)
plt.text(text_, text_8, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ucan_c_ii, 1)),    transform=ax.transAxes, color='orange',  fontsize=font_size)
plt.title('(b)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CLUSTER II', loc='center', fontsize=font_size)
plt.yscale('log')
plt.ylim(yvmin, yvmax)
plt.xlim(xvmin, xvmax)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
ax.grid(True, which='major', axis='x', linestyle='--', linewidth=0.8)
ax.grid(True, which='major', axis='y', linestyle='--', linewidth=0.8)
ax.grid(True, which='minor', axis='y', linestyle='--', linewidth=0.4)

ax = fig.add_subplot(2, 3, 3)
plt.plot(pdf_inmet_smn_c_iii,   bins_inmet_smn_c_iii,   marker='.', markersize=4, markerfacecolor='None', markeredgecolor='black',   markeredgewidth=0.5, linestyle='None', label='WS')
plt.plot(pdf_era5_c_iii,        bins_era5_c_iii,        marker='.', markersize=4, markerfacecolor='None', markeredgecolor='red',     markeredgewidth=0.5, linestyle='None', label='ERA5')
plt.plot(pdf_reg_usp_c_iii,     bins_reg_usp_c_iii,     marker='.', markersize=4, markerfacecolor='None', markeredgecolor='blue',    markeredgewidth=0.5, linestyle='None', label='Reg4')
plt.plot(pdf_reg_ictp_c_iii,    bins_reg_ictp_c_iii,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='magenta', markeredgewidth=0.5, linestyle='None', label='Reg5-Holt3')
plt.plot(pdf_reg_ictp_i_c_iii,  bins_reg_ictp_i_c_iii,  marker='.', markersize=4, markerfacecolor='None', markeredgecolor='gray',    markeredgewidth=0.5, linestyle='None', label='Reg5-Holt')
plt.plot(pdf_reg_ictp_ii_c_iii, bins_reg_ictp_ii_c_iii, marker='.', markersize=4, markerfacecolor='None', markeredgecolor='brown',   markeredgewidth=0.5, linestyle='None', label='Reg5-UW')
plt.plot(pdf_wrf_ncar_c_iii,    bins_wrf_ncar_c_iii,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='green',   markeredgewidth=0.5, linestyle='None', label='WRF-NCAR')
plt.plot(pdf_wrf_ucan_c_iii,    bins_wrf_ucan_c_iii,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='orange',  markeredgewidth=0.5, linestyle='None', label='WRF-UCAN')
plt.axvline(perc_inmet_smn_c_iii,   color='black',   linestyle='-', linewidth=0.75)
plt.axvline(perc_era5_c_iii,        color='red',     linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_usp_c_iii,     color='blue',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_c_iii,    color='magenta', linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_i_c_iii,  color='gray',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_ii_c_iii, color='brown',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ncar_c_iii,    color='green',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ucan_c_iii,    color='orange',  linestyle='-', linewidth=0.75)
plt.text(text_, text_1, r'$Pr_{{\max}} = {0}$'.format(round(max_inmet_smn_c_iii, 1)),   transform=ax.transAxes, color='black',   fontsize=font_size)
plt.text(text_, text_2, r'$Pr_{{\max}} = {0}$'.format(round(max_era5_c_iii, 1)),        transform=ax.transAxes, color='red',     fontsize=font_size)
plt.text(text_, text_3, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_usp_c_iii, 1)),     transform=ax.transAxes, color='blue',    fontsize=font_size)
plt.text(text_, text_4, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_c_iii, 1)),    transform=ax.transAxes, color='magenta', fontsize=font_size)
plt.text(text_, text_5, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_i_c_iii, 1)),  transform=ax.transAxes, color='gray',    fontsize=font_size)
plt.text(text_, text_6, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_ii_c_iii, 1)), transform=ax.transAxes, color='brown',   fontsize=font_size)
plt.text(text_, text_7, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ncar_c_iii, 1)),    transform=ax.transAxes, color='green',   fontsize=font_size)
plt.text(text_, text_8, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ucan_c_iii, 1)),    transform=ax.transAxes, color='orange',  fontsize=font_size)
plt.title('(c)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CLUSTER III', loc='center', fontsize=font_size)
plt.xlabel('{0}'.format(legend), fontsize=font_size)
plt.yscale('log')
plt.ylim(yvmin, yvmax)
plt.xlim(xvmin, xvmax)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
ax.grid(True, which='major', axis='x', linestyle='--', linewidth=0.8)
ax.grid(True, which='major', axis='y', linestyle='--', linewidth=0.8)
ax.grid(True, which='minor', axis='y', linestyle='--', linewidth=0.4)

ax = fig.add_subplot(2, 3, 4)
plt.plot(pdf_inmet_smn_c_iv,   bins_inmet_smn_c_iv,   marker='.', markersize=4, markerfacecolor='None', markeredgecolor='black',   markeredgewidth=0.5, linestyle='None', label='WS')
plt.plot(pdf_era5_c_iv,        bins_era5_c_iv,        marker='.', markersize=4, markerfacecolor='None', markeredgecolor='red',     markeredgewidth=0.5, linestyle='None', label='ERA5')
plt.plot(pdf_reg_usp_c_iv,     bins_reg_usp_c_iv,     marker='.', markersize=4, markerfacecolor='None', markeredgecolor='blue',    markeredgewidth=0.5, linestyle='None', label='Reg4')
plt.plot(pdf_reg_ictp_c_iv,    bins_reg_ictp_c_iv,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='magenta', markeredgewidth=0.5, linestyle='None', label='Reg5-Holt3')
plt.plot(pdf_reg_ictp_i_c_iv,  bins_reg_ictp_i_c_iv,  marker='.', markersize=4, markerfacecolor='None', markeredgecolor='gray',    markeredgewidth=0.5, linestyle='None', label='Reg5-Holt')
plt.plot(pdf_reg_ictp_ii_c_iv, bins_reg_ictp_ii_c_iv, marker='.', markersize=4, markerfacecolor='None', markeredgecolor='brown',   markeredgewidth=0.5, linestyle='None', label='Reg5-UW')
plt.plot(pdf_wrf_ncar_c_iv,    bins_wrf_ncar_c_iv,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='green',   markeredgewidth=0.5, linestyle='None', label='WRF-NCAR')
plt.plot(pdf_wrf_ucan_c_iv,    bins_wrf_ucan_c_iv,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='orange',  markeredgewidth=0.5, linestyle='None', label='WRF-UCAN')
plt.axvline(perc_inmet_smn_c_iv,   color='black',   linestyle='-', linewidth=0.75)
plt.axvline(perc_era5_c_iv,        color='red',     linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_usp_c_iv,     color='blue',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_c_iv,    color='magenta', linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_i_c_iv,  color='gray',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_ii_c_iv, color='brown',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ncar_c_iv,    color='green',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ucan_c_iv,    color='orange',  linestyle='-', linewidth=0.75)
plt.text(text_, text_1, r'$Pr_{{\max}} = {0}$'.format(round(max_inmet_smn_c_iv, 1)),   transform=ax.transAxes, color='black',   fontsize=font_size)
plt.text(text_, text_2, r'$Pr_{{\max}} = {0}$'.format(round(max_era5_c_iv, 1)),        transform=ax.transAxes, color='red',     fontsize=font_size)
plt.text(text_, text_3, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_usp_c_iv, 1)),     transform=ax.transAxes, color='blue',    fontsize=font_size)
plt.text(text_, text_4, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_c_iv, 1)),    transform=ax.transAxes, color='magenta', fontsize=font_size)
plt.text(text_, text_5, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_i_c_iv, 1)),  transform=ax.transAxes, color='gray',    fontsize=font_size)
plt.text(text_, text_6, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_ii_c_iv, 1)), transform=ax.transAxes, color='brown',   fontsize=font_size)
plt.text(text_, text_7, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ncar_c_iv, 1)),    transform=ax.transAxes, color='green',   fontsize=font_size)
plt.text(text_, text_8, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ucan_c_iv, 1)),    transform=ax.transAxes, color='orange',  fontsize=font_size)
plt.title('(d)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CLUSTER IV', loc='center', fontsize=font_size)
plt.ylabel('PDF (#)', fontsize=font_size)
plt.xlabel('{0}'.format(legend), fontsize=font_size)
plt.yscale('log')
plt.ylim(yvmin, yvmax)
plt.xlim(xvmin, xvmax)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
ax.grid(True, which='major', axis='x', linestyle='--', linewidth=0.8)
ax.grid(True, which='major', axis='y', linestyle='--', linewidth=0.8)
ax.grid(True, which='minor', axis='y', linestyle='--', linewidth=0.4)

ax = fig.add_subplot(2, 3, 5)
plt.plot(pdf_inmet_smn_c_v,   bins_inmet_smn_c_v,   marker='.', markersize=4, markerfacecolor='None', markeredgecolor='black',   markeredgewidth=0.5, linestyle='None', label='WS')
plt.plot(pdf_era5_c_v,        bins_era5_c_v,        marker='.', markersize=4, markerfacecolor='None', markeredgecolor='red',     markeredgewidth=0.5, linestyle='None', label='ERA5')
plt.plot(pdf_reg_usp_c_v,     bins_reg_usp_c_v,     marker='.', markersize=4, markerfacecolor='None', markeredgecolor='blue',    markeredgewidth=0.5, linestyle='None', label='Reg4')
plt.plot(pdf_reg_ictp_c_v,    bins_reg_ictp_c_v,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='magenta', markeredgewidth=0.5, linestyle='None', label='Reg5-Holt3')
plt.plot(pdf_reg_ictp_i_c_v,  bins_reg_ictp_i_c_v,  marker='.', markersize=4, markerfacecolor='None', markeredgecolor='gray',    markeredgewidth=0.5, linestyle='None', label='Reg5-Holt')
plt.plot(pdf_reg_ictp_ii_c_v, bins_reg_ictp_ii_c_v, marker='.', markersize=4, markerfacecolor='None', markeredgecolor='brown',   markeredgewidth=0.5, linestyle='None', label='Reg5-UW')
plt.plot(pdf_wrf_ncar_c_v,    bins_wrf_ncar_c_v,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='green',   markeredgewidth=0.5, linestyle='None', label='WRF-NCAR')
plt.plot(pdf_wrf_ucan_c_v,    bins_wrf_ucan_c_v,    marker='.', markersize=4, markerfacecolor='None', markeredgecolor='orange',  markeredgewidth=0.5, linestyle='None', label='WRF-UCAN')
plt.axvline(perc_inmet_smn_c_v,   color='black',   linestyle='-', linewidth=0.75)
plt.axvline(perc_era5_c_v,        color='red',     linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_usp_c_v,     color='blue',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_c_v,    color='magenta', linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_i_c_v,  color='gray',    linestyle='-', linewidth=0.75)
plt.axvline(perc_reg_ictp_ii_c_v, color='brown',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ncar_c_v,    color='green',   linestyle='-', linewidth=0.75)
plt.axvline(perc_wrf_ucan_c_v,    color='orange',  linestyle='-', linewidth=0.75)
plt.text(text_, text_1, r'$Pr_{{\max}} = {0}$'.format(round(max_inmet_smn_c_v, 1)),   transform=ax.transAxes, color='black',   fontsize=font_size)
plt.text(text_, text_2, r'$Pr_{{\max}} = {0}$'.format(round(max_era5_c_v, 1)),        transform=ax.transAxes, color='red',     fontsize=font_size)
plt.text(text_, text_3, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_usp_c_v, 1)),     transform=ax.transAxes, color='blue',    fontsize=font_size)
plt.text(text_, text_4, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_c_v, 1)),    transform=ax.transAxes, color='magenta', fontsize=font_size)
plt.text(text_, text_5, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_i_c_v, 1)),  transform=ax.transAxes, color='gray',    fontsize=font_size)
plt.text(text_, text_6, r'$Pr_{{\max}} = {0}$'.format(round(max_reg_ictp_ii_c_v, 1)), transform=ax.transAxes, color='brown',   fontsize=font_size)
plt.text(text_, text_7, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ncar_c_v, 1)),    transform=ax.transAxes, color='green',   fontsize=font_size)
plt.text(text_, text_8, r'$Pr_{{\max}} = {0}$'.format(round(max_wrf_ucan_c_v, 1)),    transform=ax.transAxes, color='orange',  fontsize=font_size)
plt.title('(e)', loc='left', fontsize=font_size, fontweight='bold')
plt.title('CLUSTER V', loc='center', fontsize=font_size)
plt.xlabel('{0}'.format(legend), fontsize=font_size)
plt.yscale('log')
plt.ylim(yvmin, yvmax)
plt.xlim(xvmin, xvmax)
plt.xticks(fontsize=font_size)
plt.yticks(fontsize=font_size)
ax.grid(True, which='major', axis='x', linestyle='--', linewidth=0.8)
ax.grid(True, which='major', axis='y', linestyle='--', linewidth=0.8)
ax.grid(True, which='minor', axis='y', linestyle='--', linewidth=0.4)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_graph_pdf_stats_{0}_mmd_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
