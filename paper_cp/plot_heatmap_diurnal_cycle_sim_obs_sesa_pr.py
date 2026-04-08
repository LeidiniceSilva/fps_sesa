# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Dec 08, 2025"
__description__ = "This script plot diurnal cycle"

import os
import argparse
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from scipy.stats import gaussian_kde
from matplotlib.patches import Polygon
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i

parser = argparse.ArgumentParser(description='Process variable')
parser.add_argument('--var', required=True, choices=['pr'], help='Variable name')
args = parser.parse_args()
var = args.var

path = '/home/mda_silv/users/FPS_SESA'

skip_list_inmet_i = [15,23,47,105,112,117,124,137,149,158,174,183,335,343,359,398,399,413,417,422,426,444,453,457,458,479,490,495,505,529,566] 
	
skip_list_inmet_ii = [2, 3, 4, 14, 19, 20, 21, 24, 25, 26, 27, 28, 32, 33, 34, 35, 38, 40, 41, 44, 45, 48, 52, 54, 55, 56, 59, 60, 62, 64, 68, 
70, 77, 79, 80, 82, 83, 92, 93, 96, 100, 106, 107, 111, 113, 120, 127, 130, 133, 135, 136, 140, 141, 144, 152, 154, 155, 160, 161, 163, 167, 168, 
173, 177, 180, 181, 182, 184, 186, 187, 188, 193, 197, 199, 204, 206, 207, 210, 212, 215, 216, 219, 220, 224, 225, 226, 229, 233, 237, 239, 240, 
241, 243, 248, 249, 251, 253, 254, 256, 261, 262, 264, 266, 269, 275, 276, 277, 280, 281, 282, 293, 295, 296, 298, 300, 303, 306, 308, 314, 315, 
316, 317, 319, 322, 325, 330, 331, 334, 337, 341, 344, 347, 348, 350, 353, 354, 357, 358, 360, 361, 362, 364, 370, 383, 384, 385, 389, 390, 392, 
393, 395, 396, 400, 401, 402, 404, 405, 408, 415, 416, 418, 423, 424, 427, 434, 440, 441, 443, 446, 448, 450, 451, 454, 455, 459, 465, 467, 471, 
474, 477, 481, 483, 488, 489, 492, 496, 504, 509, 513, 514, 516, 518, 519, 520, 523, 526, 528, 534, 538, 541, 544, 546, 552, 553, 557, 559]


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
			print(station_code, station_name)
		
			# Reading inmet 
			d_i = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/hourly/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(station_code))
			d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.values
			mean_i.append(np.concatenate([d_i[3:], d_i[:3]]))

			# Reading era5 
			d_ii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/tp/' + 'tp_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
			d_ii = d_ii.tp.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.values 
			mean_ii.append(np.concatenate([d_ii[3:], d_ii[:3]]))

			# Reading regcm usp
			d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
			d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.values / 24
			mean_iii.append(np.concatenate([d_iii[3:], d_iii[:3]]))
			
			# Reading regcm ictp 
			d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
			d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.values 
			mean_iv.append(np.concatenate([d_iv[3:], d_iv[:3]]))
	
			# Reading regcm ictp pbl 1
			d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
			d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.values / 24
			mean_v.append(np.concatenate([d_v[3:], d_v[:3]]))

			# Reading regcm ictp pbl 2
			d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
			d_vi = d_vi.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.values / 24
			mean_vi.append(np.concatenate([d_vi[3:], d_vi[:3]]))

			# Reading wrf ncar
			d_vii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ncar/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
			d_vii = d_vii.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.values / 24
			mean_vii.append(np.concatenate([d_vii[3:], d_vii[:3]]))
		
			# Reading wrf ucan
			d_viii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ucan/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
			d_viii = d_viii.pr.sel(time=slice('2018-06-01','2021-05-31'))
			d_viii = d_viii.values / 24
			mean_viii.append(np.concatenate([d_viii[3:], d_viii[:3]]))
		
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii


def import_smn_i():

	mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii = [], [], [], [], [], [], [], []
	for i in range(1, 73):
		station_code = f'SMN{i:03d}'
		station_name = smn_i[i][0]
		print(station_code, station_name)
			
		# Reading smn 
		d_i = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(station_name))
		d_i = d_i.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.values
		mean_i.append(np.concatenate([d_i[3:], d_i[:3]]))

		# Reading era5 
		d_ii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/obs/era5/tp/' + 'tp_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
		d_ii = d_ii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.values 
		mean_ii.append(np.concatenate([d_ii[3:], d_ii[:3]]))

		# Reading regcm usp
		d_iii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_usp/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.values / 24
		mean_iii.append(np.concatenate([d_iii[3:], d_iii[:3]]))
			
		# Reading regcm ictp 
		d_iv = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.values 
		mean_iv.append(np.concatenate([d_iv[3:], d_iv[:3]]))
	
		# Reading regcm ictp pbl 1
		d_v = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl1/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.values / 24
		mean_v.append(np.concatenate([d_v[3:], d_v[:3]]))

		# Reading regcm ictp pbl 2
		d_vi = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/reg_ictp_pbl2/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
		d_vi = d_vi.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.values / 24
		mean_vi.append(np.concatenate([d_vi[3:], d_vi[:3]]))

		# Reading wrf ncar
		d_vii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ncar/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
		d_vii = d_vii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.values / 24
		mean_vii.append(np.concatenate([d_vii[3:], d_vii[:3]]))
		
		# Reading wrf ucan
		d_viii = xr.open_dataset('/home/mda_silv/clima-archive2-b/FPS-SESA/rcm/wrf_ucan/pr/' + 'pr_{0}_{1}_H_2018-06-01-2021-05-31.nc'.format(station_code, station_name))
		d_viii = d_viii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_viii = d_viii.values / 24
		mean_viii.append(np.concatenate([d_viii[3:], d_viii[:3]]))
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii, mean_viii
	

def mask_like(reference, target):

	reference = np.asarray(reference)
	masked_target = np.asarray(target, dtype=float).copy()
	masked_target[np.isnan(reference)] = np.nan
	
	return masked_target


def compute_stats(dataset, perc=99):

    var_ = np.asarray(dataset)

    # Remove non-finite values
    var_[~np.isfinite(var_)] = np.nan

    # Remove unrealistic precipitation values
    var_[var_ < 0] = np.nan       # negative precipitation doesn't exist
    var_[var_ > 500] = np.nan     # very high threshold (mm/h)

    # Hour index
    horas = np.arange(len(var_)) % 24

    # Global percentile (for extremes definition)
    perc_global = np.nanpercentile(var_, perc)

    # Prepare outputs
    mean_hour = np.zeros(24)
    perc_hour = np.zeros(24)
    freq = np.zeros(24)
    intens = np.zeros(24)

    for h in range(24):
        var_h = var_[horas == h]
        valid = ~np.isnan(var_h)

        if np.sum(valid) > 0:
            data = var_h[valid]

            mean_hour[h] = np.nanmean(data)                 # mean per hour
            perc_hour[h] = np.nanpercentile(data, perc)     # percentile per hour

            extreme = data[data >= perc_global]             # extremes based on global perc
            freq[h] = 100.0 * len(extreme) / len(data)      # frequency of extremes
            intens[h] = np.nanmean(extreme) if len(extreme) > 0 else np.nan

        else:
            mean_hour[h] = np.nan
            perc_hour[h] = np.nan
            freq[h] = np.nan
            intens[h] = np.nan

    return mean_hour, perc_hour, freq, intens
    
    
# Import dataset
clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x, clim_vi_x, clim_vii_x, clim_viii_x = import_inmet()			
clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y, clim_vi_y, clim_vii_y, clim_viii_y = import_smn_i()

inmet_smn    = clim_i_x    + clim_i_y
era5         = clim_ii_x   + clim_ii_y
reg_usp      = clim_iii_x  + clim_iii_y
reg_ictp     = clim_iv_x   + clim_iv_y
reg_ictp_i_  = clim_v_x    + clim_v_y
reg_ictp_ii_ = clim_vi_x   + clim_vi_y
wrf_ncar     = clim_vii_x  + clim_vii_y
wrf_ucan     = clim_viii_x + clim_viii_y

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
reg_ictp_i_c_v_  = np.concatenate(reg_ictp_i_v)

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

mean_inmet_smn_c_i,   perc_inmet_smn_c_i,   freq_inmet_smn_c_i,   int_inmet_smn_c_i   = compute_stats(inmet_smn_c_i)
mean_inmet_smn_c_ii,  perc_inmet_smn_c_ii,  freq_inmet_smn_c_ii,  int_inmet_smn_c_ii  = compute_stats(inmet_smn_c_ii)
mean_inmet_smn_c_iii, perc_inmet_smn_c_iii, freq_inmet_smn_c_iii, int_inmet_smn_c_iii = compute_stats(inmet_smn_c_iii)
mean_inmet_smn_c_iv,  perc_inmet_smn_c_iv,  freq_inmet_smn_c_iv,  int_inmet_smn_c_iv  = compute_stats(inmet_smn_c_iv)
mean_inmet_smn_c_v,   perc_inmet_smn_c_v,   freq_inmet_smn_c_v,   int_inmet_smn_c_v   = compute_stats(inmet_smn_c_v)

mean_era5_c_i,   perc_era5_c_i,   freq_era5_c_i,   int_era5_c_i   = compute_stats(era5_c_i)
mean_era5_c_ii,  perc_era5_c_ii,  freq_era5_c_ii,  int_era5_c_ii  = compute_stats(era5_c_ii)
mean_era5_c_iii, perc_era5_c_iii, freq_era5_c_iii, int_era5_c_iii = compute_stats(era5_c_iii)
mean_era5_c_iv,  perc_era5_c_iv,  freq_era5_c_iv,  int_era5_c_iv  = compute_stats(era5_c_iv)
mean_era5_c_v,   perc_era5_c_v,   freq_era5_c_v,   int_era5_c_v   = compute_stats(era5_c_v)

mean_reg_usp_c_i,   perc_reg_usp_c_i,   freq_reg_usp_c_i,   int_reg_usp_c_i   = compute_stats(reg_usp_c_i)
mean_reg_usp_c_ii,  perc_reg_usp_c_ii,  freq_reg_usp_c_ii,  int_reg_usp_c_ii  = compute_stats(reg_usp_c_ii)
mean_reg_usp_c_iii, perc_reg_usp_c_iii, freq_reg_usp_c_iii, int_reg_usp_c_iii = compute_stats(reg_usp_c_iii)
mean_reg_usp_c_iv,  perc_reg_usp_c_iv,  freq_reg_usp_c_iv,  int_reg_usp_c_iv  = compute_stats(reg_usp_c_iv)
mean_reg_usp_c_v,   perc_reg_usp_c_v,   freq_reg_usp_c_v,   int_reg_usp_c_v   = compute_stats(reg_usp_c_v)

mean_reg_ictp_c_i,   perc_reg_ictp_c_i,   freq_reg_ictp_c_i,   int_reg_ictp_c_i   = compute_stats(reg_ictp_c_i)
mean_reg_ictp_c_ii,  perc_reg_ictp_c_ii,  freq_reg_ictp_c_ii,  int_reg_ictp_c_ii  = compute_stats(reg_ictp_c_ii)
mean_reg_ictp_c_iii, perc_reg_ictp_c_iii, freq_reg_ictp_c_iii, int_reg_ictp_c_iii = compute_stats(reg_ictp_c_iii)
mean_reg_ictp_c_iv,  perc_reg_ictp_c_iv,  freq_reg_ictp_c_iv,  int_reg_ictp_c_iv  = compute_stats(reg_ictp_c_iv)
mean_reg_ictp_c_v,   perc_reg_ictp_c_v,   freq_reg_ictp_c_v,   int_reg_ictp_c_v   = compute_stats(reg_ictp_c_v)

mean_reg_ictp_i_c_i,   perc_reg_ictp_i_c_i,   freq_reg_ictp_i_c_i,   int_reg_ictp_i_c_i   = compute_stats(reg_ictp_i_c_i)
mean_reg_ictp_i_c_ii,  perc_reg_ictp_i_c_ii,  freq_reg_ictp_i_c_ii,  int_reg_ictp_i_c_ii  = compute_stats(reg_ictp_i_c_ii)
mean_reg_ictp_i_c_iii, perc_reg_ictp_i_c_iii, freq_reg_ictp_i_c_iii, int_reg_ictp_i_c_iii = compute_stats(reg_ictp_i_c_iii)
mean_reg_ictp_i_c_iv,  perc_reg_ictp_i_c_iv,  freq_reg_ictp_i_c_iv,  int_reg_ictp_i_c_iv  = compute_stats(reg_ictp_i_c_iv)
mean_reg_ictp_i_c_v,   perc_reg_ictp_i_c_v,   freq_reg_ictp_i_c_v,   int_reg_ictp_i_c_v   = compute_stats(reg_ictp_i_c_v)

mean_reg_ictp_ii_c_i,   perc_reg_ictp_ii_c_i,   freq_reg_ictp_ii_c_i,   int_reg_ictp_ii_c_i   = compute_stats(reg_ictp_ii_c_i)
mean_reg_ictp_ii_c_ii,  perc_reg_ictp_ii_c_ii,  freq_reg_ictp_ii_c_ii,  int_reg_ictp_ii_c_ii  = compute_stats(reg_ictp_ii_c_ii)
mean_reg_ictp_ii_c_iii, perc_reg_ictp_ii_c_iii, freq_reg_ictp_ii_c_iii, int_reg_ictp_ii_c_iii = compute_stats(reg_ictp_ii_c_iii)
mean_reg_ictp_ii_c_iv,  perc_reg_ictp_ii_c_iv,  freq_reg_ictp_ii_c_iv,  int_reg_ictp_ii_c_iv  = compute_stats(reg_ictp_ii_c_iv)
mean_reg_ictp_ii_c_v,   perc_reg_ictp_ii_c_v,   freq_reg_ictp_ii_c_v,   int_reg_ictp_ii_c_v   = compute_stats(reg_ictp_ii_c_v)

mean_wrf_ncar_c_i,   perc_wrf_ncar_c_i,   freq_wrf_ncar_c_i,   int_wrf_ncar_c_i   = compute_stats(wrf_ncar_c_i)
mean_wrf_ncar_c_ii,  perc_wrf_ncar_c_ii,  freq_wrf_ncar_c_ii,  int_wrf_ncar_c_ii  = compute_stats(wrf_ncar_c_ii)
mean_wrf_ncar_c_iii, perc_wrf_ncar_c_iii, freq_wrf_ncar_c_iii, int_wrf_ncar_c_iii = compute_stats(wrf_ncar_c_iii)
mean_wrf_ncar_c_iv,  perc_wrf_ncar_c_iv,  freq_wrf_ncar_c_iv,  int_wrf_ncar_c_iv  = compute_stats(wrf_ncar_c_iv)
mean_wrf_ncar_c_v,   perc_wrf_ncar_c_v,   freq_wrf_ncar_c_v,   int_wrf_ncar_c_v   = compute_stats(wrf_ncar_c_v)

mean_wrf_ucan_c_i,   perc_wrf_ucan_c_i,   freq_wrf_ucan_c_i,   int_wrf_ucan_c_i   = compute_stats(wrf_ucan_c_i)
mean_wrf_ucan_c_ii,  perc_wrf_ucan_c_ii,  freq_wrf_ucan_c_ii,  int_wrf_ucan_c_ii  = compute_stats(wrf_ucan_c_ii)
mean_wrf_ucan_c_iii, perc_wrf_ucan_c_iii, freq_wrf_ucan_c_iii, int_wrf_ucan_c_iii = compute_stats(wrf_ucan_c_iii)
mean_wrf_ucan_c_iv,  perc_wrf_ucan_c_iv,  freq_wrf_ucan_c_iv,  int_wrf_ucan_c_iv  = compute_stats(wrf_ucan_c_iv)
mean_wrf_ucan_c_v,   perc_wrf_ucan_c_v,   freq_wrf_ucan_c_v,   int_wrf_ucan_c_v   = compute_stats(wrf_ucan_c_v)

print(int_inmet_smn_c_i)
print()
print(int_era5_c_i)
print()
print(int_reg_usp_c_i)
print()
print(int_reg_ictp_c_i)
print()
print(int_reg_ictp_i_c_i)
print()
print(int_reg_ictp_ii_c_i)
print()
print(int_wrf_ncar_c_i)
print()
print(int_wrf_ucan_c_i)

data_mean_c_i   = np.vstack([mean_inmet_smn_c_i,   mean_era5_c_i,   mean_reg_usp_c_i,   mean_reg_ictp_c_i,   mean_reg_ictp_i_c_i,   mean_reg_ictp_ii_c_i,   mean_wrf_ncar_c_i,   mean_wrf_ucan_c_i])
data_mean_c_ii  = np.vstack([mean_inmet_smn_c_ii,  mean_era5_c_ii,  mean_reg_usp_c_ii,  mean_reg_ictp_c_ii,  mean_reg_ictp_i_c_ii,  mean_reg_ictp_ii_c_ii,  mean_wrf_ncar_c_ii,  mean_wrf_ucan_c_ii])
data_mean_c_iii = np.vstack([mean_inmet_smn_c_iii, mean_era5_c_iii, mean_reg_usp_c_iii, mean_reg_ictp_c_iii, mean_reg_ictp_i_c_iii, mean_reg_ictp_ii_c_iii, mean_wrf_ncar_c_iii, mean_wrf_ucan_c_iii])
data_mean_c_iv  = np.vstack([mean_inmet_smn_c_iv,  mean_era5_c_iv,  mean_reg_usp_c_iv,  mean_reg_ictp_c_iv,  mean_reg_ictp_i_c_iv,  mean_reg_ictp_ii_c_iv,  mean_wrf_ncar_c_iv,  mean_wrf_ucan_c_iv])
data_mean_c_v   = np.vstack([mean_inmet_smn_c_v,   mean_era5_c_v,   mean_reg_usp_c_v,   mean_reg_ictp_c_v,   mean_reg_ictp_i_c_v,   mean_reg_ictp_ii_c_v,   mean_wrf_ncar_c_v,   mean_wrf_ucan_c_v])

data_perc_c_i   = np.vstack([perc_inmet_smn_c_i,   perc_era5_c_i,   perc_reg_usp_c_i,	perc_reg_ictp_c_i,   perc_reg_ictp_i_c_i,   perc_reg_ictp_ii_c_i,   perc_wrf_ncar_c_i,   perc_wrf_ucan_c_i])
data_perc_c_ii  = np.vstack([perc_inmet_smn_c_ii,  perc_era5_c_ii,  perc_reg_usp_c_ii,  perc_reg_ictp_c_ii,  perc_reg_ictp_i_c_ii,  perc_reg_ictp_ii_c_ii,  perc_wrf_ncar_c_ii,  perc_wrf_ucan_c_ii])
data_perc_c_iii = np.vstack([perc_inmet_smn_c_iii, perc_era5_c_iii, perc_reg_usp_c_iii, perc_reg_ictp_c_iii, perc_reg_ictp_i_c_iii, perc_reg_ictp_ii_c_iii, perc_wrf_ncar_c_iii, perc_wrf_ucan_c_iii])
data_perc_c_iv  = np.vstack([perc_inmet_smn_c_iv,  perc_era5_c_iv,  perc_reg_usp_c_iv,  perc_reg_ictp_c_iv,  perc_reg_ictp_i_c_iv,  perc_reg_ictp_ii_c_iv,  perc_wrf_ncar_c_iv,  perc_wrf_ucan_c_iv])
data_perc_c_v   = np.vstack([perc_inmet_smn_c_v,   perc_era5_c_v,   perc_reg_usp_c_v,	perc_reg_ictp_c_v,   perc_reg_ictp_i_c_v,   perc_reg_ictp_ii_c_v,   perc_wrf_ncar_c_v,   perc_wrf_ucan_c_v])

data_freq_c_i   = np.vstack([freq_inmet_smn_c_i,   freq_era5_c_i,   freq_reg_usp_c_i,	freq_reg_ictp_c_i,   freq_reg_ictp_i_c_i,   freq_reg_ictp_ii_c_i,   freq_wrf_ncar_c_i,   freq_wrf_ucan_c_i])
data_freq_c_ii  = np.vstack([freq_inmet_smn_c_ii,  freq_era5_c_ii,  freq_reg_usp_c_ii,  freq_reg_ictp_c_ii,  freq_reg_ictp_i_c_ii,  freq_reg_ictp_ii_c_ii,  freq_wrf_ncar_c_ii,  freq_wrf_ucan_c_ii])
data_freq_c_iii = np.vstack([freq_inmet_smn_c_iii, freq_era5_c_iii, freq_reg_usp_c_iii, freq_reg_ictp_c_iii, freq_reg_ictp_i_c_iii, freq_reg_ictp_ii_c_iii, freq_wrf_ncar_c_iii, freq_wrf_ucan_c_iii])
data_freq_c_iv  = np.vstack([freq_inmet_smn_c_iv,  freq_era5_c_iv,  freq_reg_usp_c_iv,  freq_reg_ictp_c_iv,  freq_reg_ictp_i_c_iv,  freq_reg_ictp_ii_c_iv,  freq_wrf_ncar_c_iv,  freq_wrf_ucan_c_iv])
data_freq_c_v   = np.vstack([freq_inmet_smn_c_v,   freq_era5_c_v,   freq_reg_usp_c_v,	freq_reg_ictp_c_v,   freq_reg_ictp_i_c_v,   freq_reg_ictp_ii_c_v,   freq_wrf_ncar_c_v,   freq_wrf_ucan_c_v])

data_int_c_i   = np.vstack([int_inmet_smn_c_i,   int_era5_c_i,   int_reg_usp_c_i,   int_reg_ictp_c_i,   int_reg_ictp_i_c_i,   int_reg_ictp_ii_c_i,   int_wrf_ncar_c_i,   int_wrf_ucan_c_i])
data_int_c_ii  = np.vstack([int_inmet_smn_c_ii,  int_era5_c_ii,  int_reg_usp_c_ii,  int_reg_ictp_c_ii,  int_reg_ictp_i_c_ii,  int_reg_ictp_ii_c_ii,  int_wrf_ncar_c_ii,  int_wrf_ucan_c_ii])
data_int_c_iii = np.vstack([int_inmet_smn_c_iii, int_era5_c_iii, int_reg_usp_c_iii, int_reg_ictp_c_iii, int_reg_ictp_i_c_iii, int_reg_ictp_ii_c_iii, int_wrf_ncar_c_iii, int_wrf_ucan_c_iii])
data_int_c_iv  = np.vstack([int_inmet_smn_c_iv,  int_era5_c_iv,  int_reg_usp_c_iv,  int_reg_ictp_c_iv,  int_reg_ictp_i_c_iv,  int_reg_ictp_ii_c_iv,  int_wrf_ncar_c_iv,  int_wrf_ucan_c_iv])
data_int_c_v   = np.vstack([int_inmet_smn_c_v,   int_era5_c_v,   int_reg_usp_c_v,   int_reg_ictp_c_v,	int_reg_ictp_i_c_v,   int_reg_ictp_ii_c_v,   int_wrf_ncar_c_v,   int_wrf_ucan_c_v])

# Plot figure
fig = plt.figure(figsize=(12, 10))
fig.subplots_adjust(wspace=0.04, hspace=0.04)

font_size = 8

labels = ['WS', 'ERA5', 'Reg4', 'Reg5-Holt3', 'Reg5-Holt', 'Reg5-UW', 'WRF-NCAR', 'WRF-UCAN']
	
legend_i = 'mm h⁻¹'
cmap_i = 'terrain_r'
mvmin = 0
mvmax = 0.8
pvmin = 0
pvmax = 10
fvmin = 0
fvmax = 6
ivmin = 0
ivmax = 15

ax = fig.add_subplot(5, 4, 1)
im = ax.imshow(data_mean_c_i, aspect='1.5', origin='lower', cmap=cmap_i, vmin=mvmin, vmax=mvmax)
ax.set_title('(a)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(labels, fontsize=font_size)
ax.set_title(f'MEAN ({legend_i})', fontsize=font_size)
ax.set_ylabel('CLUSTER I', fontsize=font_size)
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)
cax = fig.add_axes([0.12, 0.08, 0.185, 0.015])  
cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(5, 4, 2)
im = ax.imshow(data_perc_c_i, aspect='1.5', origin='lower', cmap=cmap_i, vmin=pvmin, vmax=pvmax)
ax.set_title('(b)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_title(f'P99 ({legend_i})', fontsize=font_size)
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)
cax = fig.add_axes([0.32, 0.08, 0.185, 0.015])  
cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(5, 4, 3)
im = ax.imshow(data_freq_c_i, aspect='1.5', origin='lower', cmap=cmap_i, vmin=fvmin, vmax=fvmax)
ax.set_title('(c)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_title('FREQUENCY (%)', fontsize=font_size)
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)
cax = fig.add_axes([0.515, 0.08, 0.185, 0.015])  
cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(5, 4, 4)
im = ax.imshow(data_int_c_i, aspect='1.5', origin='lower', cmap=cmap_i, vmin=ivmin, vmax=ivmax)
ax.set_title('(d)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_title(f'INTENSITY ({legend_i})', fontsize=font_size)
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)
cax = fig.add_axes([0.71, 0.08, 0.185, 0.015])  
cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
cbar.ax.tick_params(labelsize=font_size)

ax = fig.add_subplot(5, 4, 5)
im = ax.imshow(data_mean_c_ii, aspect='1.5', origin='lower', cmap=cmap_i, vmin=mvmin, vmax=mvmax)
ax.set_title('(e)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(labels, fontsize=font_size)
ax.set_ylabel('CLUSTER II', fontsize=font_size)
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 6)
im = ax.imshow(data_perc_c_ii, aspect='1.5', origin='lower', cmap=cmap_i, vmin=pvmin, vmax=pvmax)
ax.set_title('(f)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 7)
im = ax.imshow(data_freq_c_ii, aspect='1.5', origin='lower', cmap=cmap_i, vmin=fvmin, vmax=fvmax)
ax.set_title('(g)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 8)
im = ax.imshow(data_int_c_ii, aspect='1.5', origin='lower', cmap=cmap_i, vmin=ivmin, vmax=ivmax)
ax.set_title('(h)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 9)
im = ax.imshow(data_mean_c_iii, aspect='1.5', origin='lower', cmap=cmap_i, vmin=mvmin, vmax=mvmax)
ax.set_title('(i)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(labels, fontsize=font_size)
ax.set_ylabel('CLUSTER III', fontsize=font_size)
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 10)
im = ax.imshow(data_perc_c_iii, aspect='1.5', origin='lower', cmap=cmap_i, vmin=pvmin, vmax=pvmax)
ax.set_title('(j)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 11)
im = ax.imshow(data_freq_c_iii, aspect='1.5', origin='lower', cmap=cmap_i, vmin=fvmin, vmax=fvmax)
ax.set_title('(k)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 12)
im = ax.imshow(data_int_c_iii, aspect='1.5', origin='lower', cmap=cmap_i, vmin=ivmin, vmax=ivmax)
ax.set_title('(l)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 13)
im = ax.imshow(data_mean_c_iv, aspect='1.5', origin='lower', cmap=cmap_i, vmin=mvmin, vmax=mvmax)
ax.set_title('(m)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(labels, fontsize=font_size)
ax.set_ylabel('CLUSTER IV', fontsize=font_size)
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 14)
im = ax.imshow(data_perc_c_iv, aspect='1.5', origin='lower', cmap=cmap_i, vmin=pvmin, vmax=pvmax)
ax.set_title('(n)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 15)
im = ax.imshow(data_freq_c_iv, aspect='1.5', origin='lower', cmap=cmap_i, vmin=fvmin, vmax=fvmax)
ax.set_title('(o)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 16)
im = ax.imshow(data_int_c_iv, aspect='1.5', origin='lower', cmap=cmap_i, vmin=ivmin, vmax=ivmax)
ax.set_title('(p)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 17)
im = ax.imshow(data_mean_c_v, aspect='1.5', origin='lower', cmap=cmap_i, vmin=mvmin, vmax=mvmax)
ax.set_title('(q)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels(labels, fontsize=font_size)
ax.set_ylabel('CLUSTER V', fontsize=font_size)
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 18)
im = ax.imshow(data_perc_c_v, aspect='1.5', origin='lower', cmap=cmap_i, vmin=pvmin, vmax=pvmax)
ax.set_title('(r)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 19)
im = ax.imshow(data_freq_c_v, aspect='1.5', origin='lower', cmap=cmap_i, vmin=fvmin, vmax=fvmax)
ax.set_title('(s)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

ax = fig.add_subplot(5, 4, 20)
im = ax.imshow(data_int_c_v, aspect='1.5', origin='lower', cmap=cmap_i, vmin=ivmin, vmax=ivmax)
ax.set_title('(t)', loc='left', fontsize=font_size, fontweight='bold')
ax.set_xticks(np.arange(24))
ax.set_xticklabels(('00', '', '02', '', '04', '', '06', '', '08', '', '10', '', '12', '', '14', '', '16', '', '18', '', '20', '', '22', ''), fontsize=font_size)
ax.set_yticks(np.arange(len(labels)))
ax.set_yticklabels([])
ax.set_xticks(np.arange(-.5, 24, 1), minor=True)
ax.set_yticks(np.arange(-.5, len(labels), 1), minor=True)
ax.grid(which='minor', linestyle='--', alpha=0.25)
ax.tick_params(which='minor', bottom=False, left=False)

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_graph_diurnal_cycle_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
exit()
