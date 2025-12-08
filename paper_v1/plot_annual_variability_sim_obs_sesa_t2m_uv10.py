# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Sept 22, 2025"
__description__ = "This script plot annual variability"

import os
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

var = 'uv10' # t2m or uv10
path = '/home/mda_silv/users/FPS_SESA'

def import_inmet():
	
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []
	mean_v = []
	mean_vi = []
	mean_vii = []

	# Select lat and lon 
	for i in range(1, 99):
		yy=inmet[i][2]
		xx=inmet[i][3]
		
		print('Reading weather station:', i, inmet[i][0])
		if var == 't2m':		
			# reading regcm usp 
			d_i = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.values
			mean_i.append(d_i-273.15)
			
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.values
			mean_ii.append(d_ii-273.15)
		
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.values
			mean_iii.append(d_iii-273.15)
						
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'tas_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.values
			mean_iv.append(d_iv-273.15)
			
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'tas_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.tas.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.values
			mean_v.append(d_v-273.15)
		
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/tmp/'.format(path) + 'tmp_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.tmp.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.resample(time='1ME').mean()
			d_vi = d_vi.values
			mean_vi.append(d_vi)
		
			# reading era5 
			d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 't2m_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.t2m.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.values
			mean_vii.append(d_vii-273.15)
		
		else:
			# reading regcm usp 
			d_i = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
			d_i = d_i.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_i = d_i.values
			mean_i.append(d_i)
			
			# reading regcm ictp pbl 1 
			d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
			d_ii = d_ii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_ii = d_ii.values
			mean_ii.append(d_ii)
		
			# reading regcm ictp pbl 2
			d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
			d_iii = d_iii.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iii = d_iii.values
			mean_iii.append(d_iii)
						
			# reading wrf ncar 
			d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'sfcWind_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
			d_iv = d_iv.uas.sel(time=slice('2018-06-01','2021-05-31'))
			d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_iv = d_iv.values
			mean_iv.append(d_iv)
			
			# reading wrf ucan 
			d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'sfcWind_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
			d_v = d_v.sfcWind.sel(time=slice('2018-06-01','2021-05-31'))
			d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_v = d_v.values
			mean_v.append(d_v)
		
			# Reading inmet 
			d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/uv/'.format(path) + 'uv_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
			d_vi = d_vi.uv.sel(time=slice('2018-06-01','2021-05-31'))
			d_vi = d_vi.resample(time='1ME').mean()
			d_vi = d_vi.values
			mean_vi.append(d_vi)
		
			# reading era5 
			d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'uv10_era5_csam_4km_mon_20180101-20211231.nc')
			d_vii = d_vii.u10.sel(time=slice('2018-06-01','2021-05-31'))
			d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
			d_vii = d_vii.values
			mean_vii.append(d_vii)
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


# Import dataset
clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x, clim_vi_x, clim_vii_x = import_inmet()			

reg_usp     = clim_i_x  
reg_ictp_i  = clim_ii_x  
reg_ictp_ii = clim_iii_x 
wrf_ncar    = clim_iv_x  
wrf_ucan    = clim_v_x   
inmet_smn   = clim_vi_x  
era5        = clim_vii_x 

list_hc = [4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 0, 4, 0, 4, 0, 0, 0, 0, 4, 0, 4, 0, 4, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 2, 0, 2, 3, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 4, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 1, 4, 4, 3, 2, 3, 3, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 3, 3, 3, 3, 2, 1, 2, 3, 1, 2, 0, 0, 3, 0, 1, 1, 0, 2, 1, 3, 2, 0, 2, 2, 1, 1, 2, 1, 1, 2, 1, 2, 0, 1, 2, 1, 1, 0, 4, 3, 3, 4, 3, 3, 4, 0, 3, 0, 4, 3, 4, 3, 3, 3]
list_hc = list_hc[:98]

print(len(reg_usp))
print(len(reg_ictp_i))
print(len(reg_ictp_ii))
print(len(wrf_ncar))
print(len(wrf_ucan))
print(len(inmet_smn))
print(len(era5))
print(len(list_hc))

count_i, count_iii, count_v = [], [], []
for count, idx in enumerate(list_hc):
	if idx == 0:
		count_i.append(count)
	if idx == 2:
		count_iii.append(count)
	if idx == 4:
		count_v.append(count)

reg_usp_i,     reg_usp_iii,     reg_usp_v     = [], [], []
reg_ictp_i_i,  reg_ictp_i_iii,  reg_ictp_i_v  = [], [], []
reg_ictp_ii_i, reg_ictp_ii_iii, reg_ictp_ii_v = [], [], []
wrf_ncar_i,    wrf_ncar_iii,    wrf_ncar_v    = [], [], []
wrf_ucan_i,    wrf_ucan_iii,    wrf_ucan_v    = [], [], []
inmet_smn_i,   inmet_smn_iii,   inmet_smn_v   = [], [], []
era5_i,        era5_iii,        era5_v        = [], [], []

for c_i in count_i:
	reg_usp_i.append(reg_usp[c_i])
	reg_ictp_i_i.append(reg_ictp_i[c_i])
	reg_ictp_ii_i.append(reg_ictp_ii[c_i])
	wrf_ncar_i.append(wrf_ncar[c_i])
	wrf_ucan_i.append(wrf_ucan[c_i])
	inmet_smn_i.append(inmet_smn[c_i])
	era5_i.append(era5[c_i])
	
for c_iii in count_iii:
	reg_usp_iii.append(reg_usp[c_iii])
	reg_ictp_i_iii.append(reg_ictp_i[c_iii])
	reg_ictp_ii_iii.append(reg_ictp_ii[c_iii])
	wrf_ncar_iii.append(wrf_ncar[c_iii])
	wrf_ucan_iii.append(wrf_ucan[c_iii])
	inmet_smn_iii.append(inmet_smn[c_iii])
	era5_iii.append(era5[c_iii])
	
for c_v in count_v:
	reg_usp_v.append(reg_usp[c_v])
	reg_ictp_i_v.append(reg_ictp_i[c_v])
	reg_ictp_ii_v.append(reg_ictp_ii[c_v])
	wrf_ncar_v.append(wrf_ncar[c_v])
	wrf_ucan_v.append(wrf_ucan[c_v])
	inmet_smn_v.append(inmet_smn[c_v])
	era5_v.append(era5[c_v])

# Group I
# Average
reg_usp_c_i     = np.nanmean(reg_usp_i, axis=0)
reg_ictp_i_c_i  = np.nanmean(reg_ictp_i_i, axis=0)
reg_ictp_ii_c_i = np.nanmean(reg_ictp_ii_i, axis=0)
wrf_ncar_c_i    = np.nanmean(wrf_ncar_i, axis=0)
wrf_ucan_c_i    = np.nanmean(wrf_ucan_i, axis=0)
inmet_smn_c_i   = np.nanmean(inmet_smn_i, axis=0)
era5_c_i        = np.nanmean(era5_i, axis=0)

# Bias
b_reg_usp_inmet_ci     = np.nanmean(reg_usp_c_i - inmet_smn_c_i)
b_reg_ictp_i_inmet_ci  = np.nanmean(reg_ictp_i_c_i - inmet_smn_c_i)
b_reg_ictp_ii_inmet_ci = np.nanmean(reg_ictp_ii_c_i - inmet_smn_c_i)
b_wrf_ncar_inmet_ci    = np.nanmean(wrf_ncar_c_i - inmet_smn_c_i)
b_wrf_ucan_inmet_ci    = np.nanmean(wrf_ucan_c_i - inmet_smn_c_i)

b_reg_usp_era5_ci     = np.nanmean(reg_usp_c_i - era5_c_i)
b_reg_ictp_i_era5_ci  = np.nanmean(reg_ictp_i_c_i - era5_c_i)
b_reg_ictp_ii_era5_ci = np.nanmean(reg_ictp_ii_c_i - era5_c_i)
b_wrf_ncar_era5_ci    = np.nanmean(wrf_ncar_c_i - era5_c_i)
b_wrf_ucan_era5_ci    = np.nanmean(wrf_ucan_c_i - era5_c_i)

# Correlation
r_reg_usp_inmet_ci     = np.corrcoef(reg_usp_c_i, inmet_smn_c_i)[0, 1]
r_reg_ictp_i_inmet_ci  = np.corrcoef(reg_ictp_i_c_i, inmet_smn_c_i)[0, 1]
r_reg_ictp_ii_inmet_ci = np.corrcoef(reg_ictp_ii_c_i, inmet_smn_c_i)[0, 1]
r_wrf_ncar_inmet_ci    = np.corrcoef(wrf_ncar_c_i, inmet_smn_c_i)[0, 1]
r_wrf_ucan_inmet_ci    = np.corrcoef(wrf_ucan_c_i, inmet_smn_c_i)[0, 1]

r_reg_usp_era5_ci     = np.corrcoef(reg_usp_c_i, era5_c_i)[0, 1]
r_reg_ictp_i_era5_ci   = np.corrcoef(reg_ictp_i_c_i, era5_c_i)[0, 1]
r_reg_ictp_ii_era5_ci = np.corrcoef(reg_ictp_ii_c_i, era5_c_i)[0, 1]
r_wrf_ncar_era5_ci    = np.corrcoef(wrf_ncar_c_i, era5_c_i)[0, 1]
r_wrf_ucan_era5_ci    = np.corrcoef(wrf_ucan_c_i, era5_c_i)[0, 1]

# Group III
# Average
reg_usp_c_iii     = np.nanmean(reg_usp_iii, axis=0)
reg_ictp_i_c_iii  = np.nanmean(reg_ictp_i_iii, axis=0)
reg_ictp_ii_c_iii = np.nanmean(reg_ictp_ii_iii, axis=0)
wrf_ncar_c_iii    = np.nanmean(wrf_ncar_iii, axis=0)
wrf_ucan_c_iii    = np.nanmean(wrf_ucan_iii, axis=0)
inmet_smn_c_iii   = np.nanmean(inmet_smn_iii, axis=0)
era5_c_iii        = np.nanmean(era5_iii, axis=0)

# Bias
b_reg_usp_inmet_ciii     = np.nanmean(reg_usp_c_iii - inmet_smn_c_iii)
b_reg_ictp_i_inmet_ciii  = np.nanmean(reg_ictp_i_c_iii - inmet_smn_c_iii)
b_reg_ictp_ii_inmet_ciii = np.nanmean(reg_ictp_ii_c_iii - inmet_smn_c_iii)
b_wrf_ncar_inmet_ciii    = np.nanmean(wrf_ncar_c_iii - inmet_smn_c_iii)
b_wrf_ucan_inmet_ciii    = np.nanmean(wrf_ucan_c_iii - inmet_smn_c_iii)

b_reg_usp_era5_ciii     = np.nanmean(reg_usp_c_iii - era5_c_iii)
b_reg_ictp_i_era5_ciii  = np.nanmean(reg_ictp_i_c_iii - era5_c_iii)
b_reg_ictp_ii_era5_ciii = np.nanmean(reg_ictp_ii_c_iii - era5_c_iii)
b_wrf_ncar_era5_ciii    = np.nanmean(wrf_ncar_c_iii - era5_c_iii)
b_wrf_ucan_era5_ciii    = np.nanmean(wrf_ucan_c_iii - era5_c_iii)

# Correlation
r_reg_usp_inmet_ciii     = np.corrcoef(reg_usp_c_iii, inmet_smn_c_iii)[0, 1]
r_reg_ictp_i_inmet_ciii  = np.corrcoef(reg_ictp_i_c_iii, inmet_smn_c_iii)[0, 1]
r_reg_ictp_ii_inmet_ciii = np.corrcoef(reg_ictp_ii_c_iii, inmet_smn_c_iii)[0, 1]
r_wrf_ncar_inmet_ciii  = np.corrcoef(wrf_ncar_c_iii, inmet_smn_c_iii)[0, 1]
r_wrf_ucan_inmet_ciii  = np.corrcoef(wrf_ucan_c_iii, inmet_smn_c_iii)[0, 1]

r_reg_usp_era5_ciii     = np.corrcoef(reg_usp_c_iii, era5_c_iii)[0, 1]
r_reg_ictp_i_era5_ciii  = np.corrcoef(reg_ictp_i_c_iii, era5_c_iii)[0, 1]
r_reg_ictp_ii_era5_ciii = np.corrcoef(reg_ictp_ii_c_iii, era5_c_iii)[0, 1]
r_wrf_ncar_era5_ciii    = np.corrcoef(wrf_ncar_c_iii, era5_c_iii)[0, 1]
r_wrf_ucan_era5_ciii    = np.corrcoef(wrf_ucan_c_iii, era5_c_iii)[0, 1]

# Group V
# Average
reg_usp_c_v     = np.nanmean(reg_usp_v, axis=0)
reg_ictp_i_c_v  = np.nanmean(reg_ictp_i_v, axis=0)
reg_ictp_ii_c_v = np.nanmean(reg_ictp_ii_v, axis=0)
wrf_ncar_c_v    = np.nanmean(wrf_ncar_v, axis=0)
wrf_ucan_c_v    = np.nanmean(wrf_ucan_v, axis=0)
inmet_smn_c_v   = np.nanmean(inmet_smn_v, axis=0)
era5_c_v        = np.nanmean(era5_v, axis=0)

# Bias
b_reg_usp_inmet_cv     = np.nanmean(reg_usp_c_v - inmet_smn_c_v)
b_reg_ictp_i_inmet_cv  = np.nanmean(reg_ictp_i_c_v - inmet_smn_c_v)
b_reg_ictp_ii_inmet_cv = np.nanmean(reg_ictp_ii_c_v - inmet_smn_c_v)
b_wrf_ncar_inmet_cv    = np.nanmean(wrf_ncar_c_v - inmet_smn_c_v)
b_wrf_ucan_inmet_cv    = np.nanmean(wrf_ucan_c_v - inmet_smn_c_v)

b_reg_usp_era5_cv     = np.nanmean(reg_usp_c_v - era5_c_v)
b_reg_ictp_i_era5_cv  = np.nanmean(reg_ictp_i_c_v - era5_c_v)
b_reg_ictp_ii_era5_cv = np.nanmean(reg_ictp_ii_c_v - era5_c_v)
b_wrf_ncar_era5_cv    = np.nanmean(wrf_ncar_c_v - era5_c_v)
b_wrf_ucan_era5_cv    = np.nanmean(wrf_ucan_c_v - era5_c_v)

# Correlation
r_reg_usp_inmet_cv     = np.corrcoef(reg_usp_c_v, inmet_smn_c_v)[0, 1]
r_reg_ictp_i_inmet_cv  = np.corrcoef(reg_ictp_i_c_v, inmet_smn_c_v)[0, 1]
r_reg_ictp_ii_inmet_cv = np.corrcoef(reg_ictp_ii_c_v, inmet_smn_c_v)[0, 1]
r_wrf_ncar_inmet_cv    = np.corrcoef(wrf_ncar_c_v, inmet_smn_c_v)[0, 1]
r_wrf_ucan_inmet_cv    = np.corrcoef(wrf_ucan_c_v, inmet_smn_c_v)[0, 1]

r_reg_usp_era5_cv     = np.corrcoef(reg_usp_c_v, era5_c_v)[0, 1]
r_reg_ictp_i_era5_cv   = np.corrcoef(reg_ictp_i_c_v, era5_c_v)[0, 1]
r_reg_ictp_ii_era5_cv = np.corrcoef(reg_ictp_ii_c_v, era5_c_v)[0, 1]
r_wrf_ncar_era5_cv    = np.corrcoef(wrf_ncar_c_v, era5_c_v)[0, 1]
r_wrf_ucan_era5_cv    = np.corrcoef(wrf_ucan_c_v, era5_c_v)[0, 1]

# All groups 
# Average
reg_usp_c     = np.nanmean(reg_usp, axis=0)
reg_ictp_i_c  = np.nanmean(reg_ictp_i, axis=0)
reg_ictp_ii_c = np.nanmean(reg_ictp_ii, axis=0)
wrf_ncar_c    = np.nanmean(wrf_ncar, axis=0)
wrf_ucan_c    = np.nanmean(wrf_ucan, axis=0)
inmet_smn_c   = np.nanmean(inmet_smn, axis=0)
era5_c        = np.nanmean(era5, axis=0)

# Bias
b_reg_usp_inmet_c     = np.nanmean(reg_usp_c - inmet_smn_c)
b_reg_ictp_i_inmet_c  = np.nanmean(reg_ictp_i_c - inmet_smn_c)
b_reg_ictp_ii_inmet_c = np.nanmean(reg_ictp_ii_c - inmet_smn_c)
b_wrf_ncar_inmet_c    = np.nanmean(wrf_ncar_c - inmet_smn_c)
b_wrf_ucan_inmet_c    = np.nanmean(wrf_ucan_c - inmet_smn_c)

b_reg_usp_era5_c     = np.nanmean(reg_usp_c - era5_c)
b_reg_ictp_i_era5_c  = np.nanmean(reg_ictp_i_c - era5_c)
b_reg_ictp_ii_era5_c = np.nanmean(reg_ictp_ii_c - era5_c)
b_wrf_ncar_era5_c    = np.nanmean(wrf_ncar_c - era5_c)
b_wrf_ucan_era5_c    = np.nanmean(wrf_ucan_c - era5_c)

# Correlation
r_reg_usp_inmet_c     = np.corrcoef(reg_usp_c, inmet_smn_c)[0, 1]
r_reg_ictp_i_inmet_c  = np.corrcoef(reg_ictp_i_c, inmet_smn_c)[0, 1]
r_reg_ictp_ii_inmet_c = np.corrcoef(reg_ictp_ii_c, inmet_smn_c)[0, 1]
r_wrf_ncar_inmet_c    = np.corrcoef(wrf_ncar_c, inmet_smn_c)[0, 1]
r_wrf_ucan_inmet_c    = np.corrcoef(wrf_ucan_c, inmet_smn_c)[0, 1]

r_reg_usp_era5_c     = np.corrcoef(reg_usp_c, era5_c)[0, 1]
r_reg_ictp_i_era5_c   = np.corrcoef(reg_ictp_i_c, era5_c)[0, 1]
r_reg_ictp_ii_era5_c = np.corrcoef(reg_ictp_ii_c, era5_c)[0, 1]
r_wrf_ncar_era5_c    = np.corrcoef(wrf_ncar_c, era5_c)[0, 1]
r_wrf_ucan_era5_c    = np.corrcoef(wrf_ucan_c, era5_c)[0, 1]

# Plot figure
fig = plt.figure(figsize=(12, 6))
font_size = 8

if var == 't2m':
	legend = 'Temperature (°C)'
	vmin = 12
	vmax = 30
	vmax_ = 32
	int_ = 2
	tx1, tx2, tx3, tx4, tx5 = 0.95, 0.85, 0.75, 0.65, 0.55
	ty1 = 0.1
	ty2 = 0.5
else:
	legend = 'Wind 10m (m s⁻¹)'
	vmin = 1
	vmax = 6
	vmax_ = 6.5
	int_ = 0.5
	tx1, tx2, tx3, tx4, tx5 = 0.95, 0.85, 0.75, 0.65, 0.55
	ty1 = 0.1
	ty2 = 0.5

dt = pd.date_range(start="20180601", end="20210531", freq="ME")

ax = fig.add_subplot(2, 2, 1)
reg_usp_c_i_dt     = pd.Series(data=reg_usp_c_i, index=dt)
reg_ictp_i_c_i_dt  = pd.Series(data=reg_ictp_i_c_i, index=dt)
reg_ictp_ii_c_i_dt = pd.Series(data=reg_ictp_i_c_i, index=dt)
wrf_ncar_c_i_dt    = pd.Series(data=wrf_ncar_c_i, index=dt)
wrf_ucan_c_i_dt    = pd.Series(data=wrf_ucan_c_i, index=dt)
inmet_smn_c_i_dt   = pd.Series(data=inmet_smn_c_i, index=dt)
era5_c_i_dt        = pd.Series(data=era5_c_i, index=dt)
plt.plot(inmet_smn_c_i_dt,   linewidth=1, color='black', label = 'INMET')
plt.plot(era5_c_i_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_i_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_i_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_i_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_i_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_i_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(ty1, tx1, f"Reg4 = {b_reg_usp_inmet_ci:.2f}({r_reg_usp_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx2, f"Reg5-Holt = {b_reg_ictp_i_inmet_ci:.2f}({r_reg_ictp_i_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx3, f"Reg5-UW = {b_reg_ictp_ii_inmet_ci:.2f}({r_reg_ictp_ii_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx4, f"WRF-NCAR = {b_wrf_ncar_inmet_ci:.2f}({r_wrf_ncar_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx5, f"WRF-UCAN = {b_wrf_ucan_inmet_ci:.2f}({r_wrf_ucan_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty2, tx1, f"Reg4 = {b_reg_usp_era5_ci:.2f}({r_reg_usp_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx2, f"Reg5-Holt = {b_reg_ictp_i_era5_ci:.2f}({r_reg_ictp_i_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx3, f"Reg5-UW = {b_reg_ictp_ii_era5_ci:.2f}({r_reg_ictp_ii_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx4, f"WRF-NCAR = {b_wrf_ncar_era5_ci:.2f}({r_wrf_ncar_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx5, f"WRF-UCAN = {b_wrf_ucan_era5_ci:.2f}({r_wrf_ucan_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(2, 2, 2)
reg_usp_c_iii_dt     = pd.Series(data=reg_usp_c_iii, index=dt)
reg_ictp_i_c_iii_dt  = pd.Series(data=reg_ictp_i_c_iii, index=dt)
reg_ictp_ii_c_iii_dt = pd.Series(data=reg_ictp_ii_c_iii, index=dt)
wrf_ncar_c_iii_dt    = pd.Series(data=wrf_ncar_c_iii, index=dt)
wrf_ucan_c_iii_dt    = pd.Series(data=wrf_ucan_c_iii, index=dt)
inmet_smn_c_iii_dt   = pd.Series(data=inmet_smn_c_iii, index=dt)
era5_c_iii_dt        = pd.Series(data=era5_c_iii, index=dt)
plt.plot(inmet_smn_c_iii_dt,   linewidth=1, color='black', label = 'INMET')
plt.plot(era5_c_iii_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_iii_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_iii_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_iii_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_iii_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_iii_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(ty1, tx1, f"Reg4 = {b_reg_usp_inmet_ci:.2f}({r_reg_usp_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx2, f"Reg5-Holt = {b_reg_ictp_i_inmet_ci:.2f}({r_reg_ictp_i_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx3, f"Reg5-UW = {b_reg_ictp_ii_inmet_ci:.2f}({r_reg_ictp_ii_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx4, f"WRF-NCAR = {b_wrf_ncar_inmet_ci:.2f}({r_wrf_ncar_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx5, f"WRF-UCAN = {b_wrf_ucan_inmet_ci:.2f}({r_wrf_ucan_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty2, tx1, f"Reg4 = {b_reg_usp_era5_ci:.2f}({r_reg_usp_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx2, f"Reg5-Holt = {b_reg_ictp_i_era5_ci:.2f}({r_reg_ictp_i_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx3, f"Reg5-UW = {b_reg_ictp_ii_era5_ci:.2f}({r_reg_ictp_ii_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx4, f"WRF-NCAR = {b_wrf_ncar_era5_ci:.2f}({r_wrf_ncar_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx5, f"WRF-UCAN = {b_wrf_ucan_era5_ci:.2f}({r_wrf_ucan_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(b) Cluster III', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(2, 2, 3)
reg_usp_c_v_dt     = pd.Series(data=reg_usp_c_v, index=dt)
reg_ictp_i_c_v_dt  = pd.Series(data=reg_ictp_i_c_v, index=dt)
reg_ictp_ii_c_v_dt = pd.Series(data=reg_ictp_ii_c_v, index=dt)
wrf_ncar_c_v_dt    = pd.Series(data=wrf_ncar_c_v, index=dt)
wrf_ucan_c_v_dt    = pd.Series(data=wrf_ucan_c_v, index=dt)
inmet_smn_c_v_dt   = pd.Series(data=inmet_smn_c_v, index=dt)
era5_c_v_dt        = pd.Series(data=era5_c_v, index=dt)
plt.plot(inmet_smn_c_v_dt,   linewidth=1, color='black', label = 'INMET')
plt.plot(era5_c_v_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_v_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_v_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_v_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_v_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_v_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(ty1, tx1, f"Reg4 = {b_reg_usp_inmet_ci:.2f}({r_reg_usp_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx2, f"Reg5-Holt = {b_reg_ictp_i_inmet_ci:.2f}({r_reg_ictp_i_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx3, f"Reg5-UW = {b_reg_ictp_ii_inmet_ci:.2f}({r_reg_ictp_ii_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx4, f"WRF-NCAR = {b_wrf_ncar_inmet_ci:.2f}({r_wrf_ncar_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx5, f"WRF-UCAN = {b_wrf_ucan_inmet_ci:.2f}({r_wrf_ucan_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty2, tx1, f"Reg4 = {b_reg_usp_era5_ci:.2f}({r_reg_usp_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx2, f"Reg5-Holt = {b_reg_ictp_i_era5_ci:.2f}({r_reg_ictp_i_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx3, f"Reg5-UW = {b_reg_ictp_ii_era5_ci:.2f}({r_reg_ictp_ii_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx4, f"WRF-NCAR = {b_wrf_ncar_era5_ci:.2f}({r_wrf_ncar_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx5, f"WRF-UCAN = {b_wrf_ucan_era5_ci:.2f}({r_wrf_ucan_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(c) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Period 2018/06 - 2021/05', fontsize=font_size, fontweight='bold')
plt.ylabel('{0}'.format(legend), fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(2, 2, 4)
reg_usp_c_v_dt     = pd.Series(data=reg_usp_c_v, index=dt)
reg_ictp_i_c_v_dt  = pd.Series(data=reg_ictp_i_c_v, index=dt)
reg_ictp_ii_c_v_dt = pd.Series(data=reg_ictp_ii_c_v, index=dt)
wrf_ncar_c_v_dt    = pd.Series(data=wrf_ncar_c_v, index=dt)
wrf_ucan_c_v_dt    = pd.Series(data=wrf_ucan_c_v, index=dt)
inmet_smn_c_v_dt   = pd.Series(data=inmet_smn_c_v, index=dt)
era5_c_v_dt        = pd.Series(data=era5_c_v, index=dt)
plt.plot(inmet_smn_c_v_dt,   linewidth=1, color='black', label = 'INMET')
plt.plot(era5_c_v_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_v_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_v_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_v_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_v_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_v_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(ty1, tx1, f"Reg4 = {b_reg_usp_inmet_ci:.2f}({r_reg_usp_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx2, f"Reg5-Holt = {b_reg_ictp_i_inmet_ci:.2f}({r_reg_ictp_i_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx3, f"Reg5-UW = {b_reg_ictp_ii_inmet_ci:.2f}({r_reg_ictp_ii_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx4, f"WRF-NCAR = {b_wrf_ncar_inmet_ci:.2f}({r_wrf_ncar_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty1, tx5, f"WRF-UCAN = {b_wrf_ucan_inmet_ci:.2f}({r_wrf_ucan_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(ty2, tx1, f"Reg4 = {b_reg_usp_era5_ci:.2f}({r_reg_usp_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx2, f"Reg5-Holt = {b_reg_ictp_i_era5_ci:.2f}({r_reg_ictp_i_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx3, f"Reg5-UW = {b_reg_ictp_ii_era5_ci:.2f}({r_reg_ictp_ii_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx4, f"WRF-NCAR = {b_wrf_ncar_era5_ci:.2f}({r_wrf_ncar_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(ty2, tx5, f"WRF-UCAN = {b_wrf_ucan_era5_ci:.2f}({r_wrf_ucan_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(d) All clusters', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Period 2018/06 - 2021/05', fontsize=font_size, fontweight='bold')
plt.ylim(vmin, vmax)
plt.yticks(np.arange(vmin, vmax_, int_), fontsize=font_size)
plt.xticks(fontsize=7)

plt.legend(ncol=7, fontsize=font_size, loc=(-1.15, -0.35))

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_annual_variability_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


