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

var = 'pr'
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
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.values
		mean_i.append(d_i*86400)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.values
		mean_v.append(d_v*86400)
		
		# Reading inmet 
		d_vi = xr.open_dataset('{0}/database/obs/inmet/inmet_br/inmet_nc/daily/pre/'.format(path) + 'pre_{0}_D_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.resample(time='1ME').mean()
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.values
		mean_vii.append(d_vii)
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def import_smn_i():
	
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []
	mean_v = []
	mean_vi = []
	mean_vii = []
	
	# Select lat and lon 
	for i in range(1, 72):
		yy=smn_i[i][1]
		xx=smn_i[i][2]
		
		print('Reading weather station:', i, smn_i[i][0])	
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.values
		mean_i.append(d_i*86400)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.values
		mean_v.append(d_v*86400)
						
		# Reading smn 
		d_vi = xr.open_dataset('{0}/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.resample(time='1ME').mean()
		d_vi = d_vi.values
		mean_vi.append(d_vi*24)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.values
		mean_vii.append(d_vii)
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii


def import_smn_ii():
	
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []
	mean_v = []
	mean_vi = []
	mean_vii = []
	
	# Select lat and lon 
	for i in range(1, 86):
		yy=smn_ii[i][1]
		xx=smn_ii[i][2]
		
		print('Reading weather station:', i, smn_ii[i][0])	
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.values
		mean_i.append(d_i*86400)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_mon_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_mon_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.values
		mean_v.append(d_v*86400)
						
		# Reading smn 
		d_vi = xr.open_dataset('{0}/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.resample(time='1ME').mean()
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.values
		mean_vii.append(d_vii)
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii
	
	
# Import dataset
clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x, clim_vi_x, clim_vii_x = import_inmet()			
clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y, clim_vi_y, clim_vii_y = import_smn_i()			
clim_i_z, clim_ii_z, clim_iii_z, clim_iv_z, clim_v_z, clim_vi_z, clim_vii_z = import_smn_ii()			

reg_usp     = clim_i_x   + clim_i_y   + clim_i_z
reg_ictp_i  = clim_ii_x  + clim_ii_y  + clim_ii_z
reg_ictp_ii = clim_iii_x + clim_iii_y + clim_iii_z
wrf_ncar    = clim_iv_x  + clim_iv_y  + clim_iv_z
wrf_ucan    = clim_v_x   + clim_v_y   + clim_v_z
inmet_smn   = clim_vi_x  + clim_vi_y  + clim_vi_z
era5        = clim_vii_x + clim_vii_y + clim_vii_z

list_hc = [4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 0, 4, 0, 4, 0, 0, 0, 0, 4, 0, 4, 0, 4, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 2, 0, 2, 3, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 1, 4, 3, 4, 4, 4, 3, 4, 4, 4, 3, 4, 4, 1, 4, 4, 3, 2, 3, 3, 3, 3, 0, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 1, 1, 2, 2, 1, 1, 1, 2, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 3, 3, 3, 3, 2, 1, 2, 3, 1, 2, 0, 0, 3, 0, 1, 1, 0, 2, 1, 3, 2, 0, 2, 2, 1, 1, 2, 1, 1, 2, 1, 2, 0, 1, 2, 1, 1, 0, 4, 3, 3, 4, 3, 3, 4, 0, 3, 0, 4, 3, 4, 3, 3, 3]
 
print(len(reg_usp))
print(len(reg_ictp_i))
print(len(reg_ictp_ii))
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
reg_ictp_i_i,  reg_ictp_i_ii,  reg_ictp_i_iii,  reg_ictp_i_iv,  reg_ictp_i_v  = [], [], [], [], []
reg_ictp_ii_i, reg_ictp_ii_ii, reg_ictp_ii_iii, reg_ictp_ii_iv, reg_ictp_ii_v = [], [], [], [], []
wrf_ncar_i,    wrf_ncar_ii,    wrf_ncar_iii,    wrf_ncar_iv,    wrf_ncar_v    = [], [], [], [], []
wrf_ucan_i,    wrf_ucan_ii,    wrf_ucan_iii,    wrf_ucan_iv,    wrf_ucan_v    = [], [], [], [], []
inmet_smn_i,   inmet_smn_ii,   inmet_smn_iii,   inmet_smn_iv,   inmet_smn_v   = [], [], [], [], []
era5_i,        era5_ii,        era5_iii,        era5_iv,        era5_v        = [], [], [], [], []

for c_i in count_i:
	reg_usp_i.append(reg_usp[c_i])
	reg_ictp_i_i.append(reg_ictp_i[c_i])
	reg_ictp_ii_i.append(reg_ictp_ii[c_i])
	wrf_ncar_i.append(wrf_ncar[c_i])
	wrf_ucan_i.append(wrf_ucan[c_i])
	inmet_smn_i.append(inmet_smn[c_i])
	era5_i.append(era5[c_i])

for c_ii in count_ii:
	reg_usp_ii.append(reg_usp[c_ii])
	reg_ictp_i_ii.append(reg_ictp_i[c_ii])
	reg_ictp_ii_ii.append(reg_ictp_ii[c_ii])
	wrf_ncar_ii.append(wrf_ncar[c_ii])
	wrf_ucan_ii.append(wrf_ucan[c_ii])
	inmet_smn_ii.append(inmet_smn[c_ii])
	era5_ii.append(era5[c_ii])
	
for c_iii in count_iii:
	reg_usp_iii.append(reg_usp[c_iii])
	reg_ictp_i_iii.append(reg_ictp_i[c_iii])
	reg_ictp_ii_iii.append(reg_ictp_ii[c_iii])
	wrf_ncar_iii.append(wrf_ncar[c_iii])
	wrf_ucan_iii.append(wrf_ucan[c_iii])
	inmet_smn_iii.append(inmet_smn[c_iii])
	era5_iii.append(era5[c_iii])
	
for c_iv in count_iv:
	reg_usp_iv.append(reg_usp[c_iv])
	reg_ictp_i_iv.append(reg_ictp_i[c_iv])
	reg_ictp_ii_iv.append(reg_ictp_ii[c_iv])
	wrf_ncar_iv.append(wrf_ncar[c_iv])
	wrf_ucan_iv.append(wrf_ucan[c_iv])
	inmet_smn_iv.append(inmet_smn[c_iv])
	era5_iv.append(era5[c_iv])
	
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
r_reg_ictp_i_era5_ci  = np.corrcoef(reg_ictp_i_c_i, era5_c_i)[0, 1]
r_reg_ictp_ii_era5_ci = np.corrcoef(reg_ictp_ii_c_i, era5_c_i)[0, 1]
r_wrf_ncar_era5_ci    = np.corrcoef(wrf_ncar_c_i, era5_c_i)[0, 1]
r_wrf_ucan_era5_ci    = np.corrcoef(wrf_ucan_c_i, era5_c_i)[0, 1]

# Group II
# Average
reg_usp_c_ii     = np.nanmean(reg_usp_ii, axis=0)
reg_ictp_i_c_ii  = np.nanmean(reg_ictp_i_ii, axis=0)
reg_ictp_ii_c_ii = np.nanmean(reg_ictp_ii_ii, axis=0)
wrf_ncar_c_ii    = np.nanmean(wrf_ncar_ii, axis=0)
wrf_ucan_c_ii    = np.nanmean(wrf_ucan_ii, axis=0)
inmet_smn_c_ii   = np.nanmean(inmet_smn_ii, axis=0)
era5_c_ii        = np.nanmean(era5_ii, axis=0)

# Bias
b_reg_usp_inmet_cii     = np.nanmean(reg_usp_c_ii - inmet_smn_c_ii)
b_reg_ictp_i_inmet_cii  = np.nanmean(reg_ictp_i_c_ii - inmet_smn_c_ii)
b_reg_ictp_ii_inmet_cii = np.nanmean(reg_ictp_ii_c_ii - inmet_smn_c_ii)
b_wrf_ncar_inmet_cii    = np.nanmean(wrf_ncar_c_ii - inmet_smn_c_ii)
b_wrf_ucan_inmet_cii    = np.nanmean(wrf_ucan_c_ii - inmet_smn_c_ii)

b_reg_usp_era5_cii     = np.nanmean(reg_usp_c_ii - era5_c_ii)
b_reg_ictp_i_era5_cii  = np.nanmean(reg_ictp_i_c_ii - era5_c_ii)
b_reg_ictp_ii_era5_cii = np.nanmean(reg_ictp_ii_c_ii - era5_c_ii)
b_wrf_ncar_era5_cii    = np.nanmean(wrf_ncar_c_ii - era5_c_ii)
b_wrf_ucan_era5_cii    = np.nanmean(wrf_ucan_c_ii - era5_c_ii)

# Correlation
r_reg_usp_inmet_cii     = np.corrcoef(reg_usp_c_ii, inmet_smn_c_ii)[0, 1]
r_reg_ictp_i_inmet_cii  = np.corrcoef(reg_ictp_i_c_ii, inmet_smn_c_ii)[0, 1]
r_reg_ictp_ii_inmet_cii = np.corrcoef(reg_ictp_ii_c_ii, inmet_smn_c_ii)[0, 1]
r_wrf_ncar_inmet_cii    = np.corrcoef(wrf_ncar_c_ii, inmet_smn_c_ii)[0, 1]
r_wrf_ucan_inmet_cii    = np.corrcoef(wrf_ucan_c_ii, inmet_smn_c_ii)[0, 1]

r_reg_usp_era5_cii     = np.corrcoef(reg_usp_c_ii, era5_c_ii)[0, 1]
r_reg_ictp_i_era5_cii   = np.corrcoef(reg_ictp_i_c_ii, era5_c_ii)[0, 1]
r_reg_ictp_ii_era5_cii = np.corrcoef(reg_ictp_ii_c_ii, era5_c_ii)[0, 1]
r_wrf_ncar_era5_cii    = np.corrcoef(wrf_ncar_c_ii, era5_c_ii)[0, 1]
r_wrf_ucan_era5_cii    = np.corrcoef(wrf_ucan_c_ii, era5_c_ii)[0, 1]

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
r_wrf_ncar_inmet_ciii    = np.corrcoef(wrf_ncar_c_iii, inmet_smn_c_iii)[0, 1]
r_wrf_ucan_inmet_ciii    = np.corrcoef(wrf_ucan_c_iii, inmet_smn_c_iii)[0, 1]

r_reg_usp_era5_ciii     = np.corrcoef(reg_usp_c_iii, era5_c_iii)[0, 1]
r_reg_ictp_i_era5_ciii  = np.corrcoef(reg_ictp_i_c_iii, era5_c_iii)[0, 1]
r_reg_ictp_ii_era5_ciii = np.corrcoef(reg_ictp_ii_c_iii, era5_c_iii)[0, 1]
r_wrf_ncar_era5_ciii    = np.corrcoef(wrf_ncar_c_iii, era5_c_iii)[0, 1]
r_wrf_ucan_era5_ciii    = np.corrcoef(wrf_ucan_c_iii, era5_c_iii)[0, 1]

# Group IV
# Average
reg_usp_c_iv     = np.nanmean(reg_usp_iv, axis=0)
reg_ictp_i_c_iv  = np.nanmean(reg_ictp_i_iv, axis=0)
reg_ictp_ii_c_iv = np.nanmean(reg_ictp_ii_iv, axis=0)
wrf_ncar_c_iv    = np.nanmean(wrf_ncar_iv, axis=0)
wrf_ucan_c_iv    = np.nanmean(wrf_ucan_iv, axis=0)
inmet_smn_c_iv   = np.nanmean(inmet_smn_iv, axis=0)
era5_c_iv        = np.nanmean(era5_iv, axis=0)

# Bias
b_reg_usp_inmet_civ     = np.nanmean(reg_usp_c_iv - inmet_smn_c_iv)
b_reg_ictp_i_inmet_civ  = np.nanmean(reg_ictp_i_c_iv - inmet_smn_c_iv)
b_reg_ictp_ii_inmet_civ = np.nanmean(reg_ictp_ii_c_iv - inmet_smn_c_iv)
b_wrf_ncar_inmet_civ    = np.nanmean(wrf_ncar_c_iv - inmet_smn_c_iv)
b_wrf_ucan_inmet_civ    = np.nanmean(wrf_ucan_c_iv - inmet_smn_c_iv)

b_reg_usp_era5_civ     = np.nanmean(reg_usp_c_iv - era5_c_iv)
b_reg_ictp_i_era5_civ  = np.nanmean(reg_ictp_i_c_iv - era5_c_iv)
b_reg_ictp_ii_era5_civ = np.nanmean(reg_ictp_ii_c_iv - era5_c_iv)
b_wrf_ncar_era5_civ    = np.nanmean(wrf_ncar_c_iv - era5_c_iv)
b_wrf_ucan_era5_civ    = np.nanmean(wrf_ucan_c_iv - era5_c_iv)

# Correlation
r_reg_usp_inmet_civ     = np.corrcoef(reg_usp_c_iv, inmet_smn_c_iv)[0, 1]
r_reg_ictp_i_inmet_civ  = np.corrcoef(reg_ictp_i_c_iv, inmet_smn_c_iv)[0, 1]
r_reg_ictp_ii_inmet_civ = np.corrcoef(reg_ictp_ii_c_iv, inmet_smn_c_iv)[0, 1]
r_wrf_ncar_inmet_civ    = np.corrcoef(wrf_ncar_c_iv, inmet_smn_c_iv)[0, 1]
r_wrf_ucan_inmet_civ    = np.corrcoef(wrf_ucan_c_iv, inmet_smn_c_iv)[0, 1]

r_reg_usp_era5_civ     = np.corrcoef(reg_usp_c_iv, era5_c_iv)[0, 1]
r_reg_ictp_i_era5_civ  = np.corrcoef(reg_ictp_i_c_iv, era5_c_iv)[0, 1]
r_reg_ictp_ii_era5_civ = np.corrcoef(reg_ictp_ii_c_iv, era5_c_iv)[0, 1]
r_wrf_ncar_era5_civ    = np.corrcoef(wrf_ncar_c_iv, era5_c_iv)[0, 1]
r_wrf_ucan_era5_civ    = np.corrcoef(wrf_ucan_c_iv, era5_c_iv)[0, 1]

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
r_reg_ictp_i_era5_cv  = np.corrcoef(reg_ictp_i_c_v, era5_c_v)[0, 1]
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
fig = plt.figure(figsize=(12, 8))
font_size = 8

dt = pd.date_range(start="20180601", end="20210531", freq="ME")

ax = fig.add_subplot(3, 2, 1)
reg_usp_c_i_dt     = pd.Series(data=reg_usp_c_i, index=dt)
reg_ictp_i_c_i_dt  = pd.Series(data=reg_ictp_i_c_i, index=dt)
reg_ictp_ii_c_i_dt = pd.Series(data=reg_ictp_i_c_i, index=dt)
wrf_ncar_c_i_dt    = pd.Series(data=wrf_ncar_c_i, index=dt)
wrf_ucan_c_i_dt    = pd.Series(data=wrf_ucan_c_i, index=dt)
inmet_smn_c_i_dt   = pd.Series(data=inmet_smn_c_i, index=dt)
era5_c_i_dt        = pd.Series(data=era5_c_i, index=dt)
plt.plot(inmet_smn_c_i_dt,   linewidth=1, color='black', label = 'INMET+SMN')
plt.plot(era5_c_i_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_i_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_i_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_i_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_i_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_i_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(0.1, 0.95, f"Reg4 = {b_reg_usp_inmet_ci:.2f}({r_reg_usp_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.85, f"Reg5-Holt = {b_reg_ictp_i_inmet_ci:.2f}({r_reg_ictp_i_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.75, f"Reg5-UW = {b_reg_ictp_ii_inmet_ci:.2f}({r_reg_ictp_ii_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.65, f"WRF-NCAR = {b_wrf_ncar_inmet_ci:.2f}({r_wrf_ncar_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.55, f"WRF-UCAN = {b_wrf_ucan_inmet_ci:.2f}({r_wrf_ucan_inmet_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.5, 0.95, f"Reg4 = {b_reg_usp_era5_ci:.2f}({r_reg_usp_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.85, f"Reg5-Holt = {b_reg_ictp_i_era5_ci:.2f}({r_reg_ictp_i_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.75, f"Reg5-UW = {b_reg_ictp_ii_era5_ci:.2f}({r_reg_ictp_ii_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.65, f"WRF-NCAR = {b_wrf_ncar_era5_ci:.2f}({r_wrf_ncar_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.55, f"WRF-UCAN = {b_wrf_ucan_era5_ci:.2f}({r_wrf_ucan_era5_ci:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=font_size, fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1), fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(3, 2, 2)
reg_usp_c_ii_dt     = pd.Series(data=reg_usp_c_ii, index=dt)
reg_ictp_i_c_ii_dt  = pd.Series(data=reg_ictp_i_c_ii, index=dt)
reg_ictp_ii_c_ii_dt = pd.Series(data=reg_ictp_ii_c_ii, index=dt)
wrf_ncar_c_ii_dt    = pd.Series(data=wrf_ncar_c_ii, index=dt)
wrf_ucan_c_ii_dt    = pd.Series(data=wrf_ucan_c_ii, index=dt)
inmet_smn_c_ii_dt   = pd.Series(data=inmet_smn_c_ii, index=dt)
era5_c_ii_dt        = pd.Series(data=era5_c_ii, index=dt)
plt.plot(inmet_smn_c_ii_dt,   linewidth=1, color='black', label = 'INMET+SMN')
plt.plot(era5_c_ii_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_ii_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_ii_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_ii_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_ii_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_ii_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(0.1, 0.95, f"Reg4 = {b_reg_usp_inmet_cii:.2f}({r_reg_usp_inmet_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.85, f"Reg5-Holt = {b_reg_ictp_i_inmet_cii:.2f}({r_reg_ictp_i_inmet_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.75, f"Reg5-UW = {b_reg_ictp_ii_inmet_cii:.2f}({r_reg_ictp_ii_inmet_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.65, f"WRF-NCAR = {b_wrf_ncar_inmet_cii:.2f}({r_wrf_ncar_inmet_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.55, f"WRF-UCAN = {b_wrf_ucan_inmet_cii:.2f}({r_wrf_ucan_inmet_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.5, 0.95, f"Reg4 = {b_reg_usp_era5_cii:.2f}({r_reg_usp_era5_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.85, f"Reg5-Holt = {b_reg_ictp_i_era5_cii:.2f}({r_reg_ictp_i_era5_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.75, f"Reg5-UW = {b_reg_ictp_ii_era5_cii:.2f}({r_reg_ictp_ii_era5_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.65, f"WRF-NCAR = {b_wrf_ncar_era5_cii:.2f}({r_wrf_ncar_era5_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.55, f"WRF-UCAN = {b_wrf_ucan_era5_cii:.2f}({r_wrf_ucan_era5_cii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(b) Cluster II', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1), fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(3, 2, 3)
reg_usp_c_iii_dt     = pd.Series(data=reg_usp_c_iii, index=dt)
reg_ictp_i_c_iii_dt  = pd.Series(data=reg_ictp_i_c_iii, index=dt)
reg_ictp_ii_c_iii_dt = pd.Series(data=reg_ictp_ii_c_iii, index=dt)
wrf_ncar_c_iii_dt    = pd.Series(data=wrf_ncar_c_iii, index=dt)
wrf_ucan_c_iii_dt    = pd.Series(data=wrf_ucan_c_iii, index=dt)
inmet_smn_c_iii_dt   = pd.Series(data=inmet_smn_c_iii, index=dt)
era5_c_iii_dt        = pd.Series(data=era5_c_iii, index=dt)
plt.plot(inmet_smn_c_iii_dt,   linewidth=1, color='black', label = 'INMET+SMN')
plt.plot(era5_c_iii_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_iii_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_iii_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_iii_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_iii_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_iii_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(0.1, 0.95, f"Reg4 = {b_reg_usp_inmet_ciii:.2f}({r_reg_usp_inmet_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.85, f"Reg5-Holt = {b_reg_ictp_i_inmet_ciii:.2f}({r_reg_ictp_i_inmet_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.75, f"Reg5-UW = {b_reg_ictp_ii_inmet_ciii:.2f}({r_reg_ictp_ii_inmet_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.65, f"WRF-NCAR = {b_wrf_ncar_inmet_ciii:.2f}({r_wrf_ncar_inmet_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.55, f"WRF-UCAN = {b_wrf_ucan_inmet_ciii:.2f}({r_wrf_ucan_inmet_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.5, 0.95, f"Reg4 = {b_reg_usp_era5_ciii:.2f}({r_reg_usp_era5_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.85, f"Reg5-Holt = {b_reg_ictp_i_era5_ciii:.2f}({r_reg_ictp_i_era5_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.75, f"Reg5-UW = {b_reg_ictp_ii_era5_ciii:.2f}({r_reg_ictp_ii_era5_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.65, f"WRF-NCAR = {b_wrf_ncar_era5_ciii:.2f}({r_wrf_ncar_era5_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.55, f"WRF-UCAN = {b_wrf_ucan_era5_ciii:.2f}({r_wrf_ucan_era5_ciii:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(c) Cluster III', loc='left', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=font_size, fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1), fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(3, 2, 4)
reg_usp_c_iv_dt     = pd.Series(data=reg_usp_c_iv, index=dt)
reg_ictp_i_c_iv_dt  = pd.Series(data=reg_ictp_i_c_iv, index=dt)
reg_ictp_ii_c_iv_dt = pd.Series(data=reg_ictp_ii_c_iv, index=dt)
wrf_ncar_c_iv_dt    = pd.Series(data=wrf_ncar_c_iv, index=dt)
wrf_ucan_c_iv_dt    = pd.Series(data=wrf_ucan_c_iv, index=dt)
inmet_smn_c_iv_dt   = pd.Series(data=inmet_smn_c_iv, index=dt)
era5_c_iv_dt        = pd.Series(data=era5_c_iv, index=dt)
plt.plot(inmet_smn_c_iv_dt,   linewidth=1, color='black', label = 'INMET+SMN')
plt.plot(era5_c_iv_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_iv_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_iv_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_iv_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_iv_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_iv_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(0.1, 0.95, f"Reg4 = {b_reg_usp_inmet_civ:.2f}({r_reg_usp_inmet_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.85, f"Reg5-Holt = {b_reg_ictp_i_inmet_civ:.2f}({r_reg_ictp_i_inmet_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.75, f"Reg5-UW = {b_reg_ictp_ii_inmet_civ:.2f}({r_reg_ictp_ii_inmet_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.65, f"WRF-NCAR = {b_wrf_ncar_inmet_civ:.2f}({r_wrf_ncar_inmet_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.55, f"WRF-UCAN = {b_wrf_ucan_inmet_civ:.2f}({r_wrf_ucan_inmet_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.5, 0.95, f"Reg4 = {b_reg_usp_era5_civ:.2f}({r_reg_usp_era5_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.85, f"Reg5-Holt = {b_reg_ictp_i_era5_civ:.2f}({r_reg_ictp_i_era5_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.75, f"Reg5-UW = {b_reg_ictp_ii_era5_civ:.2f}({r_reg_ictp_ii_era5_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.65, f"WRF-NCAR = {b_wrf_ncar_era5_civ:.2f}({r_wrf_ncar_era5_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.55, f"WRF-UCAN = {b_wrf_ucan_era5_civ:.2f}({r_wrf_ucan_era5_civ:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(d) Cluster IV', loc='left', fontsize=font_size, fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1), fontsize=font_size)
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(3, 2, 5)
reg_usp_c_v_dt     = pd.Series(data=reg_usp_c_v, index=dt)
reg_ictp_i_c_v_dt  = pd.Series(data=reg_ictp_i_c_v, index=dt)
reg_ictp_ii_c_v_dt = pd.Series(data=reg_ictp_ii_c_v, index=dt)
wrf_ncar_c_v_dt    = pd.Series(data=wrf_ncar_c_v, index=dt)
wrf_ucan_c_v_dt    = pd.Series(data=wrf_ucan_c_v, index=dt)
inmet_smn_c_v_dt   = pd.Series(data=inmet_smn_c_v, index=dt)
era5_c_v_dt        = pd.Series(data=era5_c_v, index=dt)
plt.plot(inmet_smn_c_v_dt,   linewidth=1, color='black', label = 'INMET+SMN')
plt.plot(era5_c_v_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_v_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_v_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_v_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_v_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_v_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(0.1, 0.95, f"Reg4 = {b_reg_usp_inmet_cv:.2f}({r_reg_usp_inmet_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.85, f"Reg5-Holt = {b_reg_ictp_i_inmet_cv:.2f}({r_reg_ictp_i_inmet_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.75, f"Reg5-UW = {b_reg_ictp_ii_inmet_cv:.2f}({r_reg_ictp_ii_inmet_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.65, f"WRF-NCAR = {b_wrf_ncar_inmet_cv:.2f}({r_wrf_ncar_inmet_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.55, f"WRF-UCAN = {b_wrf_ucan_inmet_cv:.2f}({r_wrf_ucan_inmet_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.5, 0.95, f"Reg4 = {b_reg_usp_era5_cv:.2f}({r_reg_usp_era5_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.85, f"Reg5-Holt = {b_reg_ictp_i_era5_cv:.2f}({r_reg_ictp_i_era5_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.75, f"Reg5-UW = {b_reg_ictp_ii_era5_cv:.2f}({r_reg_ictp_ii_era5_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.65, f"WRF-NCAR = {b_wrf_ncar_era5_cv:.2f}({r_wrf_ncar_era5_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.55, f"WRF-UCAN = {b_wrf_ucan_era5_cv:.2f}({r_wrf_ucan_era5_cv:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(e) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Period 2018/06 - 2021/05', fontsize=font_size, fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontsize=font_size, fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1), fontsize=font_size)
plt.xticks(fontsize=7)

ax = fig.add_subplot(3, 2, 6)
reg_usp_c_dt     = pd.Series(data=reg_usp_c, index=dt)
reg_ictp_i_c_dt  = pd.Series(data=reg_ictp_i_c, index=dt)
reg_ictp_ii_c_dt = pd.Series(data=reg_ictp_ii_c, index=dt)
wrf_ncar_c_dt    = pd.Series(data=wrf_ncar_c, index=dt)
wrf_ucan_c_dt    = pd.Series(data=wrf_ucan_c, index=dt)
inmet_smn_c_dt   = pd.Series(data=inmet_smn_c, index=dt)
era5_c_dt        = pd.Series(data=era5_c, index=dt)
plt.plot(inmet_smn_c_dt,   linewidth=1, color='black', label = 'INMET+SMN')
plt.plot(era5_c_dt,        linewidth=1, color='red', label = 'ERA5')
plt.plot(reg_usp_c_dt,     linewidth=1, linestyle='--', color='blue', label = 'Reg4')
plt.plot(reg_ictp_i_c_dt,  linewidth=1, linestyle='--', color='gray', label = 'Reg5-Holt')
plt.plot(reg_ictp_ii_c_dt, linewidth=1, linestyle='--', color='brown', label = 'Reg5-UW')
plt.plot(wrf_ncar_c_dt,    linewidth=1, linestyle='--', color='green', label = 'WRF-NCAR')
plt.plot(wrf_ucan_c_dt,    linewidth=1, linestyle='--', color='orange', label = 'WRF-UCAN')
ax.text(0.1, 0.95, f"Reg4 = {b_reg_usp_inmet_c:.2f}({r_reg_usp_inmet_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.85, f"Reg5-Holt = {b_reg_ictp_i_inmet_c:.2f}({r_reg_ictp_i_inmet_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.75, f"Reg5-UW = {b_reg_ictp_ii_inmet_c:.2f}({r_reg_ictp_ii_inmet_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.65, f"WRF-NCAR = {b_wrf_ncar_inmet_c:.2f}({r_wrf_ncar_inmet_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.1, 0.55, f"WRF-UCAN = {b_wrf_ucan_inmet_c:.2f}({r_wrf_ucan_inmet_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='black')
ax.text(0.5, 0.95, f"Reg4 = {b_reg_usp_era5_c:.2f}({r_reg_usp_era5_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.85, f"Reg5-Holt = {b_reg_ictp_i_era5_c:.2f}({r_reg_ictp_i_era5_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.75, f"Reg5-UW = {b_reg_ictp_ii_era5_c:.2f}({r_reg_ictp_ii_era5_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.65, f"WRF-NCAR = {b_wrf_ncar_era5_c:.2f}({r_wrf_ncar_era5_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
ax.text(0.5, 0.55, f"WRF-UCAN = {b_wrf_ucan_era5_c:.2f}({r_wrf_ucan_era5_c:.2f})", transform=ax.transAxes, ha='left', va='top', fontsize=font_size, color='red')
plt.title('(f) All clusters', loc='left', fontsize=font_size, fontweight='bold')
plt.xlabel('Period 2018/06 - 2021/05', fontsize=font_size, fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1), fontsize=font_size)
plt.xticks(fontsize=7)

plt.legend(ncol=7, fontsize=font_size, loc=(-1., -0.35))

# Path out to save figure
path_out = '{0}/figs/paper_cp'.format(path)
name_out = 'pyplt_annual_variability_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()


