# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot cdf function"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii


path = '/afs/ictp.it/home/m/mda_silv/Documents'


def import_inmet():
	
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []
	mean_v = []
	mean_vi = []
	mean_vii = []

	# Select lat and lon 
	for i in range(1, 100):
		yy=inmet[i][2]
		xx=inmet[i][3]
		
		print('Reading weather station:', i, inmet[i][0])		
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.values
		mean_i.append(d_i*86400)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_day_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_day_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.values
		mean_v.append(d_v*86400)
		
		# Reading inmet 
		d_vi = xr.open_dataset('{0}/FPS_SESA/database/obs/inmet/inmet_nc_sesa/pre/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.resample(time='1D').sum()
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/FPS_SESA/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_day_20180101-20211231.nc')
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
	for i in range(1, 73):
		yy=smn_i[i][1]
		xx=smn_i[i][2]
		
		print('Reading weather station:', i, smn_i[i][0])	
		# reading regcm usp 
		d_i = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.values
		mean_i.append(d_i*86400)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_day_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_day_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.values
		mean_v.append(d_v*86400)
						
		# Reading smn 
		d_vi = xr.open_dataset('{0}/FPS_SESA/database/obs/smn_i/smn_nc/'.format(path) + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.resample(time='1D').sum()
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/FPS_SESA/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_day_20180101-20211231.nc')
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
		d_i = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_usp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_i = d_i.values
		mean_i.append(d_i*86400)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_day_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('{0}/FPS_SESA/database/rcm/reg_ictp/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_day_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ncar/'.format(path) + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('{0}/FPS_SESA/database/rcm/wrf_ucan/'.format(path) + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_v = d_v.values
		mean_v.append(d_v*86400)
						
		# Reading smn 
		d_vi = xr.open_dataset('{0}/FPS_SESA/database/obs/smn_ii/smn_nc/pre/'.format(path) + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('{0}/FPS_SESA/database/obs/era5/'.format(path) + 'tp_era5_csam_4km_day_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=slice(yy-0.04,yy+0.04),lon=slice(xx-0.04,xx+0.04)).mean(('lat','lon'))
		d_vii = d_vii.values
		mean_vii.append(d_vii)
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii
	
	
var = 'pr'

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

list_hc = [4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 0, 4, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0, 
4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0,
0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 3, 0, 1, 1, 0, 1, 0, 1, 3, 3, 3, 3, 1, 1, 1, 3,
1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 
3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 
3, 1, 2, 4, 0, 4, 4, 4, 3, 4, 4, 4, 1, 4, 4, 2, 4, 4, 3, 1, 3, 3, 3, 3, 0, 3,
3, 3, 3, 1, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 0, 2, 2,
1, 3, 3, 1, 3, 0, 2, 1, 3, 3, 2, 2, 2, 2, 2, 2, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2,
2, 1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 3, 0, 1, 1, 2, 1, 2, 0, 2, 2, 2, 1, 2,
2, 2, 0, 4, 3, 4, 4, 3, 4, 4, 0, 4, 0, 4, 4, 4, 3, 4, 4, 1, 2, 2, 0, 1, 1, 0]

print(len(reg_usp))
print(len(reg_ictp_i))
print(len(reg_ictp_ii))
print(len(wrf_ncar))
print(len(wrf_ucan))
print(len(inmet_smn))
print(len(era5))
print(len(list_hc))
print()

count_i = []
count_ii = []
count_iii = []
count_iv = []
count_v = []

# Select cluster
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

reg_usp_c_i     = np.nanmean(reg_usp_i, axis=0)
reg_ictp_i_c_i  = np.nanmean(reg_ictp_i_i, axis=0)
reg_ictp_ii_c_i = np.nanmean(reg_ictp_ii_i, axis=0)
wrf_ncar_c_i    = np.nanmean(wrf_ncar_i, axis=0)
wrf_ucan_c_i    = np.nanmean(wrf_ucan_i, axis=0)
inmet_smn_c_i   = np.nanmean(inmet_smn_i, axis=0)
era5_c_i        = np.nanmean(era5_i, axis=0)

reg_usp_c_ii     = np.nanmean(reg_usp_ii, axis=0)
reg_ictp_i_c_ii  = np.nanmean(reg_ictp_i_ii, axis=0)
reg_ictp_ii_c_ii = np.nanmean(reg_ictp_ii_ii, axis=0)
wrf_ncar_c_ii    = np.nanmean(wrf_ncar_ii, axis=0)
wrf_ucan_c_ii    = np.nanmean(wrf_ucan_ii, axis=0)
inmet_smn_c_ii   = np.nanmean(inmet_smn_ii, axis=0)
era5_c_ii        = np.nanmean(era5_ii, axis=0)

reg_usp_c_iii     = np.nanmean(reg_usp_iii, axis=0)
reg_ictp_i_c_iii  = np.nanmean(reg_ictp_i_iii, axis=0)
reg_ictp_ii_c_iii = np.nanmean(reg_ictp_ii_iii, axis=0)
wrf_ncar_c_iii    = np.nanmean(wrf_ncar_iii, axis=0)
wrf_ucan_c_iii    = np.nanmean(wrf_ucan_iii, axis=0)
inmet_smn_c_iii   = np.nanmean(inmet_smn_iii, axis=0)
era5_c_iii        = np.nanmean(era5_iii, axis=0)

reg_usp_c_iv     = np.nanmean(reg_usp_iv, axis=0)
reg_ictp_i_c_iv  = np.nanmean(reg_ictp_i_iv, axis=0)
reg_ictp_ii_c_iv = np.nanmean(reg_ictp_ii_iv, axis=0)
wrf_ncar_c_iv    = np.nanmean(wrf_ncar_iv, axis=0)
wrf_ucan_c_iv    = np.nanmean(wrf_ucan_iv, axis=0)
inmet_smn_c_iv   = np.nanmean(inmet_smn_iv, axis=0)
era5_c_iv        = np.nanmean(era5_iv, axis=0)

reg_usp_c_v     = np.nanmean(reg_usp_v, axis=0)
reg_ictp_i_c_v  = np.nanmean(reg_ictp_i_v, axis=0)
reg_ictp_ii_c_v = np.nanmean(reg_ictp_ii_v, axis=0)
wrf_ncar_c_v    = np.nanmean(wrf_ncar_v, axis=0)
wrf_ucan_c_v    = np.nanmean(wrf_ucan_v, axis=0)
inmet_smn_c_v   = np.nanmean(inmet_smn_v, axis=0)
era5_c_v        = np.nanmean(era5_v, axis=0)

reg_usp_c_vi     = np.nanmean(reg_usp, axis=0)
reg_ictp_i_c_vi  = np.nanmean(reg_ictp_i, axis=0)
reg_ictp_ii_c_vi = np.nanmean(reg_ictp_ii, axis=0)
wrf_ncar_c_vi    = np.nanmean(wrf_ncar, axis=0)
wrf_ucan_c_vi    = np.nanmean(wrf_ucan, axis=0)
inmet_smn_c_vi   = np.nanmean(inmet_smn, axis=0)
era5_c_vi        = np.nanmean(era5, axis=0)

# Round values to each cluster
round_reg_usp_c_i     = np.round(reg_usp_c_i,0)
round_reg_ictp_i_c_i  = np.round(reg_ictp_i_c_i,0)
round_reg_ictp_ii_c_i = np.round(reg_ictp_ii_c_i,0)
round_wrf_ncar_c_i    = np.round(wrf_ncar_c_i,0)
round_wrf_ucan_c_i    = np.round(wrf_ucan_c_i,0)
round_inmet_smn_c_i   = np.round(inmet_smn_c_i,0)
round_era5_c_i        = np.round(era5_c_i,0)

round_reg_usp_c_ii     = np.round(reg_usp_c_ii,0)
round_reg_ictp_i_c_ii  = np.round(reg_ictp_i_c_ii,0)
round_reg_ictp_ii_c_ii = np.round(reg_ictp_ii_c_ii,0)
round_wrf_ncar_c_ii    = np.round(wrf_ncar_c_ii,0)
round_wrf_ucan_c_ii    = np.round(wrf_ucan_c_ii,0)
round_inmet_smn_c_ii   = np.round(inmet_smn_c_ii,0)
round_era5_c_ii        = np.round(era5_c_ii,0)

round_reg_usp_c_iii     = np.round(reg_usp_c_iii,0)
round_reg_ictp_i_c_iii  = np.round(reg_ictp_i_c_iii,0)
round_reg_ictp_ii_c_iii = np.round(reg_ictp_ii_c_iii,0)
round_wrf_ncar_c_iii    = np.round(wrf_ncar_c_iii,0)
round_wrf_ucan_c_iii    = np.round(wrf_ucan_c_iii,0)
round_inmet_smn_c_iii   = np.round(inmet_smn_c_iii,0)
round_era5_c_iii        = np.round(era5_c_iii,0)

round_reg_usp_c_iv     = np.round(reg_usp_c_iv,0)
round_reg_ictp_i_c_iv  = np.round(reg_ictp_i_c_iv,0)
round_reg_ictp_ii_c_iv = np.round(reg_ictp_ii_c_iv,0)
round_wrf_ncar_c_iv    = np.round(wrf_ncar_c_iv,0)
round_wrf_ucan_c_iv    = np.round(wrf_ucan_c_iv,0)
round_inmet_smn_c_iv   = np.round(inmet_smn_c_iv,0)
round_era5_c_iv        = np.round(era5_c_iv,0)

round_reg_usp_c_v     = np.round(reg_usp_c_v,0)
round_reg_ictp_i_c_v  = np.round(reg_ictp_i_c_v,0)
round_reg_ictp_ii_c_v = np.round(reg_ictp_ii_c_v,0)
round_wrf_ncar_c_v    = np.round(wrf_ncar_c_v,0)
round_wrf_ucan_c_v    = np.round(wrf_ucan_c_v,0)
round_inmet_smn_c_v   = np.round(inmet_smn_c_v,0)
round_era5_c_v        = np.round(era5_c_v,0)

round_reg_usp_c_vi     = np.round(reg_usp_c_vi,0)
round_reg_ictp_i_c_vi  = np.round(reg_ictp_i_c_vi,0)
round_reg_ictp_ii_c_vi = np.round(reg_ictp_ii_c_vi,0)
round_wrf_ncar_c_vi    = np.round(wrf_ncar_c_vi,0)
round_wrf_ucan_c_vi    = np.round(wrf_ucan_c_vi,0)
round_inmet_smn_c_vi   = np.round(inmet_smn_c_vi,0)
round_era5_c_vi        = np.round(era5_c_vi,0)

# Filter 0 mm/day
filter_reg_usp_c_i     = round_reg_usp_c_i[round_reg_usp_c_i > 0.]
filter_reg_ictp_i_c_i  = round_reg_ictp_i_c_i[round_reg_ictp_i_c_i > 0.]
filter_reg_ictp_ii_c_i = round_reg_ictp_ii_c_i[round_reg_ictp_ii_c_i > 0.]
filter_wrf_ncar_c_i    = round_wrf_ncar_c_i[round_wrf_ncar_c_i > 0.]
filter_wrf_ucan_c_i    = round_wrf_ucan_c_i[round_wrf_ucan_c_i > 0.]
filter_inmet_smn_c_i   = round_inmet_smn_c_i[round_inmet_smn_c_i > 0.]
filter_era5_c_i        = round_era5_c_i[round_era5_c_i > 0.]

filter_reg_usp_c_ii     = round_reg_usp_c_ii[round_reg_usp_c_ii > 0.]
filter_reg_ictp_i_c_ii  = round_reg_ictp_i_c_ii[round_reg_ictp_i_c_ii > 0.]
filter_reg_ictp_ii_c_ii = round_reg_ictp_ii_c_ii[round_reg_ictp_ii_c_ii > 0.]
filter_wrf_ncar_c_ii    = round_wrf_ncar_c_ii[round_wrf_ncar_c_ii > 0.]
filter_wrf_ucan_c_ii    = round_wrf_ucan_c_ii[round_wrf_ucan_c_ii > 0.]
filter_inmet_smn_c_ii   = round_inmet_smn_c_ii[round_inmet_smn_c_ii > 0.]
filter_era5_c_ii        = round_era5_c_ii[round_era5_c_ii > 0.]

filter_reg_usp_c_iii     = round_reg_usp_c_iii[round_reg_usp_c_iii > 0.]
filter_reg_ictp_i_c_iii  = round_reg_ictp_i_c_iii[round_reg_ictp_i_c_iii > 0.]
filter_reg_ictp_ii_c_iii = round_reg_ictp_ii_c_iii[round_reg_ictp_ii_c_iii > 0.]
filter_wrf_ncar_c_iii    = round_wrf_ncar_c_iii[round_wrf_ncar_c_iii > 0.]
filter_wrf_ucan_c_iii    = round_wrf_ucan_c_iii[round_wrf_ucan_c_iii > 0.]
filter_inmet_smn_c_iii   = round_inmet_smn_c_iii[round_inmet_smn_c_iii > 0.]
filter_era5_c_iii        = round_era5_c_iii[round_era5_c_iii > 0.]

filter_reg_usp_c_iv     = round_reg_usp_c_iv[round_reg_usp_c_iv > 0.]
filter_reg_ictp_i_c_iv  = round_reg_ictp_i_c_iv[round_reg_ictp_i_c_iv > 0.]
filter_reg_ictp_ii_c_iv = round_reg_ictp_ii_c_iv[round_reg_ictp_ii_c_iv > 0.]
filter_wrf_ncar_c_iv    = round_wrf_ncar_c_iv[round_wrf_ncar_c_iv > 0.]
filter_wrf_ucan_c_iv    = round_wrf_ucan_c_iv[round_wrf_ucan_c_iv > 0.]
filter_inmet_smn_c_iv   = round_inmet_smn_c_iv[round_inmet_smn_c_iv > 0.]
filter_era5_c_iv        = round_era5_c_iv[round_era5_c_iv > 0.]

filter_reg_usp_c_v     = round_reg_usp_c_v[round_reg_usp_c_v > 0.]
filter_reg_ictp_i_c_v  = round_reg_ictp_i_c_v[round_reg_ictp_i_c_v > 0.]
filter_reg_ictp_ii_c_v = round_reg_ictp_ii_c_v[round_reg_ictp_ii_c_v > 0.]
filter_wrf_ncar_c_v    = round_wrf_ncar_c_v[round_wrf_ncar_c_v > 0.]
filter_wrf_ucan_c_v    = round_wrf_ucan_c_v[round_wrf_ucan_c_v > 0.]
filter_inmet_smn_c_v   = round_inmet_smn_c_v[round_inmet_smn_c_v > 0.]
filter_era5_c_v        = round_era5_c_v[round_era5_c_v > 0.]

filter_reg_usp_c_vi     = round_reg_usp_c_vi[round_reg_usp_c_vi > 0.]
filter_reg_ictp_i_c_vi  = round_reg_ictp_i_c_vi[round_reg_ictp_i_c_vi > 0.]
filter_reg_ictp_ii_c_vi = round_reg_ictp_ii_c_vi[round_reg_ictp_ii_c_vi > 0.]
filter_wrf_ncar_c_vi    = round_wrf_ncar_c_vi[round_wrf_ncar_c_vi > 0.]
filter_wrf_ucan_c_vi    = round_wrf_ucan_c_vi[round_wrf_ucan_c_vi > 0.]
filter_inmet_smn_c_vi   = round_inmet_smn_c_vi[round_inmet_smn_c_vi > 0.]
filter_era5_c_vi        = round_era5_c_vi[round_era5_c_vi > 0.]

# Compute frequency
x_freq_reg_usp_c_i,     freq_reg_usp_c_i     = np.unique(filter_reg_usp_c_i,     return_counts=True) 
x_freq_reg_ictp_i_c_i,  freq_reg_ictp_i_c_i  = np.unique(filter_reg_ictp_i_c_i,  return_counts=True) 
x_freq_reg_ictp_ii_c_i, freq_reg_ictp_ii_c_i = np.unique(filter_reg_ictp_ii_c_i, return_counts=True) 
x_freq_wrf_ncar_c_i,    freq_wrf_ncar_c_i    = np.unique(filter_wrf_ncar_c_i,    return_counts=True) 
x_freq_wrf_ucan_c_i,    freq_wrf_ucan_c_i    = np.unique(filter_wrf_ucan_c_i,    return_counts=True) 
x_freq_inmet_smn_c_i,   freq_inmet_smn_c_i   = np.unique(filter_inmet_smn_c_i,   return_counts=True) 
x_freq_era5_c_i,        freq_era5_c_i        = np.unique(filter_era5_c_i,        return_counts=True) 

x_freq_reg_usp_c_ii,     freq_reg_usp_c_ii     = np.unique(filter_reg_usp_c_ii,     return_counts=True) 
x_freq_reg_ictp_i_c_ii,  freq_reg_ictp_i_c_ii  = np.unique(filter_reg_ictp_i_c_ii,  return_counts=True) 
x_freq_reg_ictp_ii_c_ii, freq_reg_ictp_ii_c_ii = np.unique(filter_reg_ictp_ii_c_ii, return_counts=True) 
x_freq_wrf_ncar_c_ii,    freq_wrf_ncar_c_ii    = np.unique(filter_wrf_ncar_c_ii,    return_counts=True) 
x_freq_wrf_ucan_c_ii,    freq_wrf_ucan_c_ii    = np.unique(filter_wrf_ucan_c_ii,    return_counts=True) 
x_freq_inmet_smn_c_ii,   freq_inmet_smn_c_ii   = np.unique(filter_inmet_smn_c_ii,   return_counts=True) 
x_freq_era5_c_ii,        freq_era5_c_ii        = np.unique(filter_era5_c_ii,        return_counts=True) 

x_freq_reg_usp_c_iii,     freq_reg_usp_c_iii     = np.unique(filter_reg_usp_c_iii,     return_counts=True) 
x_freq_reg_ictp_i_c_iii,  freq_reg_ictp_i_c_iii  = np.unique(filter_reg_ictp_i_c_iii,  return_counts=True) 
x_freq_reg_ictp_ii_c_iii, freq_reg_ictp_ii_c_iii = np.unique(filter_reg_ictp_ii_c_iii, return_counts=True) 
x_freq_wrf_ncar_c_iii,    freq_wrf_ncar_c_iii    = np.unique(filter_wrf_ncar_c_iii,    return_counts=True) 
x_freq_wrf_ucan_c_iii,    freq_wrf_ucan_c_iii    = np.unique(filter_wrf_ucan_c_iii,    return_counts=True) 
x_freq_inmet_smn_c_iii,   freq_inmet_smn_c_iii   = np.unique(filter_inmet_smn_c_iii,   return_counts=True) 
x_freq_era5_c_iii,        freq_era5_c_iii        = np.unique(filter_era5_c_iii,        return_counts=True) 

x_freq_reg_usp_c_iv,     freq_reg_usp_c_iv     = np.unique(filter_reg_usp_c_iv,     return_counts=True) 
x_freq_reg_ictp_i_c_iv,  freq_reg_ictp_i_c_iv  = np.unique(filter_reg_ictp_i_c_iv,  return_counts=True) 
x_freq_reg_ictp_ii_c_iv, freq_reg_ictp_ii_c_iv = np.unique(filter_reg_ictp_ii_c_iv, return_counts=True) 
x_freq_wrf_ncar_c_iv,    freq_wrf_ncar_c_iv    = np.unique(filter_wrf_ncar_c_iv,    return_counts=True) 
x_freq_wrf_ucan_c_iv,    freq_wrf_ucan_c_iv    = np.unique(filter_wrf_ucan_c_iv,    return_counts=True) 
x_freq_inmet_smn_c_iv,   freq_inmet_smn_c_iv   = np.unique(filter_inmet_smn_c_iv,   return_counts=True) 
x_freq_era5_c_iv,        freq_era5_c_iv        = np.unique(filter_era5_c_iv,        return_counts=True)

x_freq_reg_usp_c_v,     freq_reg_usp_c_v     = np.unique(filter_reg_usp_c_v,     return_counts=True) 
x_freq_reg_ictp_i_c_v,  freq_reg_ictp_i_c_v  = np.unique(filter_reg_ictp_i_c_v,  return_counts=True) 
x_freq_reg_ictp_ii_c_v, freq_reg_ictp_ii_c_v = np.unique(filter_reg_ictp_ii_c_v, return_counts=True) 
x_freq_wrf_ncar_c_v,    freq_wrf_ncar_c_v    = np.unique(filter_wrf_ncar_c_v,    return_counts=True) 
x_freq_wrf_ucan_c_v,    freq_wrf_ucan_c_v    = np.unique(filter_wrf_ucan_c_v,    return_counts=True) 
x_freq_inmet_smn_c_v,   freq_inmet_smn_c_v   = np.unique(filter_inmet_smn_c_v,   return_counts=True) 
x_freq_era5_c_v,        freq_era5_c_v        = np.unique(filter_era5_c_v,        return_counts=True)

x_freq_reg_usp_c_vi,     freq_reg_usp_c_vi     = np.unique(filter_reg_usp_c_vi,     return_counts=True) 
x_freq_reg_ictp_i_c_vi,  freq_reg_ictp_i_c_vi  = np.unique(filter_reg_ictp_i_c_vi,  return_counts=True) 
x_freq_reg_ictp_ii_c_vi, freq_reg_ictp_ii_c_vi = np.unique(filter_reg_ictp_ii_c_vi, return_counts=True) 
x_freq_wrf_ncar_c_vi,    freq_wrf_ncar_c_vi    = np.unique(filter_wrf_ncar_c_vi,    return_counts=True) 
x_freq_wrf_ucan_c_vi,    freq_wrf_ucan_c_vi    = np.unique(filter_wrf_ucan_c_vi,    return_counts=True) 
x_freq_inmet_smn_c_vi,   freq_inmet_smn_c_vi   = np.unique(filter_inmet_smn_c_vi,   return_counts=True) 
x_freq_era5_c_vi,        freq_era5_c_vi        = np.unique(filter_era5_c_vi,        return_counts=True)

# Plot figure
fig = plt.figure(figsize=(10, 10))
font_size = 8

ax = fig.add_subplot(3, 2, 1)  
plt.plot(x_freq_inmet_smn_c_i,   freq_inmet_smn_c_i,   marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='INMET+SMN')
plt.plot(x_freq_era5_c_i,        freq_era5_c_i,        marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='ERA5')
plt.plot(x_freq_reg_usp_c_i,     freq_reg_usp_c_i,     marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='Reg4')
plt.plot(x_freq_reg_ictp_i_c_i,  freq_reg_ictp_i_c_i,  marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='Reg5-Holt')
plt.plot(x_freq_reg_ictp_ii_c_i, freq_reg_ictp_ii_c_i, marker='.', markersize=4, mfc='brown',  mec='brown',  linestyle='None', label='Reg5-UW')
plt.plot(x_freq_wrf_ncar_c_i,    freq_wrf_ncar_c_i,    marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='WRF-NCAR')
plt.plot(x_freq_wrf_ucan_c_i,    freq_wrf_ucan_c_i,    marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='WRF-UCAN')
plt.title('(a) Cluster I', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.legend(loc=1, nol=2, frameon=False)

ax = fig.add_subplot(3, 2, 2)
plt.plot(x_freq_inmet_smn_c_ii,   freq_inmet_smn_c_ii,   marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='INMET+SMN')
plt.plot(x_freq_era5_c_ii,        freq_era5_c_ii,        marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='ERA5')
plt.plot(x_freq_reg_usp_c_ii,     freq_reg_usp_c_ii,     marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='Reg4')
plt.plot(x_freq_reg_ictp_i_c_ii,  freq_reg_ictp_i_c_ii,  marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='Reg5-Holt')
plt.plot(x_freq_reg_ictp_ii_c_ii, freq_reg_ictp_ii_c_ii, marker='.', markersize=4, mfc='brown',  mec='brown',  linestyle='None', label='Reg5-UW')
plt.plot(x_freq_wrf_ncar_c_ii,    freq_wrf_ncar_c_ii,    marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='WRF-NCAR')
plt.plot(x_freq_wrf_ucan_c_ii,    freq_wrf_ucan_c_ii,    marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='WRF-UCAN')
plt.title('(b) Cluster II', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')

ax = fig.add_subplot(3, 2, 3)
plt.plot(x_freq_inmet_smn_c_iii,   freq_inmet_smn_c_iii,   marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='INMET+SMN')
plt.plot(x_freq_era5_c_iii,        freq_era5_c_iii,        marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='ERA5')
plt.plot(x_freq_reg_usp_c_iii,     freq_reg_usp_c_iii,     marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='Reg4')
plt.plot(x_freq_reg_ictp_i_c_iii,  freq_reg_ictp_i_c_iii,  marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='Reg5-Holt')
plt.plot(x_freq_reg_ictp_ii_c_iii, freq_reg_ictp_ii_c_iii, marker='.', markersize=4, mfc='brown',  mec='brown',  linestyle='None', label='Reg5-UW')
plt.plot(x_freq_wrf_ncar_c_iii,    freq_wrf_ncar_c_iii,    marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='WRF-NCAR')
plt.plot(x_freq_wrf_ucan_c_iii,    freq_wrf_ucan_c_iii,    marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='WRF-UCAN')
plt.title('(c) Cluster III', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.ylabel('Frequency (days)', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 2, 4)
plt.plot(x_freq_inmet_smn_c_iv,   freq_inmet_smn_c_iv,   marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='INMET+SMN')
plt.plot(x_freq_era5_c_iv,        freq_era5_c_iv,        marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='ERA5')
plt.plot(x_freq_reg_usp_c_iv,     freq_reg_usp_c_iv,     marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='Reg4')
plt.plot(x_freq_reg_ictp_i_c_iv,  freq_reg_ictp_i_c_iv,  marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='Reg5-Holt')
plt.plot(x_freq_reg_ictp_ii_c_iv, freq_reg_ictp_ii_c_iv, marker='.', markersize=4, mfc='brown',  mec='brown',  linestyle='None', label='Reg5-UW')
plt.plot(x_freq_wrf_ncar_c_iv,    freq_wrf_ncar_c_iv,    marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='WRF-NCAR')
plt.plot(x_freq_wrf_ucan_c_iv,    freq_wrf_ucan_c_iv,    marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='WRF-UCAN')
plt.title('(d) Cluster IV', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.ylabel('Frequency (days)', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 2, 5)
plt.plot(x_freq_inmet_smn_c_v,   freq_inmet_smn_c_v,   marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='INMET+SMN')
plt.plot(x_freq_era5_c_v,        freq_era5_c_v,        marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='ERA5')
plt.plot(x_freq_reg_usp_c_v,     freq_reg_usp_c_v,     marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='Reg4')
plt.plot(x_freq_reg_ictp_i_c_v,  freq_reg_ictp_i_c_v,  marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='Reg5-Holt')
plt.plot(x_freq_reg_ictp_ii_c_v, freq_reg_ictp_ii_c_v, marker='.', markersize=4, mfc='brown',  mec='brown',  linestyle='None', label='Reg5-UW')
plt.plot(x_freq_wrf_ncar_c_v,    freq_wrf_ncar_c_v,    marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='WRF-NCAR')
plt.plot(x_freq_wrf_ucan_c_v,    freq_wrf_ucan_c_v,    marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='WRF-UCAN')
plt.title('(e) Cluster V', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.xlabel('Intensity of daily precipitation (> 0 mm d⁻¹)', fontsize=font_size, fontweight='bold')

ax = fig.add_subplot(3, 2, 6)
plt.plot(x_freq_inmet_smn_c_vi,   freq_inmet_smn_c_vi,   marker='.', markersize=4, mfc='black',  mec='black',  linestyle='None', label='INMET+SMN')
plt.plot(x_freq_era5_c_vi,        freq_era5_c_vi,        marker='.', markersize=4, mfc='red',    mec='red',    linestyle='None', label='ERA5')
plt.plot(x_freq_reg_usp_c_vi,     freq_reg_usp_c_vi,     marker='.', markersize=4, mfc='blue',   mec='blue',   linestyle='None', label='Reg4')
plt.plot(x_freq_reg_ictp_i_c_vi,  freq_reg_ictp_i_c_vi,  marker='.', markersize=4, mfc='gray',   mec='gray',   linestyle='None', label='Reg5-Holt')
plt.plot(x_freq_reg_ictp_ii_c_vi, freq_reg_ictp_ii_c_vi, marker='.', markersize=4, mfc='brown',  mec='brown',  linestyle='None', label='Reg5-UW')
plt.plot(x_freq_wrf_ncar_c_vi,    freq_wrf_ncar_c_vi,    marker='.', markersize=4, mfc='green',  mec='green',  linestyle='None', label='WRF-NCAR')
plt.plot(x_freq_wrf_ucan_c_vi,    freq_wrf_ucan_c_vi,    marker='.', markersize=4, mfc='orange', mec='orange', linestyle='None', label='WRF-UCAN')
plt.title('(f) All clusters', loc='left', fontsize=font_size, fontweight='bold')
plt.yscale('log')
plt.xlabel('Intensity of daily precipitation (> 0 mm d⁻¹)', fontsize=font_size, fontweight='bold')

# Path out to save figure
path_out = '{0}/FPS_SESA/figs/paper_cp'.format(path)
name_out = 'pyplt_cdf_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
