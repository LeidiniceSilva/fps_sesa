# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot boxplot"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from pylab import setp
from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii


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
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.values
		mean_i.append(d_i*86400)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_day_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_day_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.values
		mean_v.append(d_v*86400)
		
		# Reading inmet 
		d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc/pre/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.resample(time='1D').sum()
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_day_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
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
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.values
		mean_i.append(d_i*86400)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_day_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_day_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.values
		mean_v.append(d_v*86400)
						
		# Reading smn 
		d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn_i/smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(smn_i[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.resample(time='1D').sum()
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_day_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
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
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.values
		mean_i.append(d_i*86400)
			
		# reading regcm ictp pbl 1 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl1_v0_day_20180601-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)
		
		# reading regcm ictp pbl 2
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_ictp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_ICTP-RegCM5pbl2_v0_day_20180601-20211231.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
						
		# reading wrf ncar 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_iv = d_iv.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.sel(lat=yy, lon=xx, method='nearest')
		d_iv = d_iv.values
		mean_iv.append(d_iv*86400)
			
		# reading wrf ucan 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_v = d_v.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.values
		mean_v.append(d_v*86400)
						
		# Reading smn 
		d_vi = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn_ii/smn_nc/' + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(smn_ii[i][0]))
		d_vi = d_vi.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_vi = d_vi.values
		mean_vi.append(d_vi)
		
		# reading era5 
		d_vii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_day_20180101-20211231.nc')
		d_vii = d_vii.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_vii = d_vii.sel(lat=yy, lon=xx, method='nearest')
		d_vii = d_vii.values
		mean_vii.append(d_vii)
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v, mean_vi, mean_vii
	

def setBoxColors(bp):
    setp(bp['boxes'][0], color='black')
    setp(bp['medians'][0], color='black')
    setp(bp['caps'][0], color='black')
    setp(bp['caps'][1], color='black')
    setp(bp['whiskers'][0], color='black')
    setp(bp['whiskers'][1], color='black')
        
    setp(bp['boxes'][1], color='gray')
    setp(bp['medians'][1], color='gray')
    setp(bp['caps'][2], color='gray')
    setp(bp['caps'][3], color='gray')
    setp(bp['whiskers'][2], color='gray')
    setp(bp['whiskers'][3], color='gray')
       
    setp(bp['boxes'][2], color='brown')
    setp(bp['medians'][2], color='brown')
    setp(bp['caps'][4], color='brown')
    setp(bp['caps'][5], color='brown')
    setp(bp['whiskers'][4], color='brown')
    setp(bp['whiskers'][5], color='brown')

    setp(bp['boxes'][3], color='green')
    setp(bp['medians'][3], color='green')
    setp(bp['caps'][6], color='green')
    setp(bp['caps'][7], color='green')
    setp(bp['whiskers'][6], color='green')
    setp(bp['whiskers'][7], color='green')

    setp(bp['boxes'][4], color='orange')
    setp(bp['medians'][4], color='orange')
    setp(bp['caps'][8], color='orange')
    setp(bp['caps'][9], color='orange')
    setp(bp['whiskers'][8], color='orange')
    setp(bp['whiskers'][9], color='orange')

    setp(bp['boxes'][5], color='blue')
    setp(bp['medians'][5], color='blue')
    setp(bp['caps'][10], color='blue')
    setp(bp['caps'][11], color='blue')
    setp(bp['whiskers'][10], color='blue')
    setp(bp['whiskers'][11], color='blue')

    setp(bp['boxes'][6], color='red')
    setp(bp['medians'][6], color='red')
    setp(bp['caps'][12], color='red')
    setp(bp['caps'][13], color='red')
    setp(bp['whiskers'][12], color='red')
    setp(bp['whiskers'][13], color='red')
    

var = 'pr'

# Import dataset
clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x, clim_vi_x, clim_vii_x = import_inmet()			
clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y, clim_vi_y, clim_vii_y = import_smn_i()			
clim_i_z, clim_ii_z, clim_iii_z, clim_iv_z, clim_v_z, clim_vi_z, clim_vii_z = import_smn_ii()			

reg_usp = clim_i_x + clim_i_y + clim_i_z
reg_ictp_i = clim_ii_x + clim_ii_y + clim_ii_z
reg_ictp_ii = clim_iii_x + clim_iii_y + clim_iii_z
wrf_ncar = clim_iv_x + clim_iv_y + clim_iv_z
wrf_ucan = clim_v_x + clim_v_y + clim_v_z
inmet_smn = clim_vi_x + clim_vi_y + clim_vi_z
era5 = clim_vii_x + clim_vi_y + clim_vii_z

list_hc = [4, 4, 4, 4, 4, 0, 0, 4, 4, 0, 0, 4, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0,
4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 3, 0, 0, 1, 0, 0,
0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 3, 0, 1, 1, 0, 1, 0, 1, 3, 3, 3, 3, 1, 1, 1, 3,
1, 1, 1, 1, 1, 3, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3,
3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
3, 1, 2, 4, 0, 4, 4, 4, 3, 4, 4, 4, 1, 4, 4, 2, 4, 4, 3, 1, 3, 3, 3, 3, 0, 3,
3, 3, 3, 1, 3, 3, 3, 3, 3, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 2, 0, 2, 2, 1,
3, 3, 1, 3, 0, 2, 1, 3, 3, 2, 2, 2, 2, 2, 2, 2, 0, 1, 2, 0, 2, 0, 1, 1, 2, 2,
1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 1, 3, 0, 1, 1, 2, 1, 2, 0, 2, 2, 2, 1, 2, 2,
2, 0, 4, 3, 4, 4, 3, 4, 4, 0, 4, 0, 4, 4, 4, 3, 4, 4, 1, 2, 2, 0, 1, 1, 0]

count_i = []
count_ii = []
count_iii = []
count_iv = []
count_v = []

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

reg_usp_i, reg_usp_ii, reg_usp_iii, reg_usp_iv, reg_usp_v = [], [], [], [], []
reg_ictp_i_i, reg_ictp_i_ii, reg_ictp_i_iii, reg_ictp_i_iv, reg_ictp_i_v = [], [], [], [], []
reg_ictp_ii_i, reg_ictp_ii_ii, reg_ictp_ii_iii, reg_ictp_ii_iv, reg_ictp_ii_v = [], [], [], [], []
wrf_ncar_i, wrf_ncar_ii, wrf_ncar_iii, wrf_ncar_iv, wrf_ncar_v = [], [], [], [], []
wrf_ucan_i, wrf_ucan_ii, wrf_ucan_iii, wrf_ucan_iv, wrf_ucan_v = [], [], [], [], []
inmet_smn_i, inmet_smn_ii, inmet_smn_iii, inmet_smn_iv, inmet_smn_v = [], [], [], [], []
era5_i, era5_ii, era5_iii, era5_iv, era5_v = [], [], [], [], []

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

reg_usp_c_i = np.nanmean(reg_usp_i, axis=0)
reg_ictp_i_c_i = np.nanmean(reg_ictp_i_i, axis=0)
reg_ictp_ii_c_i = np.nanmean(reg_ictp_ii_i, axis=0)
wrf_ncar_c_i = np.nanmean(wrf_ncar_i, axis=0)
wrf_ucan_c_i = np.nanmean(wrf_ucan_i, axis=0)
inmet_smn_c_i = np.nanmean(inmet_smn_i, axis=0)
era5_c_i = np.nanmean(era5_i, axis=0)

reg_usp_c_ii = np.nanmean(reg_usp_ii, axis=0)
reg_ictp_i_c_ii = np.nanmean(reg_ictp_i_ii, axis=0)
reg_ictp_ii_c_ii = np.nanmean(reg_ictp_ii_ii, axis=0)
wrf_ncar_c_ii = np.nanmean(wrf_ncar_ii, axis=0)
wrf_ucan_c_ii = np.nanmean(wrf_ucan_ii, axis=0)
inmet_smn_c_ii = np.nanmean(inmet_smn_ii, axis=0)
era5_c_ii = np.nanmean(era5_ii, axis=0)

reg_usp_c_iii = np.nanmean(reg_usp_iii, axis=0)
reg_ictp_i_c_iii = np.nanmean(reg_ictp_i_iii, axis=0)
reg_ictp_ii_c_iii = np.nanmean(reg_ictp_ii_iii, axis=0)
wrf_ncar_c_iii = np.nanmean(wrf_ncar_iii, axis=0)
wrf_ucan_c_iii = np.nanmean(wrf_ucan_iii, axis=0)
inmet_smn_c_iii = np.nanmean(inmet_smn_iii, axis=0)
era5_c_iii = np.nanmean(era5_iii, axis=0)

reg_usp_c_iv = np.nanmean(reg_usp_iv, axis=0)
reg_ictp_i_c_iv = np.nanmean(reg_ictp_i_iv, axis=0)
reg_ictp_ii_c_iv = np.nanmean(reg_ictp_ii_iv, axis=0)
wrf_ncar_c_iv = np.nanmean(wrf_ncar_iv, axis=0)
wrf_ucan_c_iv = np.nanmean(wrf_ucan_iv, axis=0)
inmet_smn_c_iv = np.nanmean(inmet_smn_iv, axis=0)
era5_c_iv = np.nanmean(era5_iv, axis=0)

reg_usp_c_v = np.nanmean(reg_usp_v, axis=0)
reg_ictp_i_c_v = np.nanmean(reg_ictp_i_v, axis=0)
reg_ictp_ii_c_v = np.nanmean(reg_ictp_ii_v, axis=0)
wrf_ncar_c_v = np.nanmean(wrf_ncar_v, axis=0)
wrf_ucan_c_v = np.nanmean(wrf_ucan_v, axis=0)
inmet_smn_c_v = np.nanmean(inmet_smn_v, axis=0)
era5_c_v = np.nanmean(era5_v, axis=0)

reg_usp_c_i = [i for i in reg_usp_c_i if i >= 0.1]	
reg_ictp_i_c_i = [i for i in reg_ictp_i_c_i if i >= 0.1]	
reg_ictp_ii_c_i = [i for i in reg_ictp_ii_c_i if i >= 0.1]	
wrf_ncar_c_i = [i for i in wrf_ncar_c_i if i >= 0.1]	
wrf_ucan_c_i = [i for i in wrf_ucan_c_i if i >= 0.1]	
inmet_smn_c_i = [i for i in inmet_smn_c_i if i >= 0.1]	
era5_c_i = [i for i in era5_c_i if i >= 0.1]	

reg_usp_c_ii = [i for i in reg_usp_c_ii if i >= 0.1]	
reg_ictp_i_c_ii = [i for i in reg_ictp_i_c_ii if i >= 0.1]	
reg_ictp_ii_c_ii = [i for i in reg_ictp_ii_c_ii if i >= 0.1]
wrf_ncar_c_ii = [i for i in wrf_ncar_c_ii if i >= 0.1]	
wrf_ucan_c_ii = [i for i in wrf_ucan_c_ii if i >= 0.1]	
inmet_smn_c_ii = [i for i in inmet_smn_c_ii if i >= 0.1]	
era5_c_ii = [i for i in era5_c_ii if i >= 0.1]	

reg_usp_c_iii = [i for i in reg_usp_c_iii if i >= 0.1]	
reg_ictp_i_c_iii = [i for i in reg_ictp_i_c_iii if i >= 0.1]	
reg_ictp_ii_c_iii = [i for i in reg_ictp_ii_c_iii if i >= 0.1]
wrf_ncar_c_iii = [i for i in wrf_ncar_c_iii if i >= 0.1]	
wrf_ucan_c_iii = [i for i in wrf_ucan_c_iii if i >= 0.1]	
inmet_smn_c_iii = [i for i in inmet_smn_c_iii if i >= 0.1]	
era5_c_iii = [i for i in era5_c_iii if i >= 0.1]	

reg_usp_c_iv = [i for i in reg_usp_c_iv if i >= 0.1]	
reg_ictp_i_c_iv = [i for i in reg_ictp_i_c_iv if i >= 0.1]	
reg_ictp_ii_c_iv = [i for i in reg_ictp_ii_c_iv if i >= 0.1]
wrf_ncar_c_iv = [i for i in wrf_ncar_c_iv if i >= 0.1]	
wrf_ucan_c_iv = [i for i in wrf_ucan_c_iv if i >= 0.1]	
inmet_smn_c_iv = [i for i in inmet_smn_c_iv if i >= 0.1]	
era5_c_iv = [i for i in era5_c_iv if i >= 0.1]	

reg_usp_c_v = [i for i in reg_usp_c_v if i >= 0.1]	
reg_ictp_i_c_v = [i for i in reg_ictp_i_c_v if i >= 0.1]	
reg_ictp_ii_c_v = [i for i in reg_ictp_ii_c_v if i >= 0.1]
wrf_ncar_c_v = [i for i in wrf_ncar_c_v if i >= 0.1]	
wrf_ucan_c_v = [i for i in wrf_ucan_c_v if i >= 0.1]	
inmet_smn_c_v = [i for i in inmet_smn_c_v if i >= 0.1]	
era5_c_v = [i for i in era5_c_v if i >= 0.1]	

cluster_i = [reg_usp_c_i, reg_ictp_i_c_i, reg_ictp_ii_c_i, wrf_ncar_c_i, wrf_ucan_c_i, inmet_smn_c_i, era5_c_i]
cluster_ii = [reg_usp_c_ii, reg_ictp_i_c_ii, reg_ictp_ii_c_ii, wrf_ncar_c_ii, wrf_ucan_c_ii, inmet_smn_c_ii, era5_c_ii]
cluster_iii = [reg_usp_c_iii, reg_ictp_i_c_iii, reg_ictp_ii_c_iii, wrf_ncar_c_iii, wrf_ucan_c_iii, inmet_smn_c_iii, era5_c_iii]
cluster_iv = [reg_usp_c_iv, reg_ictp_i_c_iv, reg_ictp_ii_c_iv, wrf_ncar_c_iv, wrf_ucan_c_iv, inmet_smn_c_iv, era5_c_iv]
cluster_v = [reg_usp_c_v, reg_ictp_i_c_v, reg_ictp_ii_c_v, wrf_ncar_c_v, wrf_ucan_c_v, inmet_smn_c_v, era5_c_v]

# Plot figure
fig = plt.figure(figsize=(10, 4))
x = np.arange(0.5, 40 + 0.5)

ax = fig.add_subplot(1, 1, 1)
bp = plt.boxplot(cluster_i, positions=[1, 2, 3, 4, 5, 6, 7], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_ii, positions=[9, 10, 11, 12, 13, 14, 15], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_iii, positions=[17, 18, 19, 20, 21, 22, 23], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_iv, positions=[25, 26, 27, 28, 29, 30, 31], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_v, positions=[33, 34, 35, 36, 37, 38, 39], sym='.')
setBoxColors(bp)

plt.xlim(0, 20)
plt.ylim(0, 80)
plt.yticks(np.arange(0, 88, 8))
plt.xticks(x, ('','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','','', ''))
plt.xlabel('Dataset', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.axhspan(72, 80, color='gray', alpha=0.5, lw=0)

plt.text(2, 74, 'Cluster I', fontweight='bold')
plt.text(10, 74, 'Cluster II', fontweight='bold')
plt.text(18, 74, 'Cluster III', fontweight='bold')
plt.text(26, 74, 'Cluster IV', fontweight='bold')
plt.text(34, 74, 'Cluster V', fontweight='bold')
plt.axhline(72, linewidth=1., linestyle='-',  color='black')
plt.axvline(8, linewidth=1., linestyle='--',  color='black')
plt.axvline(16, linewidth=1., linestyle='--',  color='black')
plt.axvline(24, linewidth=1., linestyle='--',  color='black')
plt.axvline(32, linewidth=1., linestyle='--',  color='black')

c1, = plt.plot([1,1],'black')
c2, = plt.plot([1,1],'gray')
c3, = plt.plot([1,1],'brown')
c4, = plt.plot([1,1],'green')
c5, = plt.plot([1,1],'orange')
c6, = plt.plot([1,1],'blue')
c7, = plt.plot([1,1],'red')
plt.legend((c1, c2, c3, c4, c5, c6, c7),('RegCM4', 'RegCM5 Holtslag', 'RegCM5 UW-PBL', 'WRF415', 'WRF433', 'INMET+SMN', 'ERA5'), bbox_to_anchor=(0.5, 1.09), loc=9, ncol=7, fontsize=8, frameon=False)
c1.set_visible(False)
c2.set_visible(False)
c3.set_visible(False)
c4.set_visible(False)
c5.set_visible(False)
c6.set_visible(False)
c7.set_visible(False)

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_boxplot_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
