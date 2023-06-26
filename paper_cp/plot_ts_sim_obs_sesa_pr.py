# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot annual cycle"

import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_stations import urug_smn


def import_inmet():
	
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []
	mean_v = []

	# Select lat and lon 
	for i in range(1, 101):
		yy=inmet[i][2]
		xx=inmet[i][3]
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.resample(time='1M').mean()
		d_i = d_i.values
		mean_i.append(d_i*86400)
		
		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.resample(time='1M').mean()
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.resample(time='1M').mean()
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
				
		# Reading inmet 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc_sesa/pre/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.resample(time='1M').mean()
		d_iv = d_iv.values
		mean_iv.append(d_iv*24)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.resample(time='1M').mean()
		d_v = d_v.values
		mean_v.append(d_v)
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v


def import_smn():
	
	mean_i = []
	mean_ii = []
	mean_iii = []
	mean_iv = []
	mean_v = []
	
	# Select lat and lon 
	for i in range(1, 72):
		yy=urug_smn[i][1]
		xx=urug_smn[i][2]
		
		print('Reading weather station:', i, urug_smn[i][0])	
		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_mon_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.resample(time='1M').mean()
		d_i = d_i.values
		mean_i.append(d_i*86400)
		
		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_mon_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.resample(time='1M').mean()
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_mon_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.resample(time='1M').mean()
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
				
		# Reading smn 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn/urug_smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(urug_smn[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.resample(time='1M').mean()
		d_iv = d_iv.values
		mean_iv.append(d_iv*24)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_mon_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.resample(time='1M').mean()
		d_v = d_v.values
		mean_v.append(d_v)
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v


var = 'pr'

# Import dataset
clim_i_x, clim_ii_x, clim_iii_x, clim_iv_x, clim_v_x = import_inmet()			
clim_i_y, clim_ii_y, clim_iii_y, clim_iv_y, clim_v_y = import_smn()			

reg_usp = clim_i_x + clim_i_y
wrf_ncar = clim_ii_x + clim_ii_y
wrf_ucan = clim_iii_x + clim_iii_y
inmet_smn = clim_iv_x + clim_iv_y
era5 = clim_v_x + clim_v_y

list_hc = [3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3,
3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 0, 2, 3, 3, 2, 0,
3, 0, 0, 3, 3, 0, 2, 0, 0, 3, 3, 2, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 1, 1,
1, 1, 1, 0, 1, 2, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 
1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 3, 2, 2, 2, 2, 2,
2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 4, 2, 4, 2, 4, 2, 4, 4, 4, 4,
4, 4, 2, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 2, 2, 2, 2, 2, 3, 2,
3, 2, 2, 2, 2, 1, 2]

print(len(reg_usp))
print(len(wrf_ncar))
print(len(wrf_ucan))
print(len(inmet_smn))
print(len(era5))
print(len(list_hc))

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

reg_usp_i = []
reg_usp_ii = []
reg_usp_iii = []
reg_usp_iv = []
reg_usp_v = []

wrf_ncar_i = []
wrf_ncar_ii = []
wrf_ncar_iii = []
wrf_ncar_iv = []
wrf_ncar_v = []

wrf_ucan_i = []
wrf_ucan_ii = []
wrf_ucan_iii = []
wrf_ucan_iv = []
wrf_ucan_v = []

inmet_smn_i = []
inmet_smn_ii = []
inmet_smn_iii = []
inmet_smn_iv = []
inmet_smn_v = []

era5_i = []
era5_ii = []
era5_iii = []
era5_iv = []
era5_v = []

for c_i in count_i:
	reg_usp_i.append(reg_usp[c_i])
	wrf_ncar_i.append(wrf_ncar[c_i])
	wrf_ucan_i.append(wrf_ucan[c_i])
	inmet_smn_i.append(inmet_smn[c_i])
	era5_i.append(era5[c_i])

for c_ii in count_ii:
	reg_usp_ii.append(reg_usp[c_ii])
	wrf_ncar_ii.append(wrf_ncar[c_ii])
	wrf_ucan_ii.append(wrf_ucan[c_ii])
	inmet_smn_ii.append(inmet_smn[c_ii])
	era5_ii.append(era5[c_ii])
	
for c_iii in count_iii:
	reg_usp_iii.append(reg_usp[c_iii])
	wrf_ncar_iii.append(wrf_ncar[c_iii])
	wrf_ucan_iii.append(wrf_ucan[c_iii])
	inmet_smn_iii.append(inmet_smn[c_iii])
	era5_iii.append(era5[c_iii])
	
for c_iv in count_iv:
	reg_usp_iv.append(reg_usp[c_iv])
	wrf_ncar_iv.append(wrf_ncar[c_iv])
	wrf_ucan_iv.append(wrf_ucan[c_iv])
	inmet_smn_iv.append(inmet_smn[c_iv])
	era5_iv.append(era5[c_iv])
	
for c_v in count_v:
	reg_usp_v.append(reg_usp[c_v])
	wrf_ncar_v.append(wrf_ncar[c_v])
	wrf_ucan_v.append(wrf_ucan[c_v])
	inmet_smn_v.append(inmet_smn[c_v])
	era5_v.append(era5[c_v])

reg_usp_c_i = np.nanmean(reg_usp_i, axis=0)
wrf_ncar_c_i = np.nanmean(wrf_ncar_i, axis=0)
wrf_ucan_c_i = np.nanmean(wrf_ucan_i, axis=0)
inmet_smn_c_i = np.nanmean(inmet_smn_i, axis=0)
era5_c_i = np.nanmean(era5_i, axis=0)

reg_usp_c_ii = np.nanmean(reg_usp_ii, axis=0)
wrf_ncar_c_ii = np.nanmean(wrf_ncar_ii, axis=0)
wrf_ucan_c_ii = np.nanmean(wrf_ucan_ii, axis=0)
inmet_smn_c_ii = np.nanmean(inmet_smn_ii, axis=0)
era5_c_ii = np.nanmean(era5_ii, axis=0)

reg_usp_c_iii = np.nanmean(reg_usp_iii, axis=0)
wrf_ncar_c_iii = np.nanmean(wrf_ncar_iii, axis=0)
wrf_ucan_c_iii = np.nanmean(wrf_ucan_iii, axis=0)
inmet_smn_c_iii = np.nanmean(inmet_smn_iii, axis=0)
era5_c_iii = np.nanmean(era5_iii, axis=0)

reg_usp_c_iv = np.nanmean(reg_usp_iv, axis=0)
wrf_ncar_c_iv = np.nanmean(wrf_ncar_iv, axis=0)
wrf_ucan_c_iv = np.nanmean(wrf_ucan_iv, axis=0)
inmet_smn_c_iv = np.nanmean(inmet_smn_iv, axis=0)
era5_c_iv = np.nanmean(era5_iv, axis=0)

reg_usp_c_v = np.nanmean(reg_usp_v, axis=0)
wrf_ncar_c_v = np.nanmean(wrf_ncar_v, axis=0)
wrf_ucan_c_v = np.nanmean(wrf_ucan_v, axis=0)
inmet_smn_c_v = np.nanmean(inmet_smn_v, axis=0)
era5_c_v = np.nanmean(era5_v, axis=0)

# Plot figure
fig = plt.figure(figsize=(12, 8))
dt = pd.date_range(start="20180601", end="20210531", freq="M")

ax = fig.add_subplot(3, 2, 1)
reg_usp_c_i_dt = pd.Series(data=reg_usp_c_i, index=dt)
wrf_ncar_c_i_dt = pd.Series(data=wrf_ncar_c_i, index=dt)
wrf_ucan_c_i_dt = pd.Series(data=wrf_ucan_c_i, index=dt)
inmet_smn_c_i_dt = pd.Series(data=inmet_smn_c_i, index=dt)
era5_c_i_dt = pd.Series(data=era5_c_i, index=dt)

plt.plot(reg_usp_c_i_dt, linewidth=1, linestyle='--', color='black', label = 'RegCM USP')
plt.plot(wrf_ncar_c_i_dt, linewidth=1, linestyle='--', color='green', label = 'WRF NCAR')
plt.plot(wrf_ucan_c_i_dt, linewidth=1, linestyle='--', color='orange', label = 'WRF UCAN')
plt.plot(inmet_smn_c_i_dt, linewidth=1, linestyle='--', color='blue', label = 'INMET+SMN')
plt.plot(era5_c_i_dt, linewidth=1, linestyle='--', color='red', label = 'ERA5')
plt.title('(a) Cluster I', loc='left', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1))
plt.setp(ax.get_xticklabels(), visible=False)
plt.legend(loc=9, ncol=3, fontsize=8, frameon=False)

ax = fig.add_subplot(3, 2, 2)
reg_usp_c_ii_dt = pd.Series(data=reg_usp_c_ii, index=dt)
wrf_ncar_c_ii_dt = pd.Series(data=wrf_ncar_c_ii, index=dt)
wrf_ucan_c_ii_dt = pd.Series(data=wrf_ucan_c_ii, index=dt)
inmet_smn_c_ii_dt = pd.Series(data=inmet_smn_c_ii, index=dt)
era5_c_ii_dt = pd.Series(data=era5_c_ii, index=dt)

plt.plot(reg_usp_c_ii_dt, linewidth=1, linestyle='--', color='black', label = 'RegCM USP')
plt.plot(wrf_ncar_c_ii_dt, linewidth=1, linestyle='--', color='green', label = 'WRF NCAR')
plt.plot(wrf_ucan_c_ii_dt, linewidth=1, linestyle='--', color='orange', label = 'WRF UCAN')
plt.plot(inmet_smn_c_ii_dt, linewidth=1, linestyle='--', color='blue', label = 'INMET+SMN')
plt.plot(era5_c_ii_dt, linewidth=1, linestyle='--', color='red', label = 'ERA5')
plt.title('(b) Cluster II', loc='left', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1))
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(3, 2, 3)
reg_usp_c_iii_dt = pd.Series(data=reg_usp_c_iii, index=dt)
wrf_ncar_c_iii_dt = pd.Series(data=wrf_ncar_c_iii, index=dt)
wrf_ucan_c_iii_dt = pd.Series(data=wrf_ucan_c_iii, index=dt)
inmet_smn_c_iii_dt = pd.Series(data=inmet_smn_c_iii, index=dt)
era5_c_iii_dt = pd.Series(data=era5_c_iii, index=dt)

plt.plot(reg_usp_c_iii_dt, linewidth=1, linestyle='--', color='black', label = 'RegCM USP')
plt.plot(wrf_ncar_c_iii_dt, linewidth=1, linestyle='--', color='green', label = 'WRF NCAR')
plt.plot(wrf_ucan_c_iii_dt, linewidth=1, linestyle='--', color='orange', label = 'WRF UCAN')
plt.plot(inmet_smn_c_iii_dt, linewidth=1, linestyle='--', color='blue', label = 'INMET+SMN')
plt.plot(era5_c_iii_dt, linewidth=1, linestyle='--', color='red', label = 'ERA5')
plt.title('(c) Cluster III', loc='left', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1))
plt.setp(ax.get_xticklabels(), visible=False)

ax = fig.add_subplot(3, 2, 4)
reg_usp_c_iv_dt = pd.Series(data=reg_usp_c_iv, index=dt)
wrf_ncar_c_iv_dt = pd.Series(data=wrf_ncar_c_iv, index=dt)
wrf_ucan_c_iv_dt = pd.Series(data=wrf_ucan_c_iv, index=dt)
inmet_smn_c_iv_dt = pd.Series(data=inmet_smn_c_iv, index=dt)
era5_c_iv_dt = pd.Series(data=era5_c_iv, index=dt)

plt.plot(reg_usp_c_iv_dt, linewidth=1, linestyle='--', color='black', label = 'RegCM USP')
plt.plot(wrf_ncar_c_iv_dt, linewidth=1, linestyle='--', color='green', label = 'WRF NCAR')
plt.plot(wrf_ucan_c_iv_dt, linewidth=1, linestyle='--', color='orange', label = 'WRF UCAN')
plt.plot(inmet_smn_c_iv_dt, linewidth=1, linestyle='--', color='blue', label = 'INMET+SMN')
plt.plot(era5_c_iv_dt, linewidth=1, linestyle='--', color='red', label = 'ERA5')
plt.title('(d) Cluster IV', loc='left', fontweight='bold')
plt.xlabel('Period', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1))
plt.xticks(fontsize=8)

ax = fig.add_subplot(3, 2, 5)
reg_usp_c_v_dt = pd.Series(data=reg_usp_c_v, index=dt)
wrf_ncar_c_v_dt = pd.Series(data=wrf_ncar_c_v, index=dt)
wrf_ucan_c_v_dt = pd.Series(data=wrf_ucan_c_v, index=dt)
inmet_smn_c_v_dt = pd.Series(data=inmet_smn_c_v, index=dt)
era5_c_v_dt = pd.Series(data=era5_c_v, index=dt)

plt.plot(reg_usp_c_v_dt, linewidth=1, linestyle='--', color='black', label = 'RegCM USP')
plt.plot(wrf_ncar_c_v_dt, linewidth=1, linestyle='--', color='green', label = 'WRF NCAR')
plt.plot(wrf_ucan_c_v_dt, linewidth=1, linestyle='--', color='orange', label = 'WRF UCAN')
plt.plot(inmet_smn_c_v_dt, linewidth=1, linestyle='--', color='blue', label = 'ERA5')
plt.plot(era5_c_v_dt, linewidth=1, linestyle='--', color='red', label = 'ERA5')
plt.title('(e) Cluster V', loc='left', fontweight='bold')
plt.xlabel('Period', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.ylim(0, 12)
plt.yticks(np.arange(0, 13, 1))
plt.xticks(fontsize=8)

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_ts_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
