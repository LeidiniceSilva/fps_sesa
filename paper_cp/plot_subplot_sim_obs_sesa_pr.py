# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "Jun 16, 2023"
__description__ = "This script plot subplots"

import os
import numpy as np
import xarray as xr
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt

from pylab import setp
from dict_inmet_stations import inmet
from dict_smn_stations import urug_smn


def import_inmet():
	
	clim_i = []
	mon_i = []
	day_i = []
	
	clim_ii = []
	mon_ii = []
	day_ii = []

	clim_iii = []
	mon_iii = []
	day_iii = []
	
	clim_iv = []
	mon_iv = []
	day_iv = []
	
	clim_v = []
	mon_v = []
	day_v = []

	# Select lat and lon 
	for i in range(1, 101):
		yy=inmet[i][2]
		xx=inmet[i][3]
		
		print('Reading weather station:', i, inmet[i][0], inmet[i][1])
		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		
		clim_d_i = d_i.groupby('time.month').mean('time')
		clim_d_i = clim_d_i.values
		clim_i.append(clim_d_i*86400)

		mon_d_i = d_i.resample(time='1M').mean()
		mon_d_i = mon_d_i.values
		mon_i.append(mon_d_i*86400)

		d_i = d_i.values
		day_i.append(d_i*86400)
		
		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		
		clim_d_ii = d_ii.groupby('time.month').mean('time')
		clim_d_ii = clim_d_ii.values
		clim_ii.append(clim_d_ii*86400)

		mon_d_ii = d_ii.resample(time='1M').mean()
		mon_d_ii = mon_d_ii.values
		mon_ii.append(mon_d_ii*86400)

		d_ii = d_ii.values
		day_ii.append(d_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')

		clim_d_iii = d_iii.groupby('time.month').mean('time')
		clim_d_iii = clim_d_iii.values
		clim_iii.append(clim_d_iii*86400)

		mon_d_iii = d_iii.resample(time='1M').mean()
		mon_d_iii = mon_d_iii.values
		mon_iii.append(mon_d_iii*86400)

		d_iii = d_iii.values
		day_iii.append(d_i*86400)
				
		# Reading inmet 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc_sesa/pre/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))

		clim_d_iv = d_iv.groupby('time.month').mean('time')
		clim_d_iv = clim_d_iv.values
		clim_iv.append(clim_d_iv*24)

		mon_d_iv = d_iv.resample(time='1M').mean()
		mon_d_iv = mon_d_iv.values
		mon_iv.append(mon_d_iv*24)
		
		day_d_iv = d_iv.resample(time='1D').sum()
		day_d_iv = day_d_iv.values
		day_iv.append(day_d_iv)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_day_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')

		clim_d_v = d_v.groupby('time.month').mean('time')
		clim_d_v = clim_d_v.values
		clim_v.append(clim_d_v)

		mon_d_v = d_v.resample(time='1M').mean()
		mon_d_v = mon_d_v.values
		mon_v.append(mon_d_v)

		day_d_v = d_v.values
		day_v.append(day_d_v)
				
	return clim_i, mon_i, day_i, clim_ii, mon_ii, day_ii, clim_iii, mon_iii, day_iii, clim_iv, mon_iv, day_iv, clim_v, mon_v, day_v
	
	
def import_smn():
	
	clim_i = []
	mon_i = []
	day_i = []
	
	clim_ii = []
	mon_ii = []
	day_ii = []

	clim_iii = []
	mon_iii = []
	day_iii = []
	
	clim_iv = []
	mon_iv = []
	day_iv = []
	
	clim_v = []
	mon_v = []
	day_v = []
	
	# Select lat and lon 
	for i in range(1, 72):
		yy=urug_smn[i][1]
		xx=urug_smn[i][2]
		
		print('Reading weather station:', i, urug_smn[i][0])
		
		# reading regcm usp 
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		
		clim_d_i = d_i.groupby('time.month').mean('time')
		clim_d_i = clim_d_i.values
		clim_i.append(clim_d_i*86400)

		mon_d_i = d_i.resample(time='1M').mean()
		mon_d_i = mon_d_i.values
		mon_i.append(mon_d_i*86400)

		d_i = d_i.values
		day_i.append(d_i*86400)
		
		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		
		clim_d_ii = d_ii.groupby('time.month').mean('time')
		clim_d_ii = clim_d_ii.values
		clim_ii.append(clim_d_ii*86400)

		mon_d_ii = d_ii.resample(time='1M').mean()
		mon_d_ii = mon_d_ii.values
		mon_ii.append(mon_d_ii*86400)

		d_ii = d_ii.values
		day_ii.append(d_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')

		clim_d_iii = d_iii.groupby('time.month').mean('time')
		clim_d_iii = clim_d_iii.values
		clim_iii.append(clim_d_iii*86400)

		mon_d_iii = d_iii.resample(time='1M').mean()
		mon_d_iii = mon_d_iii.values
		mon_iii.append(mon_d_iii*86400)

		d_iii = d_iii.values
		day_iii.append(d_i*86400)
		
		# Reading smn 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn/urug_smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(urug_smn[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))

		clim_d_iv = d_iv.groupby('time.month').mean('time')
		clim_d_iv = clim_d_iv.values
		clim_iv.append(clim_d_iv*24)

		mon_d_iv = d_iv.resample(time='1M').mean()
		mon_d_iv = mon_d_iv.values
		mon_iv.append(mon_d_iv*24)
		
		day_d_iv = d_iv.resample(time='1D').sum()
		day_d_iv = day_d_iv.values
		day_iv.append(day_d_iv)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_day_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')

		clim_d_v = d_v.groupby('time.month').mean('time')
		clim_d_v = clim_d_v.values
		clim_v.append(clim_d_v)

		mon_d_v = d_v.resample(time='1M').mean()
		mon_d_v = mon_d_v.values
		mon_v.append(mon_d_v)

		day_d_v = d_v.values
		day_v.append(day_d_v)
				
	return clim_i, mon_i, day_i, clim_ii, mon_ii, day_ii, clim_iii, mon_iii, day_iii, clim_iv, mon_iv, day_iv, clim_v, mon_v, day_v


def setBoxColors(bp):
    setp(bp['boxes'][0], color='black')
    setp(bp['medians'][0], color='black')
    setp(bp['caps'][0], color='black')
    setp(bp['caps'][1], color='black')
    setp(bp['whiskers'][0], color='black')
    setp(bp['whiskers'][1], color='black')
        
    setp(bp['boxes'][1], color='green')
    setp(bp['medians'][1], color='green')
    setp(bp['caps'][2], color='green')
    setp(bp['caps'][3], color='green')
    setp(bp['whiskers'][2], color='green')
    setp(bp['whiskers'][3], color='green')
       
    setp(bp['boxes'][2], color='orange')
    setp(bp['medians'][2], color='orange')
    setp(bp['caps'][4], color='orange')
    setp(bp['caps'][5], color='orange')
    setp(bp['whiskers'][4], color='orange')
    setp(bp['whiskers'][5], color='orange')

    setp(bp['boxes'][3], color='blue')
    setp(bp['medians'][3], color='blue')
    setp(bp['caps'][6], color='blue')
    setp(bp['caps'][7], color='blue')
    setp(bp['whiskers'][6], color='blue')
    setp(bp['whiskers'][7], color='blue')

    setp(bp['boxes'][4], color='red')
    setp(bp['medians'][4], color='red')
    setp(bp['caps'][9], color='red')
    setp(bp['whiskers'][8], color='red')
    setp(bp['whiskers'][9], color='red')
    

def compute_ccdf(dataset):
	
	x = np.sort(dataset)	
	cdf = x.cumsum() / x.sum()
	ccdf = 1 - cdf
		
	return x, ccdf
	
	
var = 'pr'

# Import dataset
clim_i_x, mon_i_x, day_i_x, clim_ii_x, mon_ii_x, day_ii_x, clim_iii_x, mon_iii_x, day_iii_x, clim_iv_x, mon_iv_x, day_iv_x, clim_v_x, mon_v_x, day_v_x = import_inmet()			
clim_i_y, mon_i_y, day_i_y, clim_ii_y, mon_ii_y, day_ii_y, clim_iii_y, mon_iii_y, day_iii_y, clim_iv_y, mon_iv_y, day_iv_y, clim_v_y, mon_v_y, day_v_y = import_smn()			

clim_reg_usp = clim_i_x + clim_i_y
clim_wrf_ncar = clim_ii_x + clim_ii_y
clim_wrf_ucan = clim_iii_x + clim_iii_y
clim_inmet_smn = clim_iv_x + clim_iv_y
clim_era5 = clim_v_x + clim_v_y

mon_reg_usp = mon_i_x + mon_i_y
mon_wrf_ncar = mon_ii_x + mon_ii_y
mon_wrf_ucan = mon_iii_x + mon_iii_y
mon_inmet_smn = mon_iv_x + mon_iv_y
mon_era5 = mon_v_x + mon_v_y

day_reg_usp = day_i_x + day_i_y
day_wrf_ncar = day_ii_x + day_ii_y
day_wrf_ucan = day_iii_x + day_iii_y
day_inmet_smn = day_iv_x + day_iv_y
day_era5 = day_v_x + day_v_y

clim_reg_usp = np.nanmean(clim_reg_usp, axis=0)
clim_wrf_ncar = np.nanmean(clim_wrf_ncar, axis=0)
clim_wrf_ucan = np.nanmean(clim_wrf_ucan, axis=0)
clim_inmet_smn = np.nanmean(clim_inmet_smn, axis=0)
clim_era5 = np.nanmean(clim_era5, axis=0)

mon_reg_usp = np.nanmean(mon_reg_usp, axis=0)
mon_wrf_ncar = np.nanmean(mon_wrf_ncar, axis=0)
mon_wrf_ucan = np.nanmean(mon_wrf_ucan, axis=0)
mon_inmet_smn = np.nanmean(mon_inmet_smn, axis=0)
mon_era5 = np.nanmean(mon_era5, axis=0)

day_reg_usp = np.nanmean(day_reg_usp, axis=0)
day_wrf_ncar = np.nanmean(day_wrf_ncar, axis=0)
day_wrf_ucan = np.nanmean(day_wrf_ucan, axis=0)
day_inmet_smn = np.nanmean(day_inmet_smn, axis=0)
day_era5 = np.nanmean(day_era5, axis=0)

day_reg_usp = [i for i in day_reg_usp if i >= 0.1]	
day_wrf_ncar = [i for i in day_wrf_ncar if i >= 0.1]	
day_wrf_ucan = [i for i in day_wrf_ucan if i >= 0.1]	
day_inmet_smn = [i for i in day_inmet_smn if i >= 0.1]	
day_era5 = [i for i in day_era5 if i >= 0.1]	

day_boxplot = [day_reg_usp, day_wrf_ncar, day_wrf_ucan, day_inmet_smn, day_era5]

x_reg_usp_ccdf, reg_usp_ccdf = compute_ccdf(day_reg_usp)
x_wrf_ncar_ccdf, wrf_ncar_ccdf = compute_ccdf(day_wrf_ncar)
x_wrf_ucan_ccdf, wrf_ucan_ccdf = compute_ccdf(day_wrf_ucan)
x_inmet_smn_ccdf, inmet_smn_ccdf = compute_ccdf(day_inmet_smn)
x_era5_ccdf, era5_ccdf = compute_ccdf(day_era5)

# Plot figure
fig = plt.figure(figsize=(14, 10))

ax = fig.add_subplot(2, 2, 1)
x = np.arange(1, 5 + 1)
bp = plt.boxplot(day_boxplot, positions=[1, 2, 3, 4, 5], sym='.')
setBoxColors(bp)
plt.title('(a) Daily frequency', loc='left', fontweight='bold')
plt.ylim(0, 80)
plt.yticks(np.arange(0, 88, 8))
plt.xticks(x, ('RegCM USP','WRF NCAR','WRF UCAN','INMET+SMN','ERA5'))
plt.xlabel('Dataset', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')

ax = fig.add_subplot(2, 2, 2)
y = np.arange(0, 1.25, 0.25)
plt.plot(x_reg_usp_ccdf, reg_usp_ccdf, marker='.', markersize=4, markerfacecolor='black', markeredgecolor='black', linestyle='None', label='RegCM USP')
plt.plot(x_wrf_ncar_ccdf, wrf_ncar_ccdf, marker='.', markersize=4, markerfacecolor='green', markeredgecolor='green', linestyle='None', label='WRF NCAR')
plt.plot(x_wrf_ucan_ccdf, wrf_ucan_ccdf, marker='.', markersize=4, markerfacecolor='orange', markeredgecolor='orange', linestyle='None', label='WRF UCAN')
plt.plot(x_inmet_smn_ccdf, inmet_smn_ccdf, marker='.', markersize=4, markerfacecolor='blue', markeredgecolor='blue', linestyle='None', label='INMET+SMN')
plt.plot(x_era5_ccdf, era5_ccdf, marker='.', markersize=4, markerfacecolor='red', markeredgecolor='red', linestyle='None', label='ERA5')
plt.title('(b) Daily boxplot', loc='left', fontweight='bold')
plt.yticks(y, ('10⁻⁸','10⁻⁶','10⁻⁴','10⁻²','10⁰'))
plt.ylabel('Frequency', fontweight='bold')
plt.xlabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.legend(frameon=False)

ax = fig.add_subplot(2, 2, 3)
time = np.arange(0.5, 12 + 0.5)
plt.plot(time, clim_reg_usp, linewidth=1.5, linestyle='--', markersize=3, marker='o', markerfacecolor='white', color='black', label = 'RegCM USP')
plt.plot(time, clim_wrf_ncar, linewidth=1.5, linestyle='--', markersize=3, marker='*', markerfacecolor='white', color='green', label = 'WRF NCAR')
plt.plot(time, clim_wrf_ucan, linewidth=1.5, linestyle='--', markersize=3, marker='x', markerfacecolor='white', color='orange', label = 'WRF UCAN')
plt.plot(time, clim_inmet_smn, linewidth=1.5, linestyle='--', markersize=3, marker='s', markerfacecolor='white', color='blue', label = 'INMET+SMN')
plt.plot(time, clim_era5, linewidth=1.5, linestyle='--', markersize=3, marker='^', markerfacecolor='white', color='red', label = 'ERA5')
plt.title('(c) Anual cycle', loc='left', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.ylim(0, 8)
plt.yticks(np.arange(0, 9, 1))
plt.xticks(time, ('Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'))
plt.xlabel('Months', fontweight='bold')
plt.legend(loc=9, ncol=3, fontsize=8, frameon=False)

ax = fig.add_subplot(2, 2, 4)
dt = pd.date_range(start="20180601", end="20210531", freq="M")
reg_usp_dt = pd.Series(data=mon_reg_usp, index=dt)
wrf_ncar_dt = pd.Series(data=mon_wrf_ncar, index=dt)
wrf_ucan_dt = pd.Series(data=mon_wrf_ucan, index=dt)
inmet_smn_dt = pd.Series(data=mon_inmet_smn, index=dt)
era5_dt = pd.Series(data=mon_era5, index=dt)

plt.plot(reg_usp_dt, linewidth=1, linestyle='--', color='black', label = 'RegCM USP')
plt.plot(wrf_ncar_dt, linewidth=1, linestyle='--', color='green', label = 'WRF NCAR')
plt.plot(wrf_ucan_dt, linewidth=1, linestyle='--', color='orange', label = 'WRF UCAN')
plt.plot(inmet_smn_dt, linewidth=1, linestyle='--', color='blue', label = 'INMET+SMN')
plt.plot(era5_dt, linewidth=1, linestyle='--', color='red', label = 'ERA5')
plt.title('(d) Monthly time series', loc='left', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.xlabel('Periods', fontweight='bold')
plt.ylim(0, 8)
plt.yticks(np.arange(0, 9, 1))
 
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_subplot_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
