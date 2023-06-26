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
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.values
		mean_i.append(d_i*86400)
		
		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
				
		# Reading inmet 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc_sesa/pre/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(inmet[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.resample(time='1D').sum()
		d_iv = d_iv.values
		mean_iv.append(d_iv)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_day_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
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
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/reg_usp/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1-USP-RegCM471_v0_day_20180601_20211231.nc')
		d_i = d_i.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_i = d_i.sel(lat=yy, lon=xx, method='nearest')
		d_i = d_i.values
		mean_i.append(d_i*86400)
		
		# reading wrf ncar 
		d_ii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ncar/' + 'pr_CSAM-4i_ERA5_evaluation_r1i1p1_NCAR-WRF415_v1_day_20180101-20211231.nc')
		d_ii = d_ii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_ii = d_ii.sel(lat=yy, lon=xx, method='nearest')
		d_ii = d_ii.values
		mean_ii.append(d_ii*86400)

		# reading wrf ucan 
		d_iii = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/rcm/wrf_ucan/' + 'pr_CSAM-4i_ECMWF-ERA5_evaluation_r1i1p1f1_UCAN-WRF433_v1_day_20180601-20210531.nc')
		d_iii = d_iii.pr.sel(time=slice('2018-06-01','2021-05-31'))
		d_iii = d_iii.sel(lat=yy, lon=xx, method='nearest')
		d_iii = d_iii.values
		mean_iii.append(d_iii*86400)
				
		# Reading smn 
		d_iv = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn/urug_smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(urug_smn[i][0]))
		d_iv = d_iv.pre.sel(time=slice('2018-06-01','2021-05-31'))
		d_iv = d_iv.resample(time='1D').sum()
		d_iv = d_iv.values
		mean_iv.append(d_iv)
		
		# reading era5 
		d_v = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/era5/' + 'tp_era5_csam_4km_day_20180101-20211231.nc')
		d_v = d_v.tp.sel(time=slice('2018-06-01','2021-05-31'))
		d_v = d_v.sel(lat=yy, lon=xx, method='nearest')
		d_v = d_v.values
		mean_v.append(d_v)
				
	return mean_i, mean_ii, mean_iii, mean_iv, mean_v
	

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

reg_usp_c_i = [i for i in reg_usp_c_i if i >= 0.1]	
wrf_ncar_c_i = [i for i in wrf_ncar_c_i if i >= 0.1]	
wrf_ucan_c_i = [i for i in wrf_ucan_c_i if i >= 0.1]	
inmet_smn_c_i = [i for i in inmet_smn_c_i if i >= 0.1]	
era5_c_i = [i for i in era5_c_i if i >= 0.1]	

reg_usp_c_ii = [i for i in reg_usp_c_ii if i >= 0.1]	
wrf_ncar_c_ii = [i for i in wrf_ncar_c_ii if i >= 0.1]	
wrf_ucan_c_ii = [i for i in wrf_ucan_c_ii if i >= 0.1]	
inmet_smn_c_ii = [i for i in inmet_smn_c_ii if i >= 0.1]	
era5_c_ii = [i for i in era5_c_ii if i >= 0.1]	

reg_usp_c_iii = [i for i in reg_usp_c_iii if i >= 0.1]	
wrf_ncar_c_iii = [i for i in wrf_ncar_c_iii if i >= 0.1]	
wrf_ucan_c_iii = [i for i in wrf_ucan_c_iii if i >= 0.1]	
inmet_smn_c_iii = [i for i in inmet_smn_c_iii if i >= 0.1]	
era5_c_iii = [i for i in era5_c_iii if i >= 0.1]	

reg_usp_c_iv = [i for i in reg_usp_c_iv if i >= 0.1]	
wrf_ncar_c_iv = [i for i in wrf_ncar_c_iv if i >= 0.1]	
wrf_ucan_c_iv = [i for i in wrf_ucan_c_iv if i >= 0.1]	
inmet_smn_c_iv = [i for i in inmet_smn_c_iv if i >= 0.1]	
era5_c_iv = [i for i in era5_c_iv if i >= 0.1]	

reg_usp_c_v = [i for i in reg_usp_c_v if i >= 0.1]	
wrf_ncar_c_v = [i for i in wrf_ncar_c_v if i >= 0.1]	
wrf_ucan_c_v = [i for i in wrf_ucan_c_v if i >= 0.1]	
inmet_smn_c_v = [i for i in inmet_smn_c_v if i >= 0.1]	
era5_c_v = [i for i in era5_c_v if i >= 0.1]	

cluster_i = [reg_usp_c_i, wrf_ncar_c_i, wrf_ucan_c_i, inmet_smn_c_i, era5_c_i]
cluster_ii = [reg_usp_c_ii, wrf_ncar_c_ii, wrf_ucan_c_ii, inmet_smn_c_ii, era5_c_ii]
cluster_iii = [reg_usp_c_iii, wrf_ncar_c_iii, wrf_ucan_c_iii, inmet_smn_c_iii, era5_c_iii]
cluster_iv = [reg_usp_c_iv, wrf_ncar_c_iv, wrf_ucan_c_iv, inmet_smn_c_iv, era5_c_iv]
cluster_v = [reg_usp_c_v, wrf_ncar_c_v, wrf_ucan_c_v, inmet_smn_c_v, era5_c_v]

# Plot figure
fig = plt.figure(figsize=(10, 4))
x = np.arange(0.5, 30 + 0.5)

ax = fig.add_subplot(1, 1, 1)
bp = plt.boxplot(cluster_i, positions=[1, 2, 3, 4, 5], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_ii, positions=[7, 8, 9, 10, 11], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_iii, positions=[13, 14, 15, 16, 17], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_iv, positions=[19, 20, 21, 22, 23], sym='.')
setBoxColors(bp)
bp = plt.boxplot(cluster_v, positions=[25, 26, 27, 28, 29], sym='.')
setBoxColors(bp)

plt.xlim(0, 20)
plt.ylim(0, 80)
plt.yticks(np.arange(0, 88, 8))
plt.xticks(x, ('','','','','','','','','','','','','','','','','','','','','','','','','','','','','', ''))
plt.xlabel('Dataset', fontweight='bold')
plt.ylabel('Precipitation (mm d⁻¹)', fontweight='bold')
plt.axhspan(72, 80, color='gray', alpha=0.5, lw=0)

plt.text(1.5, 74, 'Cluster I', fontweight='bold')
plt.text(7.5, 74, 'Cluster II', fontweight='bold')
plt.text(13.5, 74, 'Cluster III', fontweight='bold')
plt.text(19.5, 74, 'Cluster IV', fontweight='bold')
plt.text(25.5, 74, 'Cluster V', fontweight='bold')
plt.axhline(72, linewidth=1., linestyle='-',  color='black')
plt.axvline(6., linewidth=1., linestyle='--',  color='black')
plt.axvline(12, linewidth=1., linestyle='--',  color='black')
plt.axvline(18, linewidth=1., linestyle='--',  color='black')
plt.axvline(24, linewidth=1., linestyle='--',  color='black')

c1, = plt.plot([1,1],'black')
c2, = plt.plot([1,1],'green')
c3, = plt.plot([1,1],'orange')
c4, = plt.plot([1,1],'blue')
c5, = plt.plot([1,1],'red')
plt.legend((c1, c2, c3, c4, c5),('RegCM USP', 'WRF NCAR', 'WRF UCAN', 'INMET+SMN', 'ERA5'), bbox_to_anchor=(0.5, 1.09), loc=9, ncol=5, frameon=False)
c1.set_visible(False)
c2.set_visible(False)
c3.set_visible(False)
c4.set_visible(False)
c5.set_visible(False)

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs/paper_cp'
name_out = 'pyplt_boxplot_{0}_sesa.png'.format(var)
plt.savefig(os.path.join(path_out, name_out), dpi=400, bbox_inches='tight')
plt.show()
exit()
