# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "10/17/2022"
__description__ = "This script plot time series to each INMET automatic station"

import os
import numpy as np
import xarray as xr
import matplotlib.pyplot as plt

from dict_inmet_stations import inmet
from dict_smn_i_stations import smn_i
from dict_smn_ii_stations import smn_ii

stations = 'smn_ii' # 100, 72 and 87

# Getting the data, lat and lon 
for i in range(1, 100):
	
	if stations == 'inmet':
		yy = inmet[i][2]
		xx = inmet[i][3]
		name = inmet[i][0]
	elif stations == 'smn_i':
		yy = smn_i[i][1]
		xx = smn_i[i][2]
		name = smn_i[i][0]
	else:
		yy = smn_ii[i][1]
		xx = smn_ii[i][2]
		name = smn_ii[i][0]

	print('Reading weather station:', stations, name)
	if stations == 'inmet':
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/inmet/inmet_hr/inmet_nc/pre/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(name))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))	
		d_i = d_i.resample(time='1D').sum()
		d_i = d_i.values
	elif stations == 'smn_i':
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn_i/smn_nc/' + 'pre_{0}_H_2018-01-01_2021-12-31.nc'.format(name))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))	
		d_i = d_i.resample(time='1D').sum()
		d_i = d_i.values
	else:
		d_i = xr.open_dataset('/home/nice/Documentos/FPS_SESA/database/obs/smn_ii/smn_nc/' + 'pre_{0}_D_1979-01-01_2021-12-31.nc'.format(name))
		d_i = d_i.pre.sel(time=slice('2018-01-01','2021-12-31'))	
		d_i = d_i.values
		
	print('Plot figure')
	# Plot figure
	fig = plt.figure()
	time = np.arange(0.2, 1461 + 0.2)
	
	ax = fig.add_subplot(1, 1, 1)
	plt.plot(time, d_i, linewidth=1, color='blue')
	plt.title(u'{0} lat:{1} lon:{2}'.format(name, yy, xx), fontweight='bold')
	plt.xlabel(u'Period 2018-01-01 2021-12-31', fontweight='bold')
	plt.ylabel(u'Precipitation (mm d⁻¹)', fontweight='bold')
	plt.grid()
					
	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos/FPS_SESA/figs/ts'
	name_out = 'pyplt_ts_station_{0}_{1}_pr.png'.format(stations, name)
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

