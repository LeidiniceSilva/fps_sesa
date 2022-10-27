# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "10/17/2022"
__description__ = "This script plot time series to each inmet automatic station"

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

from dict_inmet_stations_code import codes
from dict_inmet_stations_latlon import coord
from dict_inmet_stations_name import names

dt = 'H_2018-01-01_2021-12-31'

# Getting the data 
for j in range(1, 289):
	
	if j == 4:
		continue
	if j == 12:
		continue
	if j == 45:
		continue
	if j == 55:
		continue
	if j == 77:
		continue
	if j == 98:
		continue
	if j == 99:
		continue
	if j == 118:
		continue
	if j == 122:
		continue
	if j == 130:
		continue
	if j == 135:
		continue
	if j == 151:
		continue
	if j == 155:
		continue
	if j == 159:
		continue
	if j == 160:
		continue
	if j == 163:
		continue
	if j == 164:
		continue
	if j == 181:
		continue
	if j == 183:
		continue
	if j == 186:
		continue
	if j == 187:
		continue
	if j == 188:
		continue
	if j == 209:
		continue
	if j == 216:
		continue
	if j == 228:
		continue
	if j == 236:
		continue
	if j == 268:
		continue
	if j == 287:
		continue
			
	yy = coord[j][0]
	xx = coord[j][1]
	name_i = names[j][0]
	name_ii = names[j][1]

	print('Reading inmet weather station:', j, codes[j], names[j][1])
	# Reading inmet weather station	"pr"
	pr_input = '/home/nice/Downloads/FPS_SESA/inmet/inmet_nc/pr_{0}_{1}.nc'.format(codes[j], dt)
	pr_input = nc.Dataset(pr_input)
	pre = pr_input['pre'][:]

	# Reading inmet weather station	"tp"
	tp_input = '/home/nice/Downloads/FPS_SESA/inmet/inmet_nc/tp_{0}_{1}.nc'.format(codes[j], dt)
	tp_input = nc.Dataset(tp_input)
	tmp = tp_input['tmp'][:]

	# Reading inmet weather station	"uv"
	uv_input = '/home/nice/Downloads/FPS_SESA/inmet/inmet_nc/uv_{0}_{1}.nc'.format(codes[j], dt)
	uv_input = nc.Dataset(uv_input)
	uv = uv_input['uv'][:]
	
	print('Plot figure')
	# Plot figure
	fig = plt.figure(figsize=(10, 6))
	time = np.arange(0.2, 35064 + 0.2)
	
	ax = fig.add_subplot(3, 1, 1)
	plt.plot(time, pre, linewidth=0.5, color='blue', label = 'INMET')
	plt.title('{0} lat: {1} lon: {2}'.format(name_ii, yy, xx), fontsize=8, fontweight='bold')
	plt.ylabel('Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.yticks(np.arange(0, 22, 2), fontsize=8)
	plt.ylim(-0.5, 20)
	plt.grid(linestyle='--')
	plt.legend(loc=1)
	
	ax = fig.add_subplot(3, 1, 2)
	plt.plot(time, tmp, linewidth=0.5, color='blue', label = 'INMET')
	plt.title('{0} lat: {1} lon: {2}'.format(name_ii, yy, xx), fontsize=8, fontweight='bold')
	plt.ylabel('Temperature (°C)', fontsize=8, fontweight='bold')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.yticks(np.arange(0, 33, 3), fontsize=8)
	plt.ylim(-0.5, 30)
	plt.grid(linestyle='--')
	plt.legend(loc=1)
	
	ax = fig.add_subplot(3, 1, 3)
	plt.plot(time, uv, linewidth=0.5, color='blue', label = 'INMET')
	plt.title('{0} lat: {1} lon: {2}'.format(name_ii, yy, xx), fontsize=8, fontweight='bold')
	plt.xlabel('00h 01/01/2018 - 23h 31/12/2021', fontsize=8, fontweight='bold')
	plt.ylabel('Wind speed (m s⁻¹)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 22, 2), fontsize=8)
	plt.ylim(-0.5, 20)
	plt.grid(linestyle='--')
	plt.legend(loc=1)
			
	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Downloads/FPS_SESA/figs'
	name_out = 'pyplt_nc_vars_{0}_{1}_{2}.png'.format(yy, xx, name_i)
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

