# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "10/17/2022"
__description__ = "This script plot time series to each Argentina automatic station"

import os
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

from dict_stations_arg_emas import arg_emas

dt = 'H_2018-01-01_2021-12-31'

for idx in range(1, 88):

	print()
	print('Reading Argentina weather station:', idx, arg_emas[idx][0])	
	# Reading inmet weather station	"pr"
	precip_input = '/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas_nc/precip_{1}_{2}.nc'.format(arg_emas[idx][0], dt)
	precip_input = nc.Dataset(precip_input)
	precip = precip_input['precip'][:]

	temp_input = '/home/nice/Documentos/FPS_SESA/arg_emas/arg_emas_nc/temp_{1}_{2}.nc'.format(arg_emas[idx][0], dt)
	temp_input = nc.Dataset(temp_input)
	temp = temp_input['temp'][:]
		
	print('Plot figure')
	# Plot figure
	fig = plt.figure(figsize=(10, 6))
	time = np.arange(0.2, 35064 + 0.2)
	
	ax = fig.add_subplot(2, 1, 1)
	plt.plot(time, precip, linewidth=0.5, color='blue', label = 'EMAS')
	plt.title(u'{0}'.format(arg_emas[idx][0]), fontsize=8, fontweight='bold')
	plt.xlabel(u'00h 01/01/2018 - 23h 31/12/2021', fontsize=8, fontweight='bold')
	plt.ylabel(u'Precipitation (mm h$\mathregular{^{-1}}$)', fontsize=8, fontweight='bold')
	plt.setp(ax.get_xticklabels(), visible=False)
	plt.yticks(np.arange(0, 22, 2), fontsize=8)
	plt.ylim(-0.5, 20)
	plt.legend(loc=1, fontsize=8)
	plt.grid()
	
	ax = fig.add_subplot(2, 1, 2)
	plt.plot(time, temp, linewidth=0.5, color='blue', label = 'EMAS')
	plt.title(u'{0}'.format(arg_emas[idx][0]), fontsize=8, fontweight='bold')
	plt.xlabel(u'00h 01/01/2018 - 23h 31/12/2021', fontsize=8, fontweight='bold')
	plt.ylabel(u'Temperature (°C)', fontsize=8, fontweight='bold')
	plt.yticks(np.arange(0, 33, 3), fontsize=8)
	plt.ylim(-0.5, 30)
	plt.legend(loc=1, fontsize=8)
	plt.grid()
		
	print('Path out to save figure')
	# Path out to save figure
	path_out = '/home/nice/Documentos/FPS_SESA/figs'
	name_out = 'pyplt_nc_vars_{0}.png'.format(arg_emas[idx][0])
	plt.savefig(os.path.join(path_out, name_out), dpi=100, bbox_inches='tight')
	plt.close('all')
	plt.cla()
exit()

