# -*- coding:utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "09/19/2022"
__description__ = "This script plot point for each inmet automatic station over sesa domain"

import os
import conda
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt

from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from dict_inmet_stations_code import codes
from dict_inmet_stations_latlon import coord
from dict_inmet_stations_name import names


def basemap(xx, yy):
	
	my_map = Basemap(projection='cyl', llcrnrlon=-75, llcrnrlat=-40., urcrnrlon=-35.,urcrnrlat=-10., resolution='c')
	my_map.drawmeridians(np.arange(-75.,-25.,10.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
	my_map.drawparallels(np.arange(-40.,5.,10.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 
	my_map.readshapefile('/home/nice/Documents/github_projects/shp/lim_pais/lim_pais', 'world', drawbounds=True, color='black', linewidth=0.5)
	x, y = my_map(xx, yy)
	
	return my_map, x, y
	
	
idx=1
dt = 'H_2018-01-01_2021-12-31'
xx = []		  
yy = []
bias = []
corr = []

dict_var = {1: ['pre', 'Pr_1h (mm d$\mathregular{^{-1}}$)', 'tp'],
            5: ['tmp', 'Tmp_1h (°C)', 't2m'],
            11: ['uv', 'Wind_1h (m s$\mathregular{^{-1}}$)', 'uv10']}

# Select lat and lon 
for i in range(1, 289):
	
	if i == 4:
		continue
	if i == 12:
		continue
	if i == 45:
		continue
	if i == 55:
		continue
	if i == 77:
		continue
	if i == 98:
		continue
	if i == 99:
		continue
	if i == 118:
		continue
	if i == 122:
		continue
	if i == 130:
		continue
	if i == 135:
		continue
	if i == 151:
		continue
	if i == 155:
		continue
	if i == 159:
		continue
	if i == 160:
		continue
	if i == 163:
		continue
	if i == 164:
		continue
	if i == 181:
		continue
	if i == 183:
		continue
	if i == 186:
		continue
	if i == 187:
		continue
	if i == 188:
		continue
	if i == 209:
		continue
	if i == 216:
		continue
	if i == 228:
		continue
	if i == 236:
		continue
	if i == 246:
		continue
	if i == 268:
		continue
	if i == 287:
		continue
	
	yy.append(coord[i][0])
	xx.append(coord[i][1])

	print('Reading inmet weather station:', i, codes[i], names[i][1])
	# Reading inmet weather station	
	df = pd.read_csv(os.path.join('/home/nice/Downloads/FPS_SESA/inmet/inmet_used/', 'dados_{0}_{1}.csv'.format(codes[i], dt)), sep='[:,|_]', engine='python')
	df['Data Medicao'] = pd.to_datetime(df['Data Medicao'], format='%Y-%m-%d %H:%M:%S', errors='ignore')
	df_i = df.groupby(pd.Grouper(key='Data Medicao', freq='M')).mean()

	if idx == 1:
		var_x = df_i.iloc[:,idx]
		clim = []
		for mon in range(0, 12):
			mon_x = np.nanmean(var_x[mon::12], axis=0)
			mon_xx = mon_x*24
			clim.append(mon_xx)
	else:
		var_x = df_i.iloc[:,idx]
		clim = []
		for mon in range(0, 12):
			mon_x = np.nanmean(var_x[mon::12], axis=0)
			clim.append(mon_x)
	
	# reading era5 reanalisis
	ds = xr.open_dataset('/home/nice/Documentos/FPS_SESA/era5/' + '{0}_sesa_era5_2018-2021.nc'.format(dict_var[idx][2]))

	if idx == 1:
		ds = ds.tp.sel(time=slice('2018-01-01','2021-12-31'))
		ds = ds.sel(latitude=coord[i][0], longitude=coord[i][1], method='nearest')
		var = ds.groupby('time.month').mean('time')
		var_i = var.values
		clim_i = var_i*24
	elif idx == 5:
		ds = ds.t2m.sel(time=slice('2018-01-01','2021-12-31'))
		ds = ds.sel(latitude=coord[i][0], longitude=coord[i][1], method='nearest')
		var = ds.groupby('time.month').mean('time')
		clim_i = var.values
	else:
		ds = ds.u10.sel(time=slice('2018-01-01','2021-12-31'))
		ds = ds.sel(latitude=coord[i][0], longitude=coord[i][1], method='nearest')
		var = ds.groupby('time.month').mean('time')
		clim_i = var.values
	
	corr_ = np.corrcoef(clim, clim_i)[0][1]
	corr.append(corr_)

	bias_ = np.nanmean(clim_i) - np.nanmean(clim)
	bias.append(bias_)

print('Plot figure')
# Plot figure   
fig = plt.figure()

if idx == 1:
	minb=-2
	maxb=2
	minc=-1
	maxc=1
	colorv='BrBG'
	
elif idx == 5:
	minb=-2
	maxb=2
	minc=-1
	maxc=1
	colorv='bwr'

else:
	minb=-2
	maxb=2
	minc=-1
	maxc=1
	colorv='PiYG'
	
ax = fig.add_subplot(1, 2, 1)
my_map, x, y = basemap(xx, yy)
pltfig = my_map.scatter(x, y, 5, bias, cmap=colorv, marker='o', vmin=minb, vmax=maxb)
plt.title('(a) Bias {0}'.format(dict_var[idx][1]), loc='left', fontsize=8)
plt.ylabel(u'Latitude', fontsize=6, labelpad=15)
plt.xlabel(u'Longitude', fontsize=6, labelpad=15)
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.91, 0.35, 0.019, 0.28]), extend='both')
cbar.ax.tick_params(labelsize=6)

ax = fig.add_subplot(1, 2, 2)
my_map, x, y = basemap(xx, yy)
pltfig = my_map.scatter(x, y, 5, corr, cmap=colorv, marker='o', vmin=minc, vmax=maxc)
plt.title('(b) Correlation {0}'.format(dict_var[idx][1]), loc='left', fontsize=8)
plt.xlabel(u'Longitude', fontsize=6, labelpad=15)
cbar = plt.colorbar(pltfig, cax=fig.add_axes([0.99, 0.35, 0.019, 0.28]), extend='both')
cbar.ax.tick_params(labelsize=6)

print('Path out to save figure')
# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_maps_stats_{0}_weather_station.png'.format(dict_var[idx][0])
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=300, bbox_inches='tight')
exit()

