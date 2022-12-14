# -*- coding: utf-8 -*-

__author__      = "Leidinice Silva"
__email__       = "leidinicesilva@gmail.com"
__date__        = "12/12/2022"
__description__ = "This script plot climatology maps from regcm5 and database"

import os
import cmocean
import netCDF4
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt

# mpl.use('Agg')

from mpl_toolkits.basemap import Basemap

	
def import_rcm(freq):
	
	path  = '/home/nice/Documentos/FPS_SESA/reg5'
	arq   = '{0}/pr_regcm5_sam4km-pbl1_SRF_{1}_2018-2019_lonlat.nc'.format(path, freq)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables['pr'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	rcm = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rcm


def import_rea(freq):
	
	path  = '/home/nice/Documentos/FPS_SESA/era5'
	arq   = '{0}/tp_era5_sesa_djf_2018-2019_lonlat.nc'.format(path, freq)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables['tp'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	rea = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, rea
	

def basemap(lat, lon):
	
	aux_lon1 = []
	aux_lon2 = []
	for l in lon:
		if l <= 180:
			aux_lon1.append(l)
		else:
			aux_lon2.append(l-360)
		
	lon = np.array(aux_lon1[::-1] + aux_lon2[::-1])
	new_lat = lat
	new_lon = lon[::-1]
	
	map = Basemap(projection='cyl', llcrnrlon=-75., llcrnrlat=-40., urcrnrlon=-40.,urcrnrlat=-15., resolution='c')
	lons, lats = np.meshgrid(new_lon, new_lat)
	xx, yy = map(lons,lats)
	
	path = '/home/nice/Documentos/github_projects/shp'
	map.readshapefile('{0}/shp_america_sul/america_sul'.format(path), 'world', drawbounds=True, color='black', linewidth=1.)
	
	return map, xx, yy
	

# Import regcm5 and obs database 
lat, lon, rcm_djf = import_rcm(u'djf')
lat, lon, rcm_mam = import_rcm(u'mam')
lat, lon, rcm_jja = import_rcm(u'jja')
lat, lon, rcm_son = import_rcm(u'son')

lat, lon, sat_djf = import_rea(u'djf')
lat, lon, sat_mam = import_rea(u'mam')
lat, lon, sat_jja = import_rea(u'jja')
lat, lon, sat_son = import_rea(u'son')

lat, lon, rea_djf = import_rea(u'djf')
lat, lon, rea_mam = import_rea(u'mam')
lat, lon, rea_jja = import_rea(u'jja')
lat, lon, rea_son = import_rea(u'son')

# Plot maps with the function
fig = plt.figure(figsize=(10, 6))
cmap = cmocean.cm.dense 
lev = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20]

ax = fig.add_subplot(3, 4, 1)
plt.title(u'(a) RegCM5-CP-4km DJF', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rcm_djf, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 2)
plt.title(u'(b) RegCM5-CP-4km MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rcm_mam, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 3)
plt.title(u'(c) RegCM5-CP-4km JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rcm_jja, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 4)
plt.title(u'(d) RegCM5-CP-4km SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rcm_son, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 5)
plt.title(u'(e) CMORPH DJF', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rea_djf, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 6)
plt.title(u'(f) CMORPH MAM', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rea_mam, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 7)
plt.title(u'(g) CMORPH JJA', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rea_jja, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 8)
plt.title(u'(h) CMORPH SON', loc='left', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rea_son, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 9)
plt.title(u'(i) ERA5 (4 km) DJF', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rea_djf, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 10)
plt.title(u'(j) ERA5 (4 km) MAM', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rea_mam, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 11)
plt.title(u'(k) ERA5 (4 km) JJA', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rea_jja, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 12)
plt.title(u'(l) ERA5 (4 km) SON', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, rea_son, levels=lev, latlon=True, cmap=cmap, extend='max') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

cb_ax = fig.add_axes([0.92, 0.2, 0.016, 0.6])
bounds = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 20]
cbar = fig.colorbar(plt_map, cax=cb_ax, orientation='vertical', boundaries=bounds, shrink=0.5, pad=0.5)
cbar.set_label(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
cbar.ax.tick_params(labelsize=8)  

# Path out to save figure
path_out = '/home/nice/Downloads'
name_out = 'pyplt_maps_clim_pr_regcm5_obs.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()



