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

	
def import_regcm(freq):
	
	path  = '/home/nice/Documentos/FPS_SESA/reg5'
	arq   = '{0}/pr_regcm5_sam4km-pbl1_SRF_{1}_2018-2019_lonlat.nc'.format(path, freq)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables['pr'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, mean


def import_cmorph(freq):
	
	path  = '/home/nice/Documentos/FPS_SESA/cmorph'
	arq   = '{0}/CMORPH_V1.0_ADJ_SESA_8km_{1}_2018-2019_lonlat.nc'.format(path, freq)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables['cmorph'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, mean
	
	
def import_era5(freq):
	
	path  = '/home/nice/Documentos/FPS_SESA/era5'
	arq   = '{0}/tp_era5_sesa_djf_2018-2019_lonlat.nc'.format(path, freq)	
	
	data  = netCDF4.Dataset(arq)
	var   = data.variables['tp'][:] 
	lat   = data.variables['lat'][:]
	lon   = data.variables['lon'][:]
	mean = np.nanmean(var[:][:,:,:], axis=0)

	return lat, lon, mean
	

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
lat, lon, regcm_djf = import_regcm(u'djf')
lat, lon, regcm_mam = import_regcm(u'mam')
lat, lon, regcm_jja = import_regcm(u'jja')
lat, lon, regcm_son = import_regcm(u'son')

lat, lon, cmorph_djf = import_cmorph(u'djf')
lat, lon, cmorph_mam = import_cmorph(u'mam')
lat, lon, cmorph_jja = import_cmorph(u'jja')
lat, lon, cmorph_son = import_cmorph(u'son')

lat, lon, era5_djf = import_era5(u'djf')
lat, lon, era5_mam = import_era5(u'mam')
lat, lon, era5_jja = import_era5(u'jja')
lat, lon, era5_son = import_era5(u'son')

diff_regcm_cmorph_djf = regcm_djf - cmorph_djf
diff_regcm_cmorph_mam = regcm_mam - cmorph_mam
diff_regcm_cmorph_jja = regcm_jja - cmorph_jja
diff_regcm_cmorph_son = regcm_son - cmorph_son

diff_regcm_era5_djf = regcm_djf - era5_djf
diff_regcm_era5_mam = regcm_mam - era5_mam
diff_regcm_era5_jja = regcm_jja - era5_jja
diff_regcm_era5_son = regcm_son - era5_son

# Plot maps with the function
fig = plt.figure(figsize=(10, 6))
cmap = cmocean.cm.tarn 
cmap = cmocean.tools.crop_by_percent(cmap, 30, which='both', N=None)
lev = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]

ax = fig.add_subplot(3, 4, 1)
plt.title(u'(a) RegCM5-CP_4km - CMORPH', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
plt.text(-47, -37, 'DJF', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_regcm_cmorph_djf, levels=lev, latlon=True, cmap=cmap, extend='both') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 2)
plt.title(u'(b) RegCM5-CP_4km - CMORPH', loc='left', fontsize=8, fontweight='bold')
plt.text(-47, -37, 'MAM', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_regcm_cmorph_mam, levels=lev, latlon=True, cmap=cmap, extend='both') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 3)
plt.title(u'(c) RegCM5-CP_4km - CMORPH', loc='left', fontsize=8, fontweight='bold')
plt.text(-47, -37, 'JJA', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_regcm_cmorph_jja, levels=lev, latlon=True, cmap=cmap, extend='both') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 4)
plt.title(u'(d) RegCM5-CP_4km - CMORPH', loc='left', fontsize=8, fontweight='bold')
plt.text(-47, -37, 'SON', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_regcm_cmorph_son, levels=lev, latlon=True, cmap=cmap, extend='both') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 5)
plt.title(u'(e) RegCM5-CP_4km - ERA5', loc='left', fontsize=8, fontweight='bold')
plt.ylabel(u'Latitude', fontsize=8, labelpad=20, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
plt.text(-47, -37, 'DJF', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_regcm_era5_djf, levels=lev, latlon=True, cmap=cmap, extend='both') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[1,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 6)
plt.title(u'(f) RegCM5-CP_4km - ERA5', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
plt.text(-47, -37, 'MAM', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_regcm_era5_mam, levels=lev, latlon=True, cmap=cmap, extend='both') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 7)
plt.title(u'(g) RegCM5-CP_4km - ERA5', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
plt.text(-47, -37, 'JJA', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_regcm_era5_jja, levels=lev, latlon=True, cmap=cmap, extend='both') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

ax = fig.add_subplot(3, 4, 8)
plt.title(u'(h) RegCM5-CP_4km - ERA5', loc='left', fontsize=8, fontweight='bold')
plt.xlabel(u'Longitude', fontsize=8, labelpad=20, fontweight='bold')
plt.text(-47, -37, 'SON', fontsize=8, fontweight='bold')
map, xx, yy = basemap(lat, lon)
plt_map = map.contourf(xx, yy, diff_regcm_era5_son, levels=lev, latlon=True, cmap=cmap, extend='both') 
map.drawmeridians(np.arange(-75.,-40.,5.), size=6, labels=[0,0,0,1], linewidth=0.5, color='black')
map.drawparallels(np.arange(-40.,-15.,5.), size=6, labels=[0,0,0,0], linewidth=0.5, color='black') 

cb_ax = fig.add_axes([0.92, 0.2, 0.016, 0.6])
bounds = lev = [-10, -8, -6, -4, -2, 2, 4, 6, 8, 10]
cbar = fig.colorbar(plt_map, cax=cb_ax, orientation='vertical', boundaries=bounds, shrink=0.5, pad=0.5)
cbar.set_label(u'Precipitation (mm d⁻¹)', fontsize=8, fontweight='bold')
cbar.ax.tick_params(labelsize=8)  

# Path out to save figure
path_out = '/home/nice/Documentos/FPS_SESA/figs'
name_out = 'pyplt_maps_bias_regcm5_obs.png'
if not os.path.exists(path_out):
	create_path(path_out)
plt.savefig(os.path.join(path_out, name_out), dpi=600, bbox_inches='tight')
plt.show()
exit()



