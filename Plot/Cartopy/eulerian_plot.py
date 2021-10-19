import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
import cartopy
from cartopy import geodesic
import os
import cartopy.mpl.geoaxes
import scipy
from matplotlib.colors import LinearSegmentedColormap
import pickle
from GeneralUtilities.Data.depth.depth_utilities import ETopo1Depth
from GeneralUtilities.Filepath.instance import get_base_folder


datadir = get_base_folder()+'Processed/eulerian_plot/basemap/'




def ellipse(self,geod,lon, lat, b, a, n_samples=360,phi=0):
	"""
	Return the coordinates of a geodetic ellipse of a given
	x radius of a and y radius of b around a lon/lat point.

	Radius is in meters in the geodetic's coordinate system.

	"""

	tsid = a*b
	radius = []
	phi_r = phi* np.pi / 180
	for az in (np.linspace(360, 0, n_samples)* np.pi / 180):
		A = a * np.sin(az)*np.cos(phi_r)+b*np.cos(az)*np.sin(phi_r)
		B = b * np.cos(az)*np.cos(phi_r)-b*np.sin(az)*np.sin(phi_r)
		r = tsid / (B**2. + A**2.)**0.5
		radius.append(r)

	lons, lats, back_azim = geod.fwd(np.repeat(lon, n_samples),
									 np.repeat(lat, n_samples),
									 np.linspace(360, 0, n_samples),
									 np.array(radius),
									 radians=False,
									 )
	return (lons,lats)

def streamline_plot(self):
	mat = scipy.io.loadmat(os.path.join(datadir,'agva_mean_streamfunction_1000dbar.mat'))
	XC = mat['long'][:][0]
	XC = np.array(XC.tolist())
	XC[XC>290] = XC[XC>290]-360
	YC = mat['latg'][:][0]
	streamline=mat['mean_streamfunction_1000dbar']

	XC,YC = np.meshgrid(XC,YC)
	levels = np.arange(np.nanmin(streamline),np.nanmax(streamline),10)
	self.contour(XC,YC,streamline,linewidths=2,colors='orange',animated=True,levels = levels, label='Streamlines',alpha=0.7)


def bathy(self,color=plt.cm.Oranges,contour=False):
	r_start = 0.0
	g_start = 0.5
	b_start = 0.5
	delta = 0.08
	cdict = {'red':  ((r_start, g_start, b_start),
			(r_start+0.2, g_start+delta, b_start+delta),
			(r_start+0.4, g_start+2*delta, b_start+2*delta),
			(r_start+0.6, g_start+3*delta, b_start+3*delta),
			(r_start+0.8, g_start+4*delta, b_start+4*delta),
			(r_start+1.0, g_start+5*delta, b_start+5*delta)),

	 'green':((r_start, g_start, b_start),
			(r_start+0.2, g_start+delta, b_start+delta),
			(r_start+0.4, g_start+2*delta, b_start+2*delta),
			(r_start+0.6, g_start+3*delta, b_start+3*delta),
			(r_start+0.8, g_start+4*delta, b_start+4*delta),
			(r_start+1.0, g_start+5*delta, b_start+5*delta)),

	 'blue': ((r_start, g_start, b_start),
			(r_start+0.2, g_start+delta, b_start+delta),
			(r_start+0.4, g_start+2*delta, b_start+2*delta),
			(r_start+0.6, g_start+3*delta, b_start+3*delta),
			(r_start+0.8, g_start+4*delta, b_start+4*delta),
			(r_start+1.0, g_start+5*delta, b_start+5*delta)),
	}
	bathy = LinearSegmentedColormap('bathy', cdict)
	depth = ETopo1Depth()
	plot_data = -depth.z/1000.
	XX,YY = np.meshgrid(depth.x,depth.y)
	levels = [0,1,2,3,4,5]
	if not contour:
		cf = self.contourf(XX,YY,plot_data,levels,cmap=bathy,animated=True,vmax=5,vmin=0)
		return cf
	else:
		self.contour(XX,YY,plot_data,cmap=plt.get_cmap('Greys'),alpha=.8,vmax=5,vmin=-2)


cartopy.mpl.geoaxes.GeoAxesSubplot.streamline_plot = streamline_plot
cartopy.mpl.geoaxes.GeoAxesSubplot.bathy = bathy
cartopy.mpl.geoaxes.GeoAxesSubplot.ellipse = ellipse


class BaseCartopy():
	def __init__(self,lat_grid=np.arange(-90,90.1),lon_grid=np.arange(-180,180),ax=False,projection=ccrs.PlateCarree()):
		assert max(lat_grid)<=90
		assert min(lat_grid)>=-90
		assert max(lon_grid)<=180
		assert min(lon_grid)>=-180

		self.lat_grid = lat_grid
		self.lon_grid = lon_grid

		if not ax:
			fig = plt.figure()
			self.ax = fig.add_subplot(1, 1, 1, projection=projection)
		else:
			self.ax = ax

	def meshgrid_xx_yy(self):
		return np.meshgrid(self.lon_grid,self.lat_grid)

	def finish_map(self):
		self.ax.add_feature(cfeature.LAND,zorder=10)
		self.ax.add_feature(cfeature.COASTLINE,zorder=10)
		self.ax.set_aspect('auto')
		gl = self.ax.gridlines(draw_labels=True)
		gl.xlabels_top = False
		gl.ylabels_right = False

	def get_map(self):
		XX,YY = self.meshgrid_xx_yy()
		return (XX,YY,self.ax)

class GlobalCartopy(BaseCartopy):
	def __init__(self,*args,**kwargs):
		super().__init__(*args,**kwargs)      
		print('I am plotting global region')
		llcrnrlon=-180.
		llcrnrlat=-80.
		urcrnrlon=180.
		urcrnrlat=80
		self.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat], crs=ccrs.PlateCarree())
		self.finish_map()

class HypernavCartopy(BaseCartopy):
	def __init__(self,nc,float_pos_dict,*args,pad=1,**kwargs):
		super().__init__(*args,**kwargs)     
		try:
			urlat = nc['lat'][:].max()
			lllat = nc['lat'][:].min()
			urlon = nc['lon'][:].max()
			lllon = nc['lon'][:].min()
		except TypeError:
			urlat = max(nc.lats)
			lllat = min(nc.lats)
			urlon = max(nc.lons)
			lllon = min(nc.lons)        
		llcrnrlon=(lllon-pad)
		llcrnrlat=(lllat-pad)
		urcrnrlon=(urlon+pad)
		urcrnrlat=(urlat+pad)
		self.ax.scatter(float_pos_dict['lon'],float_pos_dict['lat'],c='pink',linewidths=5,marker='x',s=80,zorder=10)
		self.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat], crs=ccrs.PlateCarree())
		self.finish_map()

class PointCartopy(BaseCartopy):
	def __init__(self,geopy_point,*args,pad=1,**kwargs):
		super().__init__(*args,**kwargs)          
		llcrnrlon=(geopy_point.longitude-pad)
		llcrnrlat=(geopy_point.latitude-pad)
		urcrnrlon=(geopy_point.longitude+pad)
		urcrnrlat=(geopy_point.latitude+pad)
		self.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat], crs=ccrs.PlateCarree())
		self.finish_map()




