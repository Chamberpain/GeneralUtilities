from __future__ import print_function
from netCDF4 import Dataset
import numpy as np 
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from GeneralUtilities.Compute.list import find_nearest
from GeneralUtilities.Compute.constants import degree_dist
from scipy.interpolate import griddata
from GeneralUtilities.Filepath.instance import get_base_folder

class DepthBase(object):
		
	def get_index_from_pos(self,pos):
		try: # this is for pos in LatLon form
			dummy_x = pos.lon.decimal_degree
			dummy_y = pos.lat.decimal_degree
		except AttributeError:	#this is for pos in tuple form
			pass
		if pos.__class__  in [list,tuple]:
			dummy_x = pos[1]
			dummy_y = pos[0]
		assert dummy_x<=180
		assert dummy_x>=-180
		assert dummy_y<=90
		assert dummy_y>=-90
		nearest_x = find_nearest(self.x,dummy_x)
		nearest_y = find_nearest(self.y,dummy_y)
		x_index = self.x.tolist().index(nearest_x)
		y_index = self.y.tolist().index(nearest_y)
		return (x_index,y_index)	

	def calculate_subsample(self,array):
		array_mask = ~np.isnan(array)
		X,Y = np.meshgrid(self.x,self.y)
		lons = np.arange(-180,180,0.1)
		lats = np.arange(-90,90,0.1)
		XX,YY = np.meshgrid(lons,lats)
		return (lons,lats,griddata(np.array(zip(X[array_mask],Y[array_mask])),array[array_mask],(XX,YY)))

	def return_z(self,pos):
		x_index,y_index = self.get_index_from_pos(pos)
		print('x index is')
		print(x_index)
		print('y index is')
		print(y_index)
		return self.z[y_index,x_index]

	def return_gradient(self,pos):
		x_index,y_index = self.get_index_from_pos(pos)
		dz_dx = self.dz_dx/(np.cos(np.deg2rad(self.y[y_index])))
		return (dz_dx[y_index,x_index],self.dz_dy[y_index,x_index])

	def guassian_smooth(self,sigma=5):
		self.z = gaussian_filter(self.z, sigma=sigma)

	def plot(self,coords=None):
#need to switch to cartopy
		if coords:
			lllon,urlon,lllat,urlat = [find_nearest(_[0],_[1]) for _ in zip([self.x,self.x,self.y,self.y],coords)]
			lllon_index,urlon_index,lllat_index,urlat_index = [_[0].tolist().index(_[1]) for _ in zip([self.x,self.x,self.y,self.y],[lllon,urlon,lllat,urlat])]
			x = self.x[lllon_index:urlon_index]
			y = self.y[urlat_index:lllat_index]
			z = self.z[urlat_index:lllat_index,lllon_index:urlon_index]
			XX,YY = np.meshgrid(x,y)
			m = Basemap(projection='cea',llcrnrlat=lllat,urcrnrlat=urlat,\
			llcrnrlon=lllon,urcrnrlon=urlon,resolution='l',lon_0=0,\
			fix_aspect=False)			
		else:
			XX,YY = np.meshgrid(self.x[::10],self.y[::10])
			z = self.z[::10,::10]
			m = Basemap(projection='cea',llcrnrlat=-90,urcrnrlat=90,\
			llcrnrlon=-180,urcrnrlon=180,resolution='l',lon_0=0,\
			fix_aspect=False)
		z = np.ma.masked_greater(z,-1000)
		m.pcolormesh(XX,YY,z,latlon=True)
		plt.colorbar()
		plt.show()

class ETopo1Depth(DepthBase):
	def __init__(self,stride=6):
		base_folder = get_base_folder()
		file_path = base_folder+'Raw/Bathymetry/ETOPO1_Bed_c_gdal.grd'
		self.nc_fid = Dataset(file_path)
		x_start,x_end = self.nc_fid['x_range'][:]
		y_start,y_end = self.nc_fid['y_range'][:]
		delta_x,delta_y = self.nc_fid['spacing'][:]
		self.x = np.arange(x_start,x_end,delta_x)
		self.y = np.arange(y_start,y_end,delta_y)[::-1]
		self.z = self.nc_fid['z'][:].reshape(len(self.y),len(self.x))
		self.x = self.x[::stride]
		self.y = self.y[::stride]
		self.z = self.z[::stride,::stride]
		self.dz_dy,self.dz_dx = np.gradient(self.z.data,self.y,self.x)


class PACIOOS(DepthBase):
	def __init__(self):
		file_path = '/Users/pchamberlain/Data/Raw/Bathymetry/hmrg_bathytopo_1km_mhi.nc'
		self.nc_fid = Dataset(file_path)
		self.x = self.nc_fid['x'][:]
		self.y = self.nc_fid['y'][:]
		self.z = self.nc_fid['z'][:]

	def too_shallow(self,pos):
		depth = self.return_z(pos)
		return depth>-500