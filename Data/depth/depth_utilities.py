from __future__ import print_function
from netCDF4 import Dataset
import numpy as np 
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from GeneralUtilities.Compute.list import LatList,LonList
from GeneralUtilities.Compute.constants import degree_dist
from scipy.interpolate import griddata
from GeneralUtilities.Filepath.instance import get_base_folder
import geopy

class DepthBase(object):
	def __init__(self,*args,lon=None,lat=None,z=None,dz_dy=None, dz_dx = None,**kwargs):
		assert isinstance(lon,LonList)
		assert isinstance(lat,LatList)
		assert z.shape[1] == len(lon)
		assert z.shape[0] == len(lat)
		assert isinstance(z,np.ma.masked_array)
		self.z = np.ma.masked_greater(z,0)
		self.lat = lat
		self.lon = lon
		self.dz_dy = dz_dy
		self.dz_dx = dz_dx

	def get_index_from_pos(self,pos):
		isinstance(pos,geopy.Point)
		lon_index = self.lon.find_nearest(pos.longitude)
		lat_index = self.lat.find_nearest(pos.latitude)
		return (lon_index,lat_index)	

	def return_z(self,pos):
		lon_index,lat_index = self.get_index_from_pos(pos)
		return self.z[lat_index,lon_index]

	def return_gradient(self,pos):
		x_index,y_index = self.get_index_from_pos(pos)
		dz_dx = self.dz_dx/(np.cos(np.deg2rad(self.y[y_index])))
		return (dz_dx[y_index,x_index],self.dz_dy[y_index,x_index])

	def get_gradient(self):
		self.dz_dy,self.dz_dx = np.gradient(self.z.data,self.lat,self.lon)

	def griddata_subsample(self,array):
		array_mask = ~np.isnan(array)
		X,Y = np.meshgrid(self.x,self.y)
		lons = np.arange(-180,180,0.1)
		lats = np.arange(-90,90,0.1)
		XX,YY = np.meshgrid(lons,lats)
		z = griddata(np.array(list(zip(X[array_mask],Y[array_mask]))),array[array_mask],(XX,YY))
		return self.__class__(lon = lons,lat = lats,z=z)

	def stride_subsample(self,stride):
		lon = self.lon[::stride]
		lat = self.lat[::stride]
		z = self.z[::stride,::stride]
		return self.__class__(lon=lon,lat=lat,z=z)

	def guassian_smooth(self,sigma=5):
		z = gaussian_filter(self.z, sigma=sigma)
		return self.__class__(lon=self.lon,lat=self.lat,z=z)

	def regional_subsample(self,urlon,lllon,urlat,lllat):
		urlon_idx = self.lon.find_nearest(urlon,idx=True)
		lllon_idx = self.lon.find_nearest(lllon,idx=True)
		lllat_idx = self.lat.find_nearest(lllat,idx=True)
		urlat_idx = self.lat.find_nearest(urlat,idx=True)
		lon = self.lon[lllon_idx:(urlon_idx+1)]
		lat = self.lat[lllat_idx:(urlat_idx+1)]
		z = self.z[lllat_idx:(urlat_idx+1),lllon_idx:(urlon_idx+1)]
		return self.__class__(lon=lon,lat=lat,z=z)

class ETopo1Depth(DepthBase):
	def __init__(self,*args,**kwargs):
		super().__init__(*args, **kwargs)

	@staticmethod
	def load():
		nc_fid_z_data = Dataset(get_base_folder()+'Raw/Bathymetry/ETOPO1_Bed_c_gdal.grd')
		nc_fid_coord = Dataset(get_base_folder()+'Raw/Bathymetry/ETOPO1_Bed_g_gmt4.grd')
		lon = LonList(nc_fid_coord['x'][:-1])
		lat = LatList(nc_fid_coord['y'][:-1])
		z = nc_fid_z_data['z'][:].reshape(len(lat),len(lon))
		z = z[::-1,:]
		return ETopo1Depth(lon=lon,lat=lat,z=z)


class PACIOOS(DepthBase):
	def __init__(self,*args,**kwargs):
		super().__init__(*args, **kwargs)

	@staticmethod
	def load():
		nc_fid = Dataset(get_base_folder()+'Raw/Bathymetry/hmrg_bathytopo_1km_mhi.nc')
		lon = LonList(nc_fid['x'][:])
		lat = LatList(nc_fid['y'][:])
		z = nc_fid['z'][:]
		return PACIOOS(lon=lon,lat=lat,z=z)
