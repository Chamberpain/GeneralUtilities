from __future__ import print_function
from data_save_utilities.depth.depth_utilities import DepthBase
from netCDF4 import Dataset
import os
import numpy as np
from compute_utilities.constants import degree_dist
from data_save_utilities.file_path_utilities import get_base_folder

class AVGAStream(DepthBase):
	def __init__(self,depth_level=18):
		#depth level 18 corresponds to 1000 meters
		base_folder = get_base_folder()
		file_path = base_folder+'Processed/AGVA/'
		data_list = []
		for file in os.listdir(file_path):
			full_path = os.path.join(file_path,file)
			nc_fid = Dataset(full_path)
			data_list.append(nc_fid['geostrophic_streamfunction'][:,depth_level,:,:])
		data = np.vstack(data_list)
		self.z = np.nanmean(data,axis=0)
		self.y = nc_fid['latitude'][:]
		self.x = nc_fid['longitude'][:]
		self.x[self.x>180] = self.x[self.x>180]-360
		self.x,self.y,self.z = self.calculate_subsample(self.z)
		self.z[np.isnan(self.z)] = np.nanmin(self.z)*10

		dz_dy,dz_dx = np.gradient(self.z,self.y,self.x)
		self.dz_dx = dz_dx/degree_dist
		self.dz_dy = dz_dy/degree_dist