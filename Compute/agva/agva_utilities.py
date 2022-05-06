from GeneralUtilities.Data.depth.depth_utilities import DepthBase
from netCDF4 import Dataset
import os
import numpy as np
from GeneralUtilities.Compute.constants import degree_dist
from GeneralUtilities.Filepath.instance import get_base_folder
from GeneralUtilities.Compute.list import LatList,LonList

class AVGAStream(DepthBase):
	def __init__(self,*args,**kwargs):
		super().__init__(*args, **kwargs)

	@classmethod
	def load(cls,depth_level=18):
		#depth level 18 corresponds to 1000 meters
		base_folder = get_base_folder()
		file_path = base_folder+'/Processed/AGVA/'
		data_list = []
		for file in os.listdir(file_path):
			full_path = os.path.join(file_path,file)
			nc_fid = Dataset(full_path)
			data_list.append(nc_fid['geostrophic_streamfunction'][:,depth_level,:,:])
		data = np.vstack(data_list)
		z = np.nanmean(data,axis=0)
		z = z-np.nanmax(z)
		y = LatList(nc_fid['latitude'][:].tolist())
		x = LonList(np.arange(-180,180).tolist())

		holder = np.zeros(z.shape)
		holder[:,:180] = z[:,180:]
		holder[:,180:] = z[:,:180]
		z = np.ma.masked_array(holder)
		z[np.isnan(z)] = np.nanmin(z)*10
		return cls(lon=x,lat=y,z=z)