from OptimalArray.Utilities.CM4Mat import CovCM4Global
from GeneralUtilities.Data.Mapped.mapped_base import MappedBase
from GeneralUtilities.Compute.list import LonList, LatList
from GeneralUtilities.__init__ import ROOT_DIR
import matplotlib.pyplot as plt
import datetime
import os
import numpy as np
from netCDF4 import Dataset
from TransitionMatrix.Utilities.Utilities import shiftgrid
from OptimalArray.Utilities.CM4Mat import CovCM4

datadir = CovCM4.data_directory

class CM4(MappedBase):

	def return_dimensions(self):
		filename = dict([(x,y) for x,y in CovCM4Global.get_filenames()])[self.variable]
		nc_fid = Dataset(filename[0])
		lats = nc_fid["lat"][:]
		lons = nc_fid["lon"][:]
		data = nc_fid[self.variable][0, 0,:,:]
		data,lons = shiftgrid(180.5, data, lons, start=False)
		return (LonList(lons),LatList(lats))

	def return_dataset(self,depth_idx=2):
		master_list = CovCM4Global.get_filenames()
		file_list = dict([(x,y) for x,y in master_list])[self.variable]
		array_variable_list = []
		for file in np.sort(file_list):
			time_list = []
			holder_list = []
			nc_fid = Dataset(file)
			array_variable_list.append(nc_fid[self.variable][:, depth_idx,:,:])
		lons = nc_fid["lon"][:]
		data = np.vstack(array_variable_list)
		data,lons = shiftgrid(180.5, data, lons, start=False)
		data = np.ma.masked_greater(data, 10 ** 19)
		data.units = nc_fid[self.variable].units
		return (data,lons)

class CM4O2(CM4):
	variable = 'o2'

class CM4ThetaO(CM4):
	variable = 'thetao'

class CM4Sal(CM4):
	variable = 'so'

class CM4PH(CM4):
	variable = 'ph'

class CM4CHL(CM4):
	variable = 'chl'