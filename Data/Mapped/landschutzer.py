from GeneralUtilities.Data.mapped.mapped_base import MappedBase
from GeneralUtilities.Compute.list import LonList, LatList
from GeneralUtilities.__init__ import ROOT_DIR
from netCDF4 import Dataset
import os
import matplotlib.pyplot as plt
from TransitionMatrix.Utilities.Utilities import shiftgrid
import datetime
import numpy as np
from GeneralUtilities.Data.Download.landschutzer_download import data_folder as datadir

class Landschutzer(MappedBase):

	def return_dimensions(self):
		filename = datadir+'/spco2_1982-2015_MPI_SOM-FFN_v2016.nc'
		nc_fid = Dataset(filename)
		lats = nc_fid["lat"][:]
		lons = nc_fid["lon"][:]
		time = [datetime.datetime(2000,1,1)+datetime.timedelta(seconds=int(x)) for x in nc_fid["time"][:]]
		return (LonList(lons),LatList(lats))

	def return_dataset(self,depth_idx=0):
		filename = datadir+'/spco2_1982-2015_MPI_SOM-FFN_v2016.nc'
		nc_fid = Dataset(filename)
		lons = nc_fid['lon'][:]
		data = nc_fid["fgco2_smoothed"][:]
		data = np.ma.masked_greater(data, 10 ** 19)
		data.units = nc_fid["fgco2_smoothed"].units
		return (data,lons)
