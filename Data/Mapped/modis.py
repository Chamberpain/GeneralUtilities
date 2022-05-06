from GeneralUtilities.Data.mapped.mapped_base import MappedBase
from GeneralUtilities.Compute.list import LonList, LatList
from GeneralUtilities.__init__ import ROOT_DIR
from netCDF4 import Dataset
import numpy as np
import os
import matplotlib.pyplot as plt
import datetime
from TransitionMatrix.Utilities.Utilities import shiftgrid
datadir = os.path.join(ROOT_DIR,'DataDir/Processed/MODIS/')[1:]

class MODIS(MappedBase):

	def return_dimensions(self):
		filename = datadir+'/GMIS_A_CHLA_01_2003.nc'
		nc_fid = Dataset(filename)
		lons = nc_fid['lon'][:][::4]
		lats = nc_fid['lat'][:][::4]
		return (LonList(lons),LatList(lats))

	def return_dataset(self,depth_idx=0):
		lons,lats = self.return_dimensions()
		holder_list = []
		for year in range(2003,2018):
			for month in range(1,13):
				time = datetime.date(year,month,1)
				filename = 'GMIS_A_CHLA_%02d_%4d.nc' % (month,year)
				holder = Dataset(datadir+filename)
				holder_list.append(holder['Chl_a'][:][::4,::4])
		data = np.stack(holder_list)
		data = np.ma.masked_equal(data, -9999)
		data.units = holder["Chl_a"].units
		return (data,lons)
