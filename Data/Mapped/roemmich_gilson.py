from GeneralUtilities.Data.mapped.mapped_base import MappedBase
from GeneralUtilities.Compute.list import LonList, LatList
from GeneralUtilities.__init__ import ROOT_DIR
from netCDF4 import Dataset
import os
import matplotlib.pyplot as plt
from TransitionMatrix.Utilities.Utilities import shiftgrid
from GeneralUtilities.Data.Download.roemmich_gilson_download import data_folder as datadir

class RoemmichGilsonSal(MappedBase):

	def return_dimensions(self):
		filename = datadir+'/RG_ArgoClim_Salinity_2019.nc'
		nc_fid = Dataset(filename)
		data,lons = self.return_dataset()
		lats = nc_fid['LATITUDE'][:]
		pres = nc_fid['PRESSURE'][:]
		return (LonList(lons),LatList(lats))

	def return_dataset(self,depth_idx=0):
		filename = datadir+'/RG_ArgoClim_Salinity_2019.nc'
		nc_fid = Dataset(filename)
		lons = nc_fid['LONGITUDE'][:]
		data = nc_fid["ARGO_SALINITY_ANOMALY"][:, depth_idx, :, :] + nc_fid["ARGO_SALINITY_MEAN"][depth_idx, :, :] 
		data,lons = shiftgrid(180.5, data, lons, start=False)
		data.units = nc_fid["ARGO_SALINITY_ANOMALY"].units
		return (data,lons)

class RoemmichGilsonTemp(MappedBase):

	def return_dimensions(self):
		filename = datadir+'/RG_ArgoClim_Temperature_2019.nc'
		nc_fid = Dataset(filename)
		data,lons = self.return_dataset()
		lats = nc_fid['LATITUDE']
		pres = nc_fid['PRESSURE']
		return (LonList(lons),LatList(lats))

	def return_dataset(self,depth_idx):
		filename = datadir+'/RG_ArgoClim_Temperature_2019.nc'
		nc_fid = Dataset(filename)
		lons = nc_fid['LONGITUDE']
		data = nc_fid["ARGO_TEMPERATURE_ANOMALY"][:, depth_idx, :, :] + nc_fid["ARGO_TEMPERATURE_MEAN"][depth_idx, :, :]
		data,lons = shiftgrid(180.5, data, lons, start=False)
		data.units = nc_fid["ARGO_TEMPERATURE_ANOMALY"].units
		return (data,lons)