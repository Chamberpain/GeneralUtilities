import os
from GeneralUtilities.Data.lagrangian.bgc.__init__ import ROOT_DIR
from GeneralUtilities.Data.lagrangian.argo.argo_read import ArgoReader,BaseDate,data_adjust,data_return,format_byte_list_to_string
from netCDF4 import Dataset
import datetime
import geopy
import re
import numpy as np 
import geopy.distance
from GeneralUtilities.Data.lagrangian.drifter_base_class import BasePosition,Speed,BaseRead
from GeneralUtilities.Data.pickle_utilities import load,save
from GeneralUtilities.Filepath.instance import FilePathHandler
import pickle

class BGCReader(ArgoReader):
	data_description = 'bgc_argo'
	def __init__(self,bgc_folder,*args,**kwargs):
		self.bgc_folder = bgc_folder
		files = os.listdir(bgc_folder)
		self.bgc_prof = self._file_reader(bgc_folder,files,'.*prof.nc',self.BGCProfClass)
		super().__init__(*args,**kwargs)

	def get_variables(self):
		return self.bgc_prof.variables

	@staticmethod
	def recent_bins_by_sensor(variable,lat_bins,lon_bins):
		date_list = BGCReader.get_recent_date_list()
		bin_list = BGCReader.get_recent_bins(lat_bins,lon_bins)
		sensor_list = BGCReader.get_sensors()
		sensor_mask = [variable in x for x in sensor_list]
		date_mask =[max(date_list)-datetime.timedelta(days=180)<x for x in date_list]
		mask = np.array(sensor_mask)&np.array(date_mask)
		return np.array(bin_list)[mask]

	@staticmethod
	def compile_classes(num):
		data_file_name = os.getenv("HOME")+'/Data/Raw/Argo'
		bgc_matches = []
		argo_matches = []
		for root, dirnames, filenames in os.walk(data_file_name):
			if any([file.endswith('Sprof.nc') for file in filenames]):
				bgc_matches.append(root)
				file_list = root.split('/')
				file_list[-2] = file_list[-2].split('_')[0]
				argo_matches.append('/'.join(file_list))
				BGCReader(bgc_matches[-1],argo_matches[-1])


	class BGCProfClass(ArgoReader.BaseReadClass):
		""" class to organize all of read in profile net cdf data
			----------
			file_path: the file path of the net cdf file you wish to read

		"""		
		name = 'BGCProfile'
		def __init__(self,file_path):
			print('I am opening BGCProf File')
			nc_fid = Dataset(file_path)
			self.variables = [data_return(x) for x in nc_fid['STATION_PARAMETERS'][:][-1]]
			self.date = self.Date(nc_fid)
			mask = self.date.return_mask()
			self.pos = self.Position(nc_fid,mask)
			self.speed = Speed(self.date,self.pos,speed_limit=5)
			nc_fid.close()

		class Date(BaseDate):
			def __init__(self,nc_fid):
				date_data = nc_fid['JULD_LOCATION'][:]
				mask = (np.array([_ in ['1','2'] for _ in format_byte_list_to_string(nc_fid['POSITION_QC'][:].tolist())]))&\
				(~nc_fid['LONGITUDE'][:].mask)&\
				(~nc_fid['LATITUDE'][:].mask)&\
				(~date_data.mask)
				
				self._list = self.julian_time_parse(date_data[mask].data, self.parse_time(nc_fid['REFERENCE_DATE_TIME']))
				self._mask = mask