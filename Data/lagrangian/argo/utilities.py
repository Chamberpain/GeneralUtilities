import datetime
import geopy
from GeneralUtilities.Data.Lagrangian.drifter_base_class import BasePosition,BaseSpeed,BaseTime
from GeneralUtilities.Data.Filepath.instance import get_data_folder
from netCDF4 import Dataset
import os
import re
import numpy as np
from GeneralUtilities.Data.Lagrangian.lagrangian_utilities import data_return,parse_time,julian_time_parse

def compile_classes(file_pattern='meta.nc'):
	#this should probably use the find files function in search utility
	data_file_name = get_data_folder()+'/Raw/Argo'
	matches = []
	for root, dirnames, filenames in os.walk(data_file_name):
		print(filenames)
		if any([file.endswith(file_pattern) for file in filenames]):
			matches.append(root)
	return matches

class MetaClass():
	""" class to organize all of read in meta net cdf data
		----------
		file_path: the file path of the net cdf file you wish to read

	"""
	name = 'Meta'
	def __init__(self,file_path):
		print('I am opening Meta File')
		nc_fid = Dataset(file_path)
		self.folder = os.path.dirname(file_path)
		self.id = data_return(nc_fid['PLATFORM_NUMBER'])
		try:
			self.platform_type = data_return(nc_fid['PLATFORM_TYPE'])
		except IndexError:
			self.platform_type = data_return(nc_fid['PLATFORM_MODEL'])
		self.project_name = data_return(nc_fid['PROJECT_NAME'])
		self.positioning_system = self._positioning_system_format(nc_fid.variables['POSITIONING_SYSTEM'])
		if bool(data_return(nc_fid.variables['LAUNCH_QC'])):
			launch_lat = data_return(nc_fid.variables['LAUNCH_LATITUDE'])
			launch_lon = data_return(nc_fid.variables['LAUNCH_LONGITUDE'])
			self.launch_loc = geopy.Point(launch_lat,launch_lon)
		self.launch_date = parse_time(nc_fid.variables['LAUNCH_DATE'])
		if bool(data_return(nc_fid.variables['START_DATE_QC'])):
			self.start_date = parse_time(nc_fid.variables['START_DATE'])

		nc_fid.close()


	def _positioning_system_format(self,data_instance):
		""" Performs necessary checks on data instance
			Parameters
			----------
			data_instance: net cdf positioning system data instance

			Returns
			-------
			positioning system string
		"""
		positioning_string = data_return(data_instance)
		if not positioning_string:
			return None
		if positioning_string in ['IRIDIUM','GTS','GPSIRIDIUM','IRIDIUMGPS','GPSIRIDIUMRAFOS','BEIDOU']:
			positioning_string = 'GPS'
		assert positioning_string in ['ARGOS','GPS']
		return positioning_string

def data_adjust(data_instance, data_adjusted_instance):
	"""
	Function replaces adjusted data in the original data array

	Parameters
	----------
	data_instance: net cdf data instance
	adjusted: adjusted net cdf data instance

	Returns
	-------
	data processed in a usable form (string, array, float)
	"""
	masked_data_array = data_instance[:].data
	masked_adjusted_instance = data_adjusted_instance[:]
	masked_data_array[~masked_adjusted_instance.mask]=masked_adjusted_instance[~masked_adjusted_instance.mask]
	masked_data_array = np.ma.masked_equal(masked_data_array,data_instance._FillValue)
	return masked_data_array



def format_byte_list_to_string(byte_list):
	return ['None' if x is None else x.decode("utf-8") for x in byte_list]

class BaseReadClass(object):
	pass
	""" base class that contains functions used by both ProfClass and TrajClass
		----------
	"""	
	class Position(BasePosition):
		""" class that is a holder of the position information

		"""
		max_distance_between = 926 #this is 500nm
		@classmethod
		def from_ncfid(cls,nc_fid,mask):
			pos_data = [geopy.Point(_[0],_[1]) for _ in zip(nc_fid['LATITUDE'][mask],nc_fid['LONGITUDE'][mask])]
			pos_list = cls(pos_data) 
			if not pos_list:
				print('the position information was empty, to prove it, take a look at')
				print(nc_fid['LATITUDE'][:])
				print(nc_fid['LONGITUDE'][:])
			return pos_list

class ArgoTime(BaseTime):
	project_start = datetime.datetime(1998,1,1)
	maximum_time_between_data = 270 # days

class ProfDate(ArgoTime):
	@classmethod
	def from_ncfid(cls,nc_fid,*args,**kwargs):
		date_data = nc_fid['JULD_LOCATION'][:]
		mask = (np.array([_ in ['1','2'] for _ in format_byte_list_to_string(nc_fid['POSITION_QC'][:].tolist())]))&\
		(~nc_fid['LONGITUDE'][:].mask)&\
		(~nc_fid['LATITUDE'][:].mask)&\
		(~date_data.mask)
		
		data_list = julian_time_parse(date_data[mask].data, parse_time(nc_fid['REFERENCE_DATE_TIME']))
		return (cls(data_list),mask)

class Speed(BaseSpeed):
	speed_limit = 2.57 # equivalent to 5kts 

def argo_file_mover(base_folder):
	#function made to make transferring argo snapshots to external harddrive easier
	import shutil
	transfer_to_base_folder_name = get_data_folder()+'/Raw/Argo/'
	matches_from = []
	matches_to = []
	for root, dirnames, filenames in os.walk(base_folder):
		for file in filenames:
			if any([file.endswith(x) for x in ['meta.nc','traj.nc','prof.nc','tech.nc']]):
				matches_from.append(os.path.join(root,file))
				matches_to.append(os.path.join(root.replace(base_folder,transfer_to_base_folder_name),file))
	for match_from,match_to in zip(matches_from,matches_to):
		print(match_from)
		print(matches_from.index(match_from))
		print('of')
		print(len(matches_from))
		try:
			shutil.copy2(match_from,match_to)
		except FileNotFoundError:
			os.makedirs(os.path.dirname(match_to))