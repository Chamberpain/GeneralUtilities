import os
from GeneralUtilities.Data.lagrangian.argo.__init__ import ROOT_DIR
from netCDF4 import Dataset
import datetime
import geopy
import re
import numpy as np 
import geopy.distance
from GeneralUtilities.Data.lagrangian.drifter_base_class import BasePosition,Speed,BaseRead,BaseDate,data_return
from GeneralUtilities.Data.pickle_utilities import load,save
from GeneralUtilities.Filepath.instance import FilePathHandler
import pickle



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



class ArgoReader(BaseRead):
	"""
	class holds all argo meta, trajectory, profile data from argo meta class objects

	Parameters
	----------
	folder containing all data file opened with the netCDF4 Dataset library
	meta: boolean object that instructs file to open meta data netcdf file
	traj: boolean object that instructs file to open traj netcdf file and read
	prof: boolean object that instructs file to open profile netcdf file and read
	** note **
	these boolean objects are only to improve efficiency. By default all data will be loaded

	Returns
	-------
	Class object of netcdf file data.
	All lat/lon information formatted to -180 to 180
	All times are returned as datetime objects
	"""
	data_description = 'argo'
	def __init__(self,nc_folder,meta=True, traj= True, prof=True,tech=True):
		self.folder = nc_folder
		files = os.listdir(nc_folder)
		if meta:
			self.meta = self._file_reader(nc_folder,files,'.*meta.nc',self.MetaClass)
		if traj:
			self.traj = self._file_reader(nc_folder,files,'.*_traj.nc',self.TrajClass)
		if prof:
			self.prof = self._file_reader(nc_folder,files,'.*prof.nc',self.ProfClass)
		if tech:
			self.tech = self._file_reader(nc_folder,files,'.*tech.nc',self.TechClass)
		super(ArgoReader,self).__init__()

	def _file_reader(self,nc_folder,files,data_format,data_class):
		match_alg = re.compile(data_format)
		file_name = next(iter(list(filter(match_alg.match,files))+[None]))
		if file_name:
			print(os.path.join(nc_folder,file_name))
			return data_class(os.path.join(nc_folder,file_name))
		else:
			return None

	@staticmethod
	def get_pos_breakdown():
		argos_list = []
		gps_list = []
		for _ in ArgoReader.all_dict.items():
			if _[1].meta.positioning_system =='ARGOS':
				argos_list.append(_[0])
			if _[1].meta.positioning_system =='GPS':
				gps_list.append(_[0])
		return (argos_list,gps_list)

	@staticmethod
	def get_pos_list():
		pos_list = []
		for ii in ArgoReader.all_dict.values():
			pos_list.append(ii.meta.positioning_system)
		return pos_list

	@staticmethod
	def get_full_lat_lon_list():
		lat_list = []
		lon_list = []
		for ii in ArgoReader.all_dict.values():
			lat,lon = zip(*[(dummy.latitude,dummy.longitude) for dummy in ii.prof.pos._list])
			lat_list += [lat]
			lon_list += [lon]
		return (lat_list,lon_list)

	def get_variables(self):
		return ['PRES', 'TEMP', 'PSAL']

	@staticmethod
	def recent_bins_by_sensor(variable,lat_bins,lon_bins):
		date_list = ArgoReader.get_recent_date_list()
		bin_list = ArgoReader.get_recent_bins(lat_bins,lon_bins)
		sensor_list = [['PRES', 'TEMP', 'PSAL'] for x in range(len(bin_list))]
		sensor_mask = [variable in x for x in sensor_list]
		date_mask =[max(date_list)-datetime.timedelta(days=180)<x for x in date_list]
		mask = np.array(sensor_mask)&np.array(date_mask)
		return np.array(bin_list)[mask]

	@staticmethod
	def compile_classes(num):
		#this should probably use the find files function in search utility
		data_file_name = os.getenv("HOME")+'/Data/Raw/Argo'
		def compile_matches():
			matches = []
			for root, dirnames, filenames in os.walk(data_file_name):
				meta_match = re.compile('.*meta.nc') # all folders will have a meta file
				if any([file.endswith('meta.nc') for file in filenames]):
					matches.append(root)
			return matches

		matches = compile_matches()
		number = 0
		while number<len(matches[:num]):
			print('I am opening number match ',number)
			print ('filename is ',matches[number])
			ArgoReader(matches[number])
			number += 1 

	class BaseReadClass(object):
		pass
		""" base class that contains functions used by both ProfClass and TrajClass
			----------

		"""	
		class Position(BasePosition):
			""" class that is a holder of the position information

			"""     
			def __init__(self,nc_fid,mask):			
				self._list = [geopy.Point(_[0],_[1]) for _ in zip(nc_fid['LATITUDE'][mask],nc_fid['LONGITUDE'][mask])]
				if not self._list:
					print('the position information was empty, to prove it, take a look at')
					print(nc_fid['LATITUDE'][:])
					print(nc_fid['LONGITUDE'][:])

	class MetaClass():
		""" class to organize all of read in meta net cdf data
			----------
			file_path: the file path of the net cdf file you wish to read

		"""
		name = 'Meta'
		def __init__(self,file_path):
			print('I am opening Meta File')
			nc_fid = Dataset(file_path)
			self.id = data_return(nc_fid['PLATFORM_NUMBER'])
			self.project_name = data_return(nc_fid['PROJECT_NAME'])
			self.date = self.Date(nc_fid)
			self.positioning_system = self._positioning_system_format(nc_fid.variables['POSITIONING_SYSTEM'])
			if bool(data_return(nc_fid.variables['LAUNCH_QC'])):
				launch_lat = data_return(nc_fid.variables['LAUNCH_LATITUDE'])
				launch_lon = data_return(nc_fid.variables['LAUNCH_LONGITUDE'])
				self.launch_loc = geopy.Point(launch_lat,launch_lon)
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

		class Date(BaseDate):
			def __init__(self,nc_fid):
				self.launch_date = self.parse_time(nc_fid.variables['LAUNCH_DATE'])
				if bool(data_return(nc_fid.variables['START_DATE_QC'])):
					self.start_date = self.parse_time(nc_fid.variables['START_DATE'])

	class TrajClass(BaseReadClass):
		""" class to organize all of read in trajectory net cdf data
			----------
			file_path: the file path of the net cdf file you wish to read

		"""	
		name = 'Traj'
		def __init__(self,file_path):
			print('I am opening Traj File')
			nc_fid = Dataset(file_path)
			self.date = self.Date(nc_fid)
			mask = self.date.return_mask()
			self.position_accuracy = nc_fid['POSITION_ACCURACY'][mask]
			self.pos = self.Position(nc_fid,mask)
			self.speed = Speed(self.date,self.pos,speed_limit=5)
			nc_fid.close()

		class Date(BaseDate):
			def __init__(self,nc_fid):
				try:
					date_data = data_adjust(nc_fid['JULD'],nc_fid['JULD_ADJUSTED'])
				except IndexError:
					date_data = nc_fid['JULD'][:]
				mask = (np.array([_ in ['1','2'] for _ in format_byte_list_to_string(nc_fid['POSITION_QC'][:].tolist())]))&\
				(~nc_fid['LONGITUDE'][:].mask)&\
				(~nc_fid['LATITUDE'][:].mask)&\
				(~nc_fid['POSITION_ACCURACY'][:].mask)&\
				(~date_data.mask)
				self._list = self.julian_time_parse(date_data[mask].data, self.parse_time(nc_fid['REFERENCE_DATE_TIME']))
				self._mask = mask


		# def assign_depth(self,depth):
		# 	try:
		# 		self.depth = [depth.return_z(_) for _ in self.pos]

	class ProfClass(BaseReadClass):
		""" class to organize all of read in profile net cdf data
			----------
			file_path: the file path of the net cdf file you wish to read

		"""		
		name = 'Profile'
		def __init__(self,file_path):
			print('I am opening Prof File')
			nc_fid = Dataset(file_path)
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

	class TechClass():
		name = 'Tech'
		def __init__(self,file_path):
			print('I am opening Tech File')
			nc_fid = Dataset(file_path)
			self.drift_depth = self.DriftDepth(nc_fid)


		class DriftDepth():
			def __init__(self,nc_fid):


				variable_name_list = [''.join(format_byte_list_to_string(i)).replace(" ","") for i in nc_fid['TECHNICAL_PARAMETER_NAME'][:].data]
				
				indices = [] 
				variable_list = ['PRESSURE_InternalVacuumProfileStart_mbar','PRES_ParkMaximum_dBAR','PRES_ParkMinimum_dBAR']
				indices += [i for i, x in enumerate(variable_name_list) if x in variable_list]
				# if not indices:
				# 	print(np.unique(variable_name_list))

				holder = [''.join(format_byte_list_to_string(nc_fid['TECHNICAL_PARAMETER_VALUE'][:].data[_,:])).replace(" ","") for _ in indices]
				self._list = []
				for dummy in holder:
					self._list.append(int(''.join([_ for _ in dummy if _.isdigit()])))
			


			def is_problem(self):
				print(self._list)
				if not self._list:
					return False
				test_1 = sum(np.array(self._list)<500)<3	# allow a window of 500 on either side of drift depth
				if not test_1:
					print('drift depth less than 500 mb')
				test_2 = sum(np.array(self._list)>1500)<3 # 
				if not test_2:
					print('drift depth greater than 1500 mb')
				return ~(test_1&test_2)

def aggregate_argo_list(read_class=ArgoReader,num=-1):
	"""
	function that returns dictionary of argo read classes

	Parameters
	----------
	num: the length of the dictionary you want (default is all files)
	list_save: whether you want to save the list (default is false)

	Returns
	-------
	dictionary of argo read classes.
	"""

	all_dict_filename = os.getenv("HOME")+'/Data/Raw/Argo/all_dict_'+read_class.data_description
	try: 
		with open(all_dict_filename,'rb') as pickle_file:
			out_data = pickle.load(pickle_file)
			BaseRead.all_dict = out_data		
	except FileNotFoundError:
		read_class.compile_classes(num)
		with open(all_dict_filename, 'wb') as pickle_file:
			pickle.dump(BaseRead.all_dict,pickle_file)
		pickle_file.close()

def full_argo_list():
	import copy
	from GeneralUtilities.Data.lagrangian.bgc.bgc_read import BGCReader
	aggregate_argo_list()
	all_dict_holder = copy.deepcopy(BaseRead.all_dict)
	aggregate_argo_list(read_class=BGCReader)
	for float_num,float_class in BGCReader.all_dict.items():
		all_dict_holder[float_num] = float_class
	BaseRead.all_dict = all_dict_holder
