import os
from GeneralUtilities.Data.lagrangian.argo.__init__ import ROOT_DIR
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

def data_return(data_instance):
	"""
	Function reads in netcdf data instance and returns sensible data

	Parameters
	----------
	data_instance: net cdf data instance
	mask: specifies whether to mask the values or return fill value
	Returns
	-------
	data processed in a usable form (string, array, float)
	"""
	data_type = data_instance.datatype
	masked_array = data_instance[:]
	if data_type in ['S1','S']:
		data = ''.join([_.decode("utf-8") for _ in masked_array[~masked_array.mask].data.ravel()])
	elif data_type == 'float64':
		data = [_ for _ in masked_array[~masked_array.mask].data.ravel()]
		if len(data)==1:
			data = data[0]
	else: 
		print(data_type)
		print('I dont know what to do with this data')
		raise
	return data

def format_byte_list_to_string(byte_list):
	return ['None' if x is None else x.decode("utf-8") for x in byte_list]

class BaseDate(object):
	name = 'Date'

	def julian_time_parse(self,array,reference_date):
		time_instance = [reference_date +datetime.timedelta(days=_) for _ in array]
		return time_instance

	def parse_time(self,variable_instance,reference_date=None):
		"""
		Function parses a netcdf time related instance and returns a time instance

		Parameters
		----------
		variable instance of the time you wish decyfered
		reference date for juld 

		Returns
		-------
		date time instance
		"""

		if variable_instance.conventions == 'YYYYMMDDHHMISS':
			time_string = data_return(variable_instance)
			if time_string: #this tests for an empty string
				try:
					time_instance = datetime.datetime.strptime(time_string,'%Y%m%d%H%M%S')
				except ValueError:
					time_instance = None
			else:
				time_instance=None
		elif variable_instance.conventions == 'Relative julian days with decimal part (as parts of day)':
			if not reference_date:
				print('I need a reference date for Julian data')
				raise
			time_instance = julian_time_parse(data_return(variable_instance),reference_date)
		else:
			print(variable_instance.conventions)
			print('I dont know what to do with this')
			raise
		return time_instance

	def return_mask(self):
		return self._mask

	def is_problem(self):
		"""
		Returns boolean value to show that all the date tests have been passed.

		Current tests:
		All time differences are greater than 0
		All profiles happen before today 
		All profiles happen after the beginning of the argo program
		The maximum time difference between profiles cannot be greater than 9 months

		Returns
		-------
		Boolean value to show if this profile is a problem
		"""

		time_diff_list = [self._list[idx+1]-self._list[idx] for idx in range(len(self._list)-1)]
		seconds_diff_list = [_.days*24+_.seconds/3600. for _ in time_diff_list]
		test_1 = (np.array(seconds_diff_list)>0).all()	# make sure all time is going in the right direction
		if not test_1:
			print('time is not going in the right direction')
		test_2 = ((np.array(self._list))<datetime.datetime.today()).all() # make sure all profiles happen before today
		if not test_2:
			print('a profile happened in the future')
		test_3 = ((np.array(self._list))>datetime.datetime(1998,1,1)).all() # make sure all profiles are after beginning of argo program
		if not test_2:
			print('a profile happened before argo existed')
		test_4 = ((np.array([_.days for _ in time_diff_list])<270)).all() # make sure maximum time difference between profiles is less than 9 months
		if not test_2:
			print('maximum time between profiles is greater than 270 days')
		return ~(test_1&test_2&test_3&test_4)


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
			if positioning_string in ['IRIDIUM','GTS','GPSIRIDIUM','IRIDIUMGPS','GPSIRIDIUMRAFOS']:
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
				test_1 = (np.array(self._list)>500).all()	# allow a window of 500 on either side of drift depth
				if not test_1:
					print('drift depth less than 500 mb')
				test_2 = (np.array(self._list)<1500).all() # 
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

	all_dict_filename = os.getenv("HOME")+'/Data/Raw/Argo/all_dict'
	try: 
		with open(all_dict_filename,'rb') as pickle_file:
			out_data = pickle.load(pickle_file)
			BaseRead.all_dict = out_data		
	except FileNotFoundError:
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
			read_class(matches[number])
			number += 1 
		with open(all_dict_filename, 'wb') as pickle_file:
			pickle.dump(BaseRead.all_dict,pickle_file)
		pickle_file.close()