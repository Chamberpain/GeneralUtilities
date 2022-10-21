import os
from GeneralUtilities.Data.Lagrangian.Argo.__init__ import ROOT_DIR
from netCDF4 import Dataset
import datetime
import geopy
import re
import numpy as np 
from GeneralUtilities.Data.Lagrangian.Argo.prof_class import BaseProfClass,prof_class_dict
from GeneralUtilities.Data.Lagrangian.Argo.traj_class import BaseTrajClass,traj_class_dict
from GeneralUtilities.Data.Lagrangian.Argo.tech_class import BaseTechClass,tech_class_dict
from GeneralUtilities.Data.Lagrangian.drifter_base_class import DrifterArrayBase,BaseRead
from GeneralUtilities.Data.Lagrangian.lagrangian_utilities import data_return,parse_time,julian_time_parse
from GeneralUtilities.Compute.list import DepthList
from GeneralUtilities.Data.pickle_utilities import load,save
from GeneralUtilities.Data.Filepath.instance import FilePathHandler
from GeneralUtilities.Data.Lagrangian.Argo.utilities import compile_classes
import pickle


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
			self.meta = self._file_reader(nc_folder,files,'.*meta.nc',MetaClass)
		if traj:
			try:
				TrajClass = traj_class_dict[self.meta.platform_type]
			except KeyError:
				TrajClass = BaseTrajClass
			self.traj = self._file_reader(nc_folder,files,'.*_traj.nc',TrajClass)
		if prof:
			try:
				ProfClass = prof_class_dict[self.meta.platform_type]
			except KeyError:
				ProfClass = BaseProfClass
			self.prof = self._file_reader(nc_folder,files,'.*prof.nc',ProfClass)
		if tech:
			try:
				TechClass = tech_class_dict[self.meta.platform_type]
			except KeyError:
				TechClass = BaseTechClass
			self.tech = self._file_reader(nc_folder,files,'.*tech.nc',TechClass)
		super(ArgoReader,self).__init__()

	def _file_reader(self,nc_folder,files,data_format,data_class):
		match_alg = re.compile(data_format)
		file_name = next(iter(list(filter(match_alg.match,files))+[None]))
		if file_name:
			print(os.path.join(nc_folder,file_name))
			return data_class(os.path.join(nc_folder,file_name))
		else:
			return None

	def get_variables(self):
		return ['PRES', 'TEMP', 'PSAL']

	def is_problem(self):
		def string_report(bool_list,dummy_class,dummy_class_subtype):
			if bool_list[-1]:
				print(dummy_class.name)
				print(dummy_class_subtype.name)
		bool_list = []
		class_list = []
		try:
			bool_list.append(self.tech.drift_depth.is_problem())
		except AttributeError:
			pass
		try:
			self.prof.pos._list
		except AttributeError:
			bool_list.append(True)
		try: 
			class_list.append(self.traj)
		except AttributeError:
			pass
		try:
			class_list.append(self.prof)
		except AttributeError:
			pass
		for dummy_class in class_list:
			if dummy_class:
				bool_list.append(dummy_class.pos.is_problem())
				string_report(bool_list,dummy_class,dummy_class.pos)
				bool_list.append(dummy_class.date.is_problem())
				string_report(bool_list,dummy_class,dummy_class.date)
				bool_list.append(dummy_class.speed.is_problem())
				string_report(bool_list,dummy_class,dummy_class.speed)
		truth_value = any(bool_list)
		if not bool_list:
			truth_value = True 
		if truth_value: 
			print('I am a problem, do not add me to anything')
		return truth_value

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
			self.platform = data_return(nc_fid['PLATFORM_TYPE'])
		except IndexError:
			self.platform = data_return(nc_fid['PLATFORM_MODEL'])
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
		compile_classes(num)
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

class ArgoArray(DrifterArrayBase):

	def __init__(self,*args,num=99999,**kwargs):
		super().__init__(*args, **kwargs)
		matches = compile_classes()
		number = 0
		for match in matches:
			if number>num:
				continue
			argo_instance = ArgoReader(match,prof=False,traj=False,tech=False)
			print('I am opening number match ',number)
			print ('filename is ',matches[number])
			# if not argo_instance.problem:
			self.update({(argo_instance.meta.id,argo_instance)})
			number += 1 


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
	def get_id_list():
		id_list = []
		for ii in ArgoReader.all_dict.values():
			id_list.append(ii.meta.id)
		return id_list

	@staticmethod
	def get_full_lat_lon_list():
		lat_list = []
		lon_list = []
		for ii in ArgoReader.all_dict.values():
			lat,lon = zip(*[(dummy.latitude,dummy.longitude) for dummy in ii.prof.pos._list])
			lat_list += [lat]
			lon_list += [lon]
		return (lat_list,lon_list)

	@staticmethod
	def recent_bins_by_sensor(variable,lat_bins,lon_bins):
		date_list = ArgoReader.get_recent_date_list()
		bin_list = ArgoReader.get_recent_bins(lat_bins,lon_bins)
		sensor_list = [['PRES', 'TEMP', 'PSAL'] for x in range(len(bin_list))]
		sensor_mask = [variable in x for x in sensor_list]
		date_mask =[max(date_list)-datetime.timedelta(days=180)<x for x in date_list]
		mask = np.array(sensor_mask)&np.array(date_mask)
		return np.array(bin_list)[mask]


