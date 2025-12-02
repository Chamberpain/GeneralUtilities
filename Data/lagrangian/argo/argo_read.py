import os
from GeneralUtilities.Data.Lagrangian.Argo.__init__ import ROOT_DIR
from netCDF4 import Dataset
import datetime
import geopy
import re
import numpy as np 
from GeneralUtilities.Data.Lagrangian.Argo.prof_class import BaseProfClass,BGCProfClass,prof_class_dict
from GeneralUtilities.Data.Lagrangian.Argo.traj_class import BaseTrajClass,traj_class_dict
from GeneralUtilities.Data.Lagrangian.Argo.tech_class import BaseTechClass,tech_class_dict
from GeneralUtilities.Data.Lagrangian.drifter_base_class import BaseRead
from GeneralUtilities.Compute.list import DepthList
from GeneralUtilities.Data.pickle_utilities import load,save
from GeneralUtilities.Data.Filepath.instance import FilePathHandler
from GeneralUtilities.Data.Lagrangian.Argo.utilities import MetaClass
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
		files = [x for x in files if not x.startswith('.')]
		try:
			if meta:
				self.meta = self._file_reader(nc_folder,files,'meta.nc',MetaClass)
			if traj:
				try:
					TrajClass = traj_class_dict[self.meta.platform_type]
				except KeyError:
					TrajClass = BaseTrajClass
				self.traj = self._file_reader(nc_folder,files,'traj.nc',TrajClass)
			if prof:
				try:
					sprof_folder = nc_folder.replace('_core','_bgc')
					sprof_files = os.listdir(sprof_folder)
					sprof_files = [x for x in sprof_files if x.endswith('Sprof.nc')]
					self.prof = self._file_reader(sprof_folder,sprof_files,'prof.nc',BGCProfClass)
					if not self.prof:
						raise ValueError
					print('******  I have opened a BGC file ******')
				except:
					try:
						ProfClass = prof_class_dict[self.meta.platform_type]
					except KeyError:
						ProfClass = BaseProfClass
					self.prof = self._file_reader(nc_folder,files,'prof.nc',ProfClass)
			if tech:
				try:
					TechClass = tech_class_dict[self.meta.platform_type]
				except KeyError:
					TechClass = BaseTechClass
				self.tech = self._file_reader(nc_folder,files,'tech.nc',TechClass)
		except AssertionError:
			print('I have spawned an assertion error during creation')
		super(ArgoReader,self).__init__()

	def _file_reader(self,nc_folder,files,data_format,data_class):
		files = [x for x in files if x.endswith(data_format)]
		if files:
			print(os.path.join(nc_folder,files[0]))
			return data_class(os.path.join(nc_folder,files[0]))
		else:
			return None

	def get_variables(self):
		return self.prof.variables

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
			class_list.append(self.prof)
		except AttributeError:
			pass
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
				try:
					bool_list.append(dummy_class.pos.is_problem())
					string_report(bool_list,dummy_class,dummy_class.pos)
					bool_list.append(dummy_class.date.is_problem())
					string_report(bool_list,dummy_class,dummy_class.date)
					bool_list.append(dummy_class.speed.is_problem())
					string_report(bool_list,dummy_class,dummy_class.speed)
				except:
					print('I have spawned an error')
					print('Definitely do not add me to anything')
					return True
		truth_value = any(bool_list)
		if not bool_list:
			truth_value = True 
		if truth_value: 
			print('I am a problem, do not add me to anything')
		return truth_value