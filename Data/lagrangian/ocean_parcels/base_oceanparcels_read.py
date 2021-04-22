from __future__ import print_function
import os
import datetime
import geopy
import numpy as np 
import geopy.distance
from data_save_utilities.lagrangian.drifter_base_class import Speed,BasePosition,BaseRead
from data_save_utilities.lagrangian.argo.argo_read import BaseDate
import xarray as xr

class Date(BaseDate):
	def __init__(self,t):
		self.name = 'Date'
		self._list = t
		self._mask = [False]*len(t)

class Position(BasePosition):
	def __init__(self,x,y):
		self.name = 'Position'
		self._list = [geopy.Point(_[1],_[0]) for _ in zip(x,y)]
		self._mask = [False]*len(x)

class OPReader(BaseRead):
	"""
	class holds all argo meta, trajectory, profile data from SOSE meta class objects

	Parameters
	----------
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
	def __init__(self,x,y,t,k,meta=True, traj= True, prof=True,tech=True):
		if meta:
			self.meta = self.MetaClass(x,y,t,k)
		if traj:
			self.traj = self.TrajClass(x,y,t,k)
		if prof:
			self.prof = self.ProfClass(x,y,t,k)
		if tech:
			self.tech = self.TechClass(x,y,t,k)
		super(OPReader,self).__init__()

	@classmethod
	def aggregate_OP_list(cls):
		max_idx = 5000
		dataset = xr.open_dataset(cls.file_path+'output.nc')
		lats = dataset.lat.data[:,:max_idx]
		lats = lats - lats.min()
		lats = lats/lats.max()
		lats = lats*30
		lons = dataset.lon.data[:,:max_idx]
		lons = lons - lons.min()
		lons = lons/lons.max()
		lons = lons*360-180

		ref_date = datetime.datetime(year=2000,month=1,day=1)
		
		time = dataset.time.data[0,:max_idx]/dataset.time.data[0,:max_idx].max()
		date = [ref_date + datetime.timedelta(days=_) for _ in time*20*365]

		k=0
		while k < dataset.lat.shape[0]:
			print(k)
			x = lons[k,:][:-2][::2]
			y = lats[k,:][:-2][::2]
			yield OPReader(x,y,date[:-2][::2],k)
			k +=1



	class MetaClass():
		""" class to organize all of read in meta data
			----------

		"""
		def __init__(self,x,y,t,k):
			self.id = k
			self.date = Date(t)
			self.positioning_system = 'model_output'


	class TrajClass():
		""" class to organize all of read in trajectory data
			----------

		"""	
		def __init__(self,x,y,t,k):
			self.name = 'Traj'
			self.date = Date(t)
			self.pos = Position(x,y)
			self.speed = Speed(self.date,self.pos,speed_limit=5)




		# def assign_depth(self,depth):
		# 	try:
		# 		self.depth = [depth.return_z(_) for _ in self.pos]

	class ProfClass(BaseRead):
		""" class to organize all of read in profile data
			----------

		"""		
		def __init__(self,x,y,t,k):
			self.name = 'Prof'
			self.date = Date(t)
			self.pos = Position(x,y)
			self.speed = Speed(self.date,self.pos,speed_limit=5)

	class TechClass():
		def __init__(self,x,y,t,k):
			self.name = 'Tech'
			self.drift_depth = 1000


