from __future__ import print_function
import os
import datetime
import geopy
import numpy as np 
import geopy.distance
from GeneralUtilities.Data.lagrangian.drifter_base_class import Speed,BasePosition,BaseRead
from GeneralUtilities.Data.lagrangian.argo.argo_read import BaseDate
import pickle


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

class SOSEReader(BaseRead):
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
	data_description = 'SOSE'
	def __init__(self,x,y,t,k,meta=True, traj= True, prof=True,tech=True):
		if meta:
			self.meta = self.MetaClass(x,y,t,k)
		if traj:
			self.traj = self.TrajClass(x,y,t,k)
		if prof:
			self.prof = self.ProfClass(x,y,t,k)
		if tech:
			self.tech = self.TechClass(x,y,t,k)
		super(SOSEReader,self).__init__()



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



def aggregate_sose_list():
	"""
	function that returns dictionary of SOSE read classes

	Parameters
	----------


	Returns
	-------
	dictionary of argo read classes.
	"""

	all_dict_filename = os.getenv("HOME")+'/Data/Raw/SOSE/all_dict'
	
	try:
		with open(all_dict_filename,'rb') as pickle_file:
			out_data = pickle.load(pickle_file)
			BaseRead.all_dict = out_data
	except:
		data_file_name = os.getenv("HOME")+'/Data/Raw/SOSE/particle_release/SO_RAND_0001.XYZ.0000000001.0003153601.data'
		npts= 10000
		ref_date = datetime.datetime(year=2010,month=1,day=1)
		opt=np.fromfile(data_file_name,'>f4').reshape(-1,3,npts)
		print("data has %i records" %(opt.shape[0]))

		k = 0
		while k < opt.shape[2]:
			print(k)
			x,y=opt[:,0,k],opt[:,1,k]#this is in model grid index coordinate, convert to lat-lon using x=x/6.0;y=y/6.0-77.875
			t = range(opt.shape[0])
			date=[ref_date+datetime.timedelta(days=10*j) for j in t]
			x=x/6.0;y=y/6.0-77.875
			x = x%360
			x[x>180]=x[x>180]-360
			SOSEReader(x,y,date,k)
			k +=1 

		with open(all_dict_filename, 'wb') as pickle_file:
			pickle.dump(BaseRead.all_dict,pickle_file)
		pickle_file.close()

