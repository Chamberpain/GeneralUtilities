import os
from netCDF4 import Dataset
import datetime
import geopy
from data_save_utilities.pickle_utilities import load,save
import numpy as np
from compute_utilities.matrix_utilities import convert_lon_to_180
from data_save_utilities.lagrangian.drifter_base_class import BasePosition,Speed,BaseRead

class MetaClass():
	""" class to organize all of read in meta net cdf data
		----------
		file_path: the file path of the net cdf file you wish to read

	"""
	def __init__(self,nc_fid):
		self.id = nc_fid['ID'][:][0]

class ProfClass():
	""" class to organize all of read in profile net cdf data
		----------
		file_path: the file path of the net cdf file you wish to read

	"""		
	def __init__(self,nc_fid,skip):
		lon_var = 'lon360'
		if np.sum(nc_fid['longitude'][:].mask)<np.sum(nc_fid['lon360'][:].mask):
			lon_var = 'logitude'
		mask = (nc_fid['latitude'][::skip].mask)|(nc_fid[lon_var][::skip].mask)|(nc_fid['time'][::skip].mask)

		self.date = self.Date(nc_fid,mask,skip)
		self.pos = self.Position(nc_fid,lon_var,mask,skip)
		self.speed = Speed(self.date,self.pos,speed_limit=20)

	class Date(object):
		def __init__(self,nc_fid,mask,skip):
			self._mask = mask
			self._list = self.parse_time(nc_fid,skip)

		def parse_time(self,variable_instance,skip,reference_date=datetime.datetime(1970,1,1)):
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
			time_instance = [reference_date+datetime.timedelta(seconds=dummy) for dummy in variable_instance['time'][::skip][~self._mask].ravel()]
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
			test_2 = ((np.array(self._list))<datetime.datetime.today()).all() # make sure all profiles happen before today
			test_3 = ((np.array(self._list))>datetime.datetime(1970,1,1)).all() # make sure all profiles are after beginning of argo program
			test_4 = ((np.array([_.days for _ in time_diff_list])<270)).all() # make sure maximum time difference between profiles is less than 9 months
			bool_value = ~(test_1&test_2&test_3&test_4)
			if bool_value: 
				print 'date is the problem'
			return bool_value


	class Position(BasePosition):
		""" class that is a holder of the position information

		"""     
		def __init__(self,nc_fid,lon_var,mask,skip):			
			self._mask = mask
			lon_array = convert_lon_to_180(nc_fid[lon_var][::skip])
			self._list = [geopy.Point(_[0],_[1]) for _ in zip(nc_fid['latitude'][::skip][~self._mask].ravel(),lon_array[~self._mask].ravel())]


class AOMLRead(BaseRead):
	data_description = 'aoml'
	def __init__(self,file_path,skip=10):
		nc_fid = Dataset(file_path)
		self.prof = ProfClass(nc_fid,skip)
		self.meta = MetaClass(nc_fid)
		nc_fid.close()
		super(AOMLRead,self).__init__()


def aggregate_aoml_list(num=-1,list_save=False,recompile=False):
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
	__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))

	def compile_matches():
		base_folder = '/Users/pchamberlain/iCloud/Data/Raw/AOML/'
		matches = []
		for root, dirnames, filenames in os.walk(base_folder):
			for filename in filenames:
				if filename.endswith('.nc'):
					matches.append(os.path.join(root,filename))
				print len(matches)
		save("aoml_matches.pkl",matches)
		return matches

	# depth = Depth()
	if recompile:
		compile_matches()
	else:
		try:
			matches = load(__location__+'/aoml_matches.pkl')
		except IOError:
			matches = compile_matches()

	num = 0
	while num<len(matches):
		print 'I am opening number match ',num
		yield AOMLRead(matches[num])
		num +=1 

