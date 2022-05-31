import pandas as pd
from data_save_utilities.lagrangian.drifter_base_class import BasePosition,Speed,BaseRead,all_drift_holder
from data_save_utilities.pickle_utilities import load,save
import datetime
import geopy
import numpy as np 


class MetaClass():
	""" class to organize all of read in meta net cdf data
		----------
		file_path: the file path of the net cdf file you wish to read

	"""
	def __init__(self,df_):
		self.id = df_.spotterId.tolist()[0]

class ProfClass():
	""" class to organize all of read in profile net cdf data
		----------
		file_path: the file path of the net cdf file you wish to read

	"""		
	def __init__(self,df_):
		self.date = self.Date(df_)
		self.pos = self.Position(df_)
		self.speed = Speed(self.date,self.pos,speed_limit=30)

	class Date(object):
		def __init__(self,df_):
			self._mask = [False]*len(df_)
			self._list = self.parse_time(df_)

		def parse_time(self,df_):
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
			fmt = '%Y-%m-%dT%H:%M:%S.000Z'
			return [datetime.datetime.strptime(dummy,fmt) for dummy in df_['timestamp'].tolist()]

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
			test_3 = ((np.array(self._list))>datetime.datetime(2016,1,1)).all() # make sure all profiles are after beginning of argo program
			test_4 = ((np.array([_.days for _ in time_diff_list])<270)).all() # make sure maximum time difference between profiles is less than 9 months
			bool_value = ~(test_1&test_2&test_3&test_4)
			if bool_value: 
				print 'date is the problem'
			return bool_value


	class Position(BasePosition):
		""" class that is a holder of the position information

		"""     
		def __init__(self,df_):			
			self._mask = [False]*len(df_)
			self._list = [geopy.Point(_[0],_[1]) for _ in zip(df_['latitude'].tolist(),df_['longitude'].tolist())]


class SpotRead(BaseRead):
	data_description = 'spotter'
	def __init__(self,df_):
		self.prof = ProfClass(df_)
		self.meta = MetaClass(df_)
		super(SpotRead,self).__init__()



def aggregate_spotter_list(num=-1):
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
	sofar_df = pd.read_csv('/Users/pchamberlain/Projects/sofar_oceanography_colab/Data/Spotter_track_data.csv',index_col=0)
	for spot in sofar_df.spotterId.unique()[:num]:
		print spot
		df_holder = sofar_df[sofar_df.spotterId==spot]
		dummy = SpotRead(df_holder[::48])
