from __future__ import print_function
import geopy
import geopy.distance
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import random
import datetime
from GeneralUtilities.Compute.list import GeoList

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
	try:
		data_type = data_instance.datatype
	except AttributeError:
		data_type = data_instance.dtype
	masked_array = data_instance[:]
	if data_type in ['S1','S']:
		data = ''.join([_.decode("utf-8",errors='ignore') for _ in masked_array[~masked_array.mask].data.ravel()])
	elif data_type == 'float64':
		data = [_ for _ in masked_array[~masked_array.mask].data.ravel()]
		if len(data)==1:
			data = data[0]
	else: 
		print(data_type)
		print('I dont know what to do with this data')
		raise
	return data

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



class BasePosition(object):
	name = 'Position'

	def plot(self):
		lat,lon = zip(*[(_.latitude,_.longitude) for _ in self._list])
		plt.plot(lon,lat)

	def return_pos_bins(self,lat_bins,lon_bins,index_values=False):
		""" Returns bins of position list in index form or actual 

			Parameters
			----------
			lat_bins: list of latitudes of the grid
			lon_bins: list of longitudes of the grid
			index_values: only set you want to include

			Returns
			-------
			list of grid indices or actual lat/lon 
		"""       
		if index_values:
			pos_list = list(filter(None,np.array(self._list)[index_values].tolist()))
		else:
			pos_list = list(filter(None,self._list))
		lat_bin_index = lat_bins.digitize([x.latitude for x in pos_list])
		lon_bin_index = lon_bins.digitize([x.longitude for x in pos_list])
		lats = [lat_bins[x-1] for x in lat_bin_index]
		lons = [lon_bins[x-1] for x in lon_bin_index]
		return [geopy.Point(x) for x in zip(lats,lons)]

	def is_problem(self):
		"""
		Returns boolean value to show that all the date tests have been passed.

		Current tests:
		position difference is greater than 500 nm
		position difference is zero
		longitude cant be greater than 180 or less than -180
		latitude cant be greater than 90 or less than -90

		Returns
		-------
		Boolean value to show if this profile is a problem
		"""

		difference_list = [geopy.distance.great_circle(self._list[idx],self._list[idx+1]).nm for idx in range(len(self._list)-1)]
		if not self._list:
			print('position list is empty')
			return True
		lat_list,lon_list = zip(*[(_.latitude,_.longitude) for _ in self._list])
		test_1 = (np.array(difference_list)<500).all()
		if not test_1:
			print('difference between positions is greater than 500 nm')				
		test_2 = sum(np.array(difference_list)==0)<3
		if not test_2:
			print('the difference in position is 0')
			print(difference_list)
		test_3 = (np.array(lon_list)<=180).all()
		if not test_3:
			print('longitude is greater than 180')
		test_4 = (np.array(lon_list)>=-180).all()
		if not test_4:
			print('longitude is less than -180')
		test_5 = (np.array(lat_list)<=90).all()
		if not test_5:
			print('latitude is greater than 90')				
		test_6 = (np.array(lat_list)>=-90).all()
		if not test_6:
			print('latitude is less than -90')

		return ~(test_1&test_2&test_3&test_4&test_5&test_6)

class BaseRead(object):
	all_dict = {}
	def __init__(self):
		if not self.is_problem():
			BaseRead.all_dict.update({(self.meta.id,self)})
	
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

	@staticmethod
	def get_full_speed_list():
		full_speed_list = []

		for k,_ in enumerate(BaseRead.all_dict.values()):
			speed = _.prof.speed._list
			full_speed_list+=speed.tolist()
		return (full_speed_list)

	@staticmethod
	def get_full_date_list():
		full_date_list = []
		for _ in BaseRead.all_dict.values():
			full_date_list+=_.prof.date._list
		return full_date_list

	@staticmethod
	def get_recent_date_list():
		recent_date_list = []
		for x in BaseRead.all_dict.values():
			recent_date_list+=[x.prof.date._list[-1]]
		return recent_date_list	

	@staticmethod
	def get_deployment_date_list():
		recent_date_list = []
		for x in BaseRead.all_dict.values():
			recent_date_list+=[x.prof.date._list[0]]
		return recent_date_list	

	@staticmethod
	def get_recent_bins(lat_bins,lon_bins):
		bin_list = []
		for x in BaseRead.all_dict.values():
			bin_list+=[x.prof.pos.return_pos_bins(lat_bins,lon_bins,index_values=False)[-1]]
		return bin_list

	@staticmethod
	def get_recent_pos():
		pos_list = [x.prof.pos._list[-1] for x in BaseRead.all_dict.values()]
		return pos_list

	@staticmethod
	def get_floats_in_box(shape):
		float_name_list = []
		for dummy_float in BaseRead.all_dict.values():
			truth_list = [shape.contains(x) for x in GeoList(dummy_float.prof.pos._list).to_shapely()]
			if any(truth_list):
				float_name_list.append(dummy_float.meta.id)
		return float_name_list

	@staticmethod
	def get_subsampled_float_dict(percent):
		N = round(len(BaseRead.all_dict)*percent)
		return dict(random.sample(BaseRead.all_dict.items(), N))

	@staticmethod
	def get_sensors():
		sensor_list = [x.get_variables() for x in BaseRead.all_dict.values()]
		return sensor_list

class Speed():
	""" class that allows dynamic logic of trajectory speed 

		Parameters
		----------
		pos: position class
		date: date class
		speed_limit: limit speed that will create a problem
	"""        
	name = 'Speed'
	def __init__(self,date,pos,speed_limit):

		dist_list = []
		for idx in range(len(pos._list)-1):
			start_pos = pos._list[idx]
			end_pos = pos._list[idx+1]
			diff = geopy.distance.great_circle(end_pos,start_pos).nm
			dist_list.append(diff)
		time_list = [date._list[idx+1]-date._list[idx] for idx in range(len(date._list)-1)]
		time_list = [(_.days)*24.+(_.seconds)/3600. for _ in time_list]
		self._list = np.array(dist_list)/np.array(time_list)
		self._speed_limit = speed_limit

	def is_problem(self):
		"""
		Returns boolean value to show that all the date tests have been passed.

		Current tests:
		All speeds less than the speed limit (default is 5 nm/hr)

		Returns
		-------
		Boolean value to show if this profile is a problem
		"""
		bool_value = (np.array(self._list)>self._speed_limit).any()
		if bool_value:
			print('speed is the problem')
		return bool_value

	def masked_speed(self):
		""" checks speeds to make sure they are sensible. Current checks require 
		that all velocities are positive and less than 5 nm/hr.

			Parameters
			----------
			list of speeds in nm/hour

			Returns
			-------
			checked list of speeds in nm/hour
		"""
		speed_array = np.ma.masked_greater(self._list,self._speed_limit)
		speed_array = np.ma.masked_less(self.array,0)
		speed = speed.tolist()+[None]
		return speed

class all_drift_holder():
	def __init__(self,drift_class,lat_spacing=2,lon_spacing=2):
		self.class_list = drift_class.all_dict
		self.lat_bins = np.arange(-90,90,lat_spacing)
		self.lon_bins = np.arange(-180,180,lon_spacing)

	def full_time_delta_parse(self,span):
		for _ in self.class_list.values():
			_.prof.time_delta_parse(span)


	def get_full_bins(self):
		bin_list = []
		for _ in self.class_list.values():
			bin_list+=_.prof.pos.return_pos_bins(self.lat_bins,self.lon_bins,index_return=False)
		return bin_list


	def data_to_gridded(self,data,apply_function):
		"""
		function that calculates applies generic function (ie. np.nanmean) to a general set of 
		class data (ie. self.get_full_speed_list()).

		Parameters
		----------
		data: any full data list from self.class_list variable
		apply_function: any generic apply function

		Returns
		-------
		gridded data (according to prescribed lat_spacing, lon_spacing) with function applied to every cell
		"""
		data = np.array(data)
		data[data==None]=np.nan
		data = data.astype(float)
		XX,YY = np.meshgrid(self.lon_bins,self.lat_bins)
		bin_list  = self.get_full_bins()
		lat_list,lon_list = zip(*bin_list)
		lon_list = np.array(lon_list)
		lat_list = np.array(lat_list)
		data_list = []
		for k,dummy_bin in enumerate(zip(XX.ravel(),YY.ravel())):
			mask = (lon_list==dummy_bin[0])&(lat_list==dummy_bin[1])
			dummy_data = data[mask]
			data_list.append(apply_function(dummy_data))
		gridded = np.array(data_list).reshape(XX.shape)
		return gridded