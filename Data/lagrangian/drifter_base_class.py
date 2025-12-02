from __future__ import print_function
import geopy
import geopy.distance
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import random
import datetime
from GeneralUtilities.Compute.list import GeoList,SpeedList,TimeList
from abc import ABC,abstractmethod

class BaseSpeed(SpeedList,ABC):
	""" class that allows dynamic logic of trajectory speed 

		Parameters
		----------
		pos: position class
		date: date class
		speed_limit: limit speed that will create a problem
	"""        
	name = 'Speed'

	def __init__(self,*args,**kwargs):
		super().__init__(*args, **kwargs)

	@property
	@abstractmethod
	def speed_limit(self):
		pass

	def is_problem(self):
		"""
		Returns boolean value to show that all the date tests have been passed.

		Current tests:
		All speeds less than the speed limit (default is 5 nm/hr)

		Returns
		-------
		Boolean value to show if this profile is a problem
		"""
		bool_value = (np.array(self)>self.speed_limit).any()
		if bool_value:
			print('speed is the problem')
			print('maximum value of speed was ',max(np.array(self)))
		return bool_value


class BasePosition(GeoList,ABC):
	name = 'Position'

	def __init__(self,*args,**kwargs):
		super().__init__(*args,**kwargs)

	def plot(self):
		lat,lon = self.lats_lons()
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
			pos_list = list(filter(None,np.array(self)[index_values].tolist()))
		else:
			pos_list = list(filter(None,self))
		lat_bin_index = lat_bins.digitize([x.latitude for x in pos_list])
		lon_bin_index = lon_bins.digitize([x.longitude for x in pos_list])
		lats = [lat_bins[x-1] for x in lat_bin_index]
		lons = [lon_bins[x-1] for x in lon_bin_index]
		return [geopy.Point(x) for x in zip(lats,lons)]

	@property
	@abstractmethod
	def max_distance_between(self):
		pass

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

		if not self:
			print('position list is empty')
			return True
		try:
			lat_list,lon_list = self.lats_lons()
		except AssertionError:
			#this is the case for lats or lons outside of geographic bounds
			return True
		difference_list = self.distance_between()
		test_1 = (np.array(difference_list)<self.max_distance_between).all()
		if not test_1:
			print('difference between positions is greater than '+str(self.max_distance_between)+' nm')				
		test_2 = sum(np.array(difference_list)==0.0)<6
		if not test_2:
			print('the difference in position is 0')
			print('print the number of 0 position occurances is ')
			print(sum(np.array(difference_list)==0.0))
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


class BaseTime(TimeList,ABC):
	name = 'Time'

	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)

	@property
	@abstractmethod
	def project_start(self):
		pass

	@property
	@abstractmethod
	def maximum_time_between_data(self):
		pass

	def is_problem(self):
		"""
		Returns boolean value to show that all the date tests have been passed.

		Current tests:
		All time differences are greater than 0
		All profiles happen before today 
		All profiles happen after the beginning of the drifter program
		The maximum time difference between profiles cannot be greater than specified interval months

		Returns
		-------
		Boolean value to show if this profile is a problem
		"""

		seconds_difference_list = self.seconds_difference()
		test_1 = (np.array(seconds_difference_list)>0.0).all()	# make sure all time is going in the right direction
		if not test_1:
			print('time is not going in the right direction')
		test_2 = ((np.array(self))<datetime.datetime.now()).all() # make sure all profiles happen before today
		if not test_2:
			print('a profile happened in the future')
		test_3 = ((np.array(self))>self.project_start).all() # make sure all profiles are after beginning of the drifter program
		if not test_2:
			print('a profile happened before project existed')
		test_4 = ((np.array([x/(3600*24) for x in seconds_difference_list])<self.maximum_time_between_data)).all() # make sure maximum time difference between profiles is less than 9 months
		if not test_2:
			print('maximum time between profiles is greater than '+str(self.maximum_time_between_data)+' days')
		return ~(test_1&test_2&test_3&test_4)



class BaseRead(object):
	def __init__(self):
		self.problem = self.is_problem()

	@property
	@abstractmethod
	def is_problem(self):
		pass
	
class DrifterArrayBase(dict):
	def __init__(self,*args,lat_spacing=2,lon_spacing=2,**kwargs):
		super().__init__(*args, **kwargs)
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

	def get_full_speed_list(self):
		full_speed_list = []
		for _ in self.values():
			full_speed_list+=_.prof.speed
		return (full_speed_list)

	def get_full_date_list(self):
		full_date_list = []
		for _ in self.values():
			full_date_list+=_.prof.date
		return full_date_list

	def get_recent_date_list(self):
		recent_date_list = []
		for x in self.values():
			recent_date_list+=[x.prof.date[-1]]
		return recent_date_list	

	def get_time_mask(self,time):
		time_mask = []
		for x in self.values():
			nearest_time = x.prof.date.find_nearest(time,test=False)
			time_mask.append(abs((nearest_time-time).days)<=90) #if the nearest time is outside of 3 months on either side, then it is not included 
		return time_mask

	def get_time_bins(self,lat_bins,lon_bins,time):
		bin_list = []
		for x in self.values():
			idx = x.prof.date.find_nearest(time, idx=True,test=False)
			bin_list+=[x.prof.pos.return_pos_bins(lat_bins,lon_bins,index_values=False)[idx]]
		return bin_list

	def get_recent_mask(self):
		date_list = self.get_recent_date_list()
		return [max(date_list)-datetime.timedelta(days=180)<x for x in date_list]

	def get_deployment_date_list(self):
		recent_date_list = []
		for x in self.values():
			recent_date_list+=[x.prof.date[0]]
		return recent_date_list	

	def get_recent_bins(self,lat_bins,lon_bins):
		bin_list = []
		for x in self.values():
			bin_list+=[x.prof.pos.return_pos_bins(lat_bins,lon_bins,index_values=False)[-1]]
		return bin_list

	def get_full_lat_lon_list(self):
		lat_list = []
		lon_list = []
		for ii in self.values():
			lat,lon = zip(*[(dummy.latitude,dummy.longitude) for dummy in ii.prof.pos])
			lat_list += [lat]
			lon_list += [lon]
		return (LatList(lat_list),LonList(lon_list))

	def get_recent_pos(self):
		pos_list = [x.prof.pos[-1] for x in self.values()]
		return pos_list

	def get_floats_in_box(self,shape):
		float_name_list = []
		for dummy_float in self.values():
			truth_list = [shape.contains(x) for x in GeoList(dummy_float.prof.pos).to_shapely()]
			if any(truth_list):
				float_name_list.append(dummy_float.meta.id)
		return float_name_list

	def get_subsampled_float_dict(self,percent):
		N = round(len(self)*percent)
		return dict(random.sample(self.items(), N))

	def get_sensors(self):
		sensor_list = [x.get_variables() for x in self.values()]
		return sensor_list