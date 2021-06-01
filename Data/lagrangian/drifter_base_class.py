from __future__ import print_function
import geopy
import geopy.distance
import numpy as np
import matplotlib.pyplot as plt
import numpy as np
import random

class BasePosition(object):
	name = 'Position'

	def plot(self):
		lat,lon = zip(*[(_.latitude,_.longitude) for _ in self._list])
		plt.plot(lon,lat)

	def return_pos_bins(self,lat_bins,lon_bins,index_values=None):
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
		lat_bin_index = np.digitize([x.latitude for x in pos_list],lat_bins,right=True)
		lon_bin_index = np.digitize([x.longitude for x in pos_list],lon_bins,right=True)
		return [geopy.Point(x) for x in zip(lat_bins[lat_bin_index],lon_bins[lon_bin_index])]

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
		test_2 = (np.array(difference_list)>0).all()
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

	def get_pos_list(self):
		argos_list = []
		gps_list = []
		for _ in self.all_dict.iteritems():
			if _[1].meta.positioning_system =='ARGOS':
				argos_list.append(_[0])
			if _[1].meta.positioning_system =='GPS':
				gps_list.append(_[0])
		return (argos_list,gps_list)


	def get_full_speed_list(self):
		full_speed_list = []

		for k,_ in enumerate(self.all_dict.itervalues()):
			speed = _.prof.speed._list
			full_speed_list+=speed.tolist()
		return (full_speed_list)


	def get_full_lat_lon_list(self):
		lat_list = []
		lon_list = []
		for _ in self.all_dict.itervalues():
			lat,lon = zip(*[(dummy.latitude,dummy.longitude) for dummy in _.prof.pos._list])
			lat_list += lat
			lon_list += lon
		return (lat_list,lon_list)

	def get_full_date_list(self):
		full_date_list = []
		for _ in self.all_dict.itervalues():
			full_date_list+=_.prof.date._list
		return full_date_list

	@staticmethod
	def get_subsampled_float_dict(percent):
		N = round(len(BaseRead.all_dict)*percent)
		return dict(random.sample(BaseRead.all_dict.items(), N))

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
		for _ in self.class_list.itervalues():
			_.prof.time_delta_parse(span)


	def get_full_bins(self):
		bin_list = []
		for _ in self.class_list.itervalues():
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