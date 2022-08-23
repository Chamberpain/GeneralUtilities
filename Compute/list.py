from __future__ import print_function
import datetime
import numpy as np
import geopy
import shapely.geometry


def flat_list(non_flat_list):
	flat_list = [item for sublist in non_flat_list for item in sublist]
	return flat_list

def find_nearest(items, pivot,test=True):
	nearest = min(items, key=lambda x: abs(x - pivot))
	if test:
		item_range = max(items)-min(items)
		assert (nearest-pivot)<0.1*item_range # only allow 10% extrapolation
	return	nearest

class BaseList(list):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)	

	def __add__(self,other):
		return self.__class__(list.__add__(self,other))

	def __mul__(self,other):
		return self.__class__(list.__mul__(self,other))

	def __getslice__(self, start, stop):
		return self.__class__(list.__getslice__(start, stop))

	def __getitem__(self, item):
		result = list.__getitem__(self, item)
		if type(item) is slice:
			return self.__class__(result)
		else:
			return result

	def find_nearest(self,value,test = True, idx = False):
		if idx:
			list_value = find_nearest(self,value,test=test)
			return self.index(list_value)
		else:
			return find_nearest(self,value,test=test)

	def digitize(self,dummy):
		return np.digitize(dummy,self,right=False)

class VariableList(BaseList):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		assert all([isinstance(x,str) for x in self]) 

class DepthList(BaseList):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		assert max(self)<=0
		assert min(self)>=-10984 # depth of marianas trench

class LonList(BaseList):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		assert max(self)<180
		assert min(self)>=-180

	def return_lon360(self):
		holder = np.array(self)
		holder[holder<0]=holder[holder<0]+360
		return holder

class LatList(BaseList):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		assert max(self)<=90
		assert min(self)>=-90

class GeoList(BaseList):
	def __init__(self, *args, lat_sep=None, lon_sep=None,**kwargs):
		super().__init__(*args, **kwargs)
		# total list must be composed of geopy.Points 
		if isinstance(self,GeoList):
			assert all([isinstance(x,geopy.Point) for x in self]) 
		self.lat_sep = lat_sep
		self.lon_sep = lon_sep

	def return_dimensions(self):
		lat_grid = LatList(np.arange(-90,90.01,self.lat_sep))
		lon_grid = LonList(np.arange(-180,179.99,self.lon_sep))
		return (lat_grid,lon_grid)

	def tuple_total_list(self):
		return [tuple(x)[:2] for x in self]

	def lats_lons(self):
		lats,lons = zip(*self.tuple_total_list())
		return (LatList(lats),LonList(lons))

	def unique_lats_lons(self):
		lats,lons = zip(*self.tuple_total_list())
		return (LatList(np.sort(np.unique(lats))),LonList(np.sort(np.unique(lons))))

	def reduced_res(self,idx,new_lat_sep,new_lon_sep):
		lat_ratio = self.lat_sep/new_lat_sep
		lon_ratio = self.lon_sep/new_lon_sep
		assert lat_ratio%1==0
		assert lon_ratio%1==0
		#latitude and longitude must be divisable by each other, or there will be holes
		lat = self[idx].latitude
		lon = self[idx].longitude
		new_lats = [lat+new_lat_sep*x for x in range(int(lat_ratio))]
		new_lons = [lon+new_lon_sep*x for x in range(int(lon_ratio))]
		lats,lons = np.meshgrid(new_lats,new_lons)
		return list(zip(lats.flatten(),lons.flatten()))

	def to_shapely(self):
		import geopandas as gp
		return gp.GeoSeries([shapely.geometry.Point(x.longitude,x.latitude) for x in self])

class TimeList(BaseList):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		assert all([isinstance(x,datetime.datetime) for x in self]) 

	@classmethod
	def time_list_from_seconds(cls,seconds_list,ref_date):
		return cls([ref_date + datetime.timedelta(seconds=x) for x in seconds_list])

	@classmethod
	def time_list_from_minutes(cls,minutes_list,ref_date):
		return cls([ref_date + datetime.timedelta(minutes=x) for x in minutes_list])

	@classmethod
	def time_list_from_hours(cls,hours_list,ref_date):
		return cls([ref_date + datetime.timedelta(hours=x) for x in hours_list])

	@classmethod
	def time_list_from_days(cls,days_list,ref_date):
		return cls([ref_date + datetime.timedelta(days=x) for x in days_list])

	def days_since(self,ref_date):
		time_delta_list = [x-ref_date for x in self]
		return [x.days for x in time_delta_list]

	def hours_since(self,ref_date):
		hours_in_day = 24
		seconds_in_hour = 3600
		time_delta_list = [x-ref_date for x in self]
		return [x.days*hours_in_day+x.seconds*seconds_in_hour for x in time_delta_list]

	def seconds_since(self,ref_date):
		seconds_in_day = 24*60*60
		time_delta_list = [x-ref_date for x in self]
		return [x.days*seconds_in_day+x.seconds for x in time_delta_list]