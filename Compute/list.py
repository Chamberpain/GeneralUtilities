from __future__ import print_function
import datetime
import geopy

def find_nearest(items, pivot,test=True):
	nearest = min(items, key=lambda x: abs(x - pivot))
	item_range = max(items)-min(items)
	if test:
		assert (nearest-pivot)<0.1*item_range # only allow 10% extrapolation
	return	nearest

def flat_list(non_flat_list):
	flat_list = [item for sublist in non_flat_list for item in sublist]
	return flat_list

class LonList(list):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		assert max(self)<180
		assert min(self)>-180

	def return_lon360(self):
		holder = np.array(self)
		holder[holder<0]=holder[holder<0]+360
		return holder

class LatList(list):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		assert max(self)<90
		assert min(self)>-90

class TupleList(list):
	def __init__(self, *args, **kwargs):
		assert all([isinstance(x,geopy.Point) for x in self]) 
		# total list must be composed of geopy.Points 
		self = list.__init__(*args, **kwargs)
		return self	

	def tuple_total_list(self):
		return [tuple(x)[:2] for x in self]

	def lats_lons(self):
		lats,lons = zip(*self.tuple_total_list())
		return (LatList(lats),LonList(lons))

class TimeList(list):
	def __init__(self, *args, **kwargs):
		super().__init__(*args, **kwargs)
		assert all([isinstance(x,datetime.datetime) for x in self]) 

	@classmethod
	def set_ref_date(cls,ref_date):
		cls.ref_date = ref_date

	@staticmethod
	def time_list_from_seconds(seconds_list):
		return TimeList([TimeList.ref_date + datetime.timedelta(seconds=x) for x in seconds_list])

	def closest_datetime(self,datetime_instance):
		return find_nearest(self,datetime_instance)

	def closest_index(self,datetime_instance):
		return self.index(self.closest_datetime(datetime_instance))

	def days_since(self):
		time_delta_list = [x-self.ref_date for x in self]
		return [x.days for x in time_delta_list]

	def seconds_since(self):
		seconds_in_day = 24*60*60
		time_delta_list = [x-self.ref_date for x in self]
		return [x.days*seconds_in_day+x.seconds for x in time_delta_list]