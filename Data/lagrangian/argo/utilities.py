import datetime
import geopy
from GeneralUtilities.Data.Lagrangian.drifter_base_class import BasePosition,BaseSpeed,BaseTime
from GeneralUtilities.Data.Filepath.instance import get_data_folder
import os
import re

def compile_classes():
	#this should probably use the find files function in search utility
	data_file_name = get_data_folder()+'/Data/Raw/Argo'
	matches = []
	for root, dirnames, filenames in os.walk(data_file_name):
		print(filenames)
		meta_match = re.compile('.*meta.nc') # all folders will have a meta file
		if any([file.endswith('meta.nc') for file in filenames]):
			matches.append(root)
	return matches


def data_adjust(data_instance, data_adjusted_instance):
	"""
	Function replaces adjusted data in the original data array

	Parameters
	----------
	data_instance: net cdf data instance
	adjusted: adjusted net cdf data instance

	Returns
	-------
	data processed in a usable form (string, array, float)
	"""
	masked_data_array = data_instance[:].data
	masked_adjusted_instance = data_adjusted_instance[:]
	masked_data_array[~masked_adjusted_instance.mask]=masked_adjusted_instance[~masked_adjusted_instance.mask]
	masked_data_array = np.ma.masked_equal(masked_data_array,data_instance._FillValue)
	return masked_data_array



def format_byte_list_to_string(byte_list):
	return ['None' if x is None else x.decode("utf-8") for x in byte_list]

class BaseReadClass(object):
	pass
	""" base class that contains functions used by both ProfClass and TrajClass
		----------
	"""	
	class Position(BasePosition):
		""" class that is a holder of the position information

		"""
		max_distance_between = 926 #this is 500nm
		@classmethod
		def from_ncfid(cls,nc_fid,mask):
			pos_data = [geopy.Point(_[0],_[1]) for _ in zip(nc_fid['LATITUDE'][mask],nc_fid['LONGITUDE'][mask])]
			pos_list = cls(pos_data) 
			if not pos_list:
				print('the position information was empty, to prove it, take a look at')
				print(nc_fid['LATITUDE'][:])
				print(nc_fid['LONGITUDE'][:])
			return pos_list

class ArgoTime(BaseTime):
	project_start = datetime.datetime(1998,1,1)
	maximum_time_between_data = 270 # days

class Speed(BaseSpeed):
	speed_limit = 2.57 # equivalent to 5kts 
