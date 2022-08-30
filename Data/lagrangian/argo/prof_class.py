import datetime
from netCDF4 import Dataset
import numpy as np
from GeneralUtilities.Data.Lagrangian.Argo.utilities import BaseReadClass,ArgoTime,Speed,ArgoTime,format_byte_list_to_string
from GeneralUtilities.Data.Lagrangian.lagrangian_utilities import julian_time_parse,parse_time


prof_class_dict = {}

class Date(ArgoTime):
	@classmethod
	def from_ncfid(cls,nc_fid,*args,**kwargs):
		date_data = nc_fid['JULD_LOCATION'][:]
		mask = (np.array([_ in ['1','2'] for _ in format_byte_list_to_string(nc_fid['POSITION_QC'][:].tolist())]))&\
		(~nc_fid['LONGITUDE'][:].mask)&\
		(~nc_fid['LATITUDE'][:].mask)&\
		(~date_data.mask)
		
		data_list = julian_time_parse(date_data[mask].data, parse_time(nc_fid['REFERENCE_DATE_TIME']))
		return (cls(data_list),mask)


class BaseProfClass(BaseReadClass):
	""" class to organize all of read in profile net cdf data
		----------
		file_path: the file path of the net cdf file you wish to read
	"""		
	name = 'Profile'
	def __init__(self,file_path):
		print('I am opening Prof File')
		nc_fid = Dataset(file_path)
		date_class,mask = Date.from_ncfid(nc_fid)
		self.date = date_class
		self.pos = self.Position.from_ncfid(nc_fid,mask)
		self.speed = Speed.from_pos_and_time_list(self.pos,self.date)
		nc_fid.close()


