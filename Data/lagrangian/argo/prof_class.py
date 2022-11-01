import datetime
from netCDF4 import Dataset
import numpy as np
from GeneralUtilities.Data.Lagrangian.Argo.utilities import ProfDate,BaseReadClass,Speed,format_byte_list_to_string
from GeneralUtilities.Data.Lagrangian.lagrangian_utilities import julian_time_parse,parse_time,data_return

class ProfDict(dict):
	pass

prof_class_dict = ProfDict({})


class BaseProfClass(BaseReadClass):
	""" class to organize all of read in profile net cdf data
		----------
		file_path: the file path of the net cdf file you wish to read
	"""		
	name = 'Core'
	def __init__(self,file_path):
		print('I am opening Prof File')
		nc_fid = Dataset(file_path)
		date_class,mask = ProfDate.from_ncfid(nc_fid)
		self.date = date_class
		self.pos = self.Position.from_ncfid(nc_fid,mask)
		self.speed = Speed.from_pos_and_time_list(self.pos,self.date)
		self.variables = ['PRES', 'TEMP', 'PSAL']
		nc_fid.close()


class BGCProfClass(BaseReadClass):
	""" class to organize all of read in profile net cdf data
		----------
		file_path: the file path of the net cdf file you wish to read

	"""	
	name = 'BGC'
	def __init__(self,file_path):
		print('I am opening BGCProf File')
		nc_fid = Dataset(file_path)
		self.variables = [data_return(x) for x in nc_fid['STATION_PARAMETERS'][:][-1]]
		date_class,mask = ProfDate.from_ncfid(nc_fid)
		self.date = date_class
		self.pos = self.Position.from_ncfid(nc_fid,mask)
		self.speed = Speed.from_pos_and_time_list(self.pos,self.date)
		nc_fid.close()

