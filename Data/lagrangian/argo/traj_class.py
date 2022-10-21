from GeneralUtilities.Data.Lagrangian.Argo.utilities import BaseReadClass,ArgoTime,Speed,ArgoTime
from netCDF4 import Dataset

class TrajDict(dict):
	pass

traj_class_dict = TrajDict({})

class Date(ArgoTime):
	@classmethod
	def from_ncfid(cls,nc_fid,*args,**kwargs):
		try:
			date_data = data_adjust(nc_fid['JULD'],nc_fid['JULD_ADJUSTED'])
		except IndexError:
			date_data = nc_fid['JULD'][:]
		mask = (np.array([_ in ['1','2'] for _ in format_byte_list_to_string(nc_fid['POSITION_QC'][:].tolist())]))&\
		(~nc_fid['LONGITUDE'][:].mask)&\
		(~nc_fid['LATITUDE'][:].mask)&\
		(~nc_fid['POSITION_ACCURACY'][:].mask)&\
		(~date_data.mask)
		data_list = julian_time_parse(date_data[mask].data, parse_time(nc_fid['REFERENCE_DATE_TIME']))
		return (cls(data_list),mask)

class BaseTrajClass(BaseReadClass):
	""" class to organize all of read in trajectory net cdf data
		----------
		file_path: the file path of the net cdf file you wish to read

	"""	
	name = 'Traj'
	def __init__(self,file_path):
		print('I am opening Traj File')
		nc_fid = Dataset(file_path)
		date_class,mask = Date.from_ncfid(nc_fid)
		self.date = date_class
		self.position_accuracy = nc_fid['POSITION_ACCURACY'][mask]
		self.pos = self.Position.from_ncfid(nc_fid,mask)
		self.speed = Speed.from_pos_and_time_list(self.pos,self.date)
		nc_fid.close()



	# def assign_depth(self,depth):
	# 	try:
	# 		self.depth = [depth.return_z(_) for _ in self.pos]