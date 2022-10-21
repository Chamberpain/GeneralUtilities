from GeneralUtilities.Data.Lagrangian.Argo.utilities import BaseReadClass,ArgoTime,Speed,ArgoTime,format_byte_list_to_string
from netCDF4 import Dataset
from GeneralUtilities.Compute.list import DepthList

class TechDict(dict):
	pass

tech_class_dict = TechDict({})

class BaseTechClass(BaseReadClass):
	name = 'Tech'
	def __init__(self,file_path):
		print('I am opening Tech File')
		nc_fid = Dataset(file_path)
		self.drift_depth = DriftDepth.from_ncfid(nc_fid)


class DriftDepth(DepthList):
	@classmethod
	def from_ncfid(cls,nc_fid):
		variable_name_list = [''.join(format_byte_list_to_string(i)).replace(" ","") for i in nc_fid['TECHNICAL_PARAMETER_NAME'][:].data]
		
		indices = [] 
		variable_list = ['PRES_ParkMaximum_dBAR','PRES_ParkMinimum_dBAR']
		indices += [i for i, x in enumerate(variable_name_list) if x in variable_list]
		# if not indices:
		# 	print(np.unique(variable_name_list))

		holder = [''.join(format_byte_list_to_string(nc_fid['TECHNICAL_PARAMETER_VALUE'][:].data[_,:])).replace(" ","") for _ in indices]
		data_list = []
		for dummy in holder:
			data_list.append(int(''.join([_ for _ in dummy if _.isdigit()])))
		if data_list:
			return cls(data_list)
		else:
			return None

	def is_problem(self):
		print(self)
		if not self:
			return False
		test_1 = sum(np.array(self)<500)<3	# allow a window of 500 on either side of drift depth
		if not test_1:
			print('drift depth less than 500 mb')
		test_2 = sum(np.array(self)>1500)<3 # 
		if not test_2:
			print('drift depth greater than 1500 mb')
		return ~(test_1&test_2)

	# def assign_depth(self,depth):
	# 	try:
	# 		self.depth = [depth.return_z(_) for _ in self.pos]
