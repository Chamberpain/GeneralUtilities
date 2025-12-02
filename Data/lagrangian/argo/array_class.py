from GeneralUtilities.Data.Lagrangian.Argo.utilities import compile_classes
from GeneralUtilities.Data.Lagrangian.Argo.argo_read import ArgoReader
from GeneralUtilities.Data.Lagrangian.drifter_base_class import DrifterArrayBase
from netCDF4 import Dataset
import numpy as np
from OptimalArray.Utilities.H import HInstance, Float
from GeneralUtilities.Compute.list import VariableList
import geopy

class ArgoArray(DrifterArrayBase):
	def __init__(self,*args,**kwargs):
		super().__init__(*args, **kwargs)

	@classmethod
	def load(cls,num=99999):
		matches = compile_classes()
		number = 0
		holder = cls()
		for match in matches:
			if number>num:
				continue
			argo_instance = ArgoReader(match,prof=True,traj=False,tech=False)
			print('I am opening number match ',number)
			print ('filename is ',matches[number])
			if not argo_instance.problem:
				holder.update({(argo_instance.meta.id,argo_instance)})
			number += 1
		return holder 

	def get_positioning_system_list(self):
		pos_list = []
		for ii in self.values():
			pos_list.append(ii.meta.positioning_system)
		return pos_list

	def get_positioning_system_breakdown(self):
		argos_list = []
		gps_list = []
		for _ in self.items():
			if _[1].meta.positioning_system =='ARGOS':
				argos_list.append(_[0])
			if _[1].meta.positioning_system =='GPS':
				gps_list.append(_[0])
		return (argos_list,gps_list)

	def get_id_list(self):
		id_list = []
		for ii in self.values():
			id_list.append(ii.meta.id)
		return id_list

	def recent_bins_by_sensor(self,variable,lat_bins,lon_bins):
		bin_list = self.get_recent_bins(lat_bins,lon_bins)
		sensor_list = self.get_sensors()
		sensor_mask = [variable in x for x in sensor_list]
		date_mask = self.get_recent_mask()
		mask = np.array(sensor_mask)&np.array(date_mask)
		return np.array(bin_list)[mask]




