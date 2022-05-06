import numpy as np

class MappedBase():
	def __init__(self):
		pass

	def compile_variance(self,depth_idx):
		data,lons = self.return_dataset(depth_idx)
		var_data = data.var(axis=0)
		var_data.units = data.units
		return var_data

	def compile_monthly_mean_and_variance(self,depth_idx,lat,lon):
		lons,lats = self.return_dimensions()
		lon_idx = lons.find_nearest(lon,idx=True)
		lat_idx = lats.find_nearest(lat,idx=True)
		data,lons = self.return_dataset(depth_idx)
		mean_list = []
		std_list = []
		for starting_idx in range(12):
			holder = data[starting_idx::12,lat_idx,lon_idx]
			mean_list.append(np.nanmean(holder))
			std_list.append(np.nanmean(holder))
		return (mean_list,std_list)