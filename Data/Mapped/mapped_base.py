import numpy as np
from GeneralUtilities.Plot.Cartopy.eulerian_plot import GlobalCartopy
import cartopy.crs as ccrs
from TransitionMatrix.Utilities.TransGeo import get_cmap

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

	def make_movie(self,depth_idx):
		data,dummy = self.return_dataset(depth_idx=depth_idx)
		lons,lats = self.return_dimensions()
		XX,YY = np.meshgrid(lons,lats)

		for k,data_holder in enumerate([data[x,:,:] for x in range(data.shape[0])]):
			data_holder
			fig = plt.figure(figsize=(14,11))
			ax = fig.add_subplot(1,1,1, projection=ccrs.PlateCarree())
			plot_holder = GlobalCartopy(ax=ax,adjustable=True)
			dummy,dummy,ax0 = plot_holder.get_map()
			ax0.pcolormesh(XX,YY,XX != np.nan,cmap=get_cmap(),alpha=0.7,zorder=-10,transform=ccrs.PlateCarree())
			pcm = ax0.pcolormesh(XX,YY,data_holder,cmap='YlOrBr',vmin = 0, vmax=10,transform=ccrs.PlateCarree())
			fig.colorbar(pcm,pad=-0.05,label=self.variable,location='bottom')
			plt.savefig(self.plot_file_handler.out_file(str(k)),bbox_inches='tight')
			plt.close('all')


