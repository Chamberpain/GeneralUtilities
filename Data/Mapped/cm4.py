from OptimalArray.Utilities.CM4Mat import CovCM4Global
from GeneralUtilities.Data.Mapped.mapped_base import MappedBase
from GeneralUtilities.Compute.list import LonList, LatList
import matplotlib.pyplot as plt
import datetime
import os
import numpy as np
from netCDF4 import Dataset
from TransitionMatrix.Utilities.Utilities import shiftgrid
from OptimalArray.Utilities.CM4Mat import CovCM4
import pickle
from GeneralUtilities.Data.Filepath.instance import FilePathHandler, make_folder_if_does_not_exist
from GeneralUtilities.Data.Mapped.__init__ import ROOT_DIR

datadir = CovCM4.data_directory

class CM4(MappedBase):
	type = 'cm4'
	max_depth_lev = 25
	def return_dimensions(self):
		filename = dict([(x,y) for x,y in CovCM4Global.get_filenames()])[self.variable]
		nc_fid = Dataset(filename[0])
		lats = nc_fid["lat"][:]
		lons = nc_fid["lon"][:]
		data = nc_fid[self.variable][0, 0,:,:]
		data,lons = shiftgrid(180.5, data, lons, start=False)
		return (LonList(lons),LatList(lats))

	def return_units(self):
		filename = dict([(x,y) for x,y in CovCM4Global.get_filenames()])[self.variable]
		nc_fid = Dataset(filename[0])
		return nc_fid.variables[self.variable].units

	def return_dataset(self,depth_idx=2):
		master_list = CovCM4Global.get_filenames()
		file_list = dict([(x,y) for x,y in master_list])[self.variable]
		array_variable_list = []
		for file in np.sort(file_list):
			nc_fid = Dataset(file)
			array_variable_list.append(nc_fid[self.variable][:, depth_idx,:,:])
		lons = nc_fid["lon"][:]
		data = np.vstack(array_variable_list)
		data,lons = shiftgrid(180.5, data, lons, start=False)
		data = np.ma.masked_greater(data, 10 ** 19)
		data.units = nc_fid[self.variable].units
		return (data,lons)

	def return_int(self):
		master_list = CovCM4Global.get_filenames()
		file_list = dict([(x,y) for x,y in master_list])[self.variable]
		array_variable_list = []
		for file in np.sort(file_list):
			nc_fid = Dataset(file)
			var_temp = nc_fid[self.variable][:,:self.max_depth_lev,:,:]
			if file ==np.sort(file_list)[0]:
				delta_z = np.diff(CovCM4Global.get_depths().data)
				delta_z = delta_z.flatten().tolist()[:self.max_depth_lev]
				reshape_dims = np.ones(var_temp.shape, dtype=int)
				for k,delta in enumerate(delta_z):
					reshape_dims[:,k,:,:] = reshape_dims[:,k,:,:]*delta
			print('The data mask sum is ',var_temp.mask.sum())
			var_temp = var_temp*reshape_dims[:var_temp.shape[0],:,:,:]
			var_temp = var_temp.sum(axis=1) #now we have computed the vertical integral
			array_variable_list.append(var_temp)
		lons = nc_fid["lon"][:]
		data = np.vstack(array_variable_list)
		data,lons = shiftgrid(180.5, data, lons, start=False)
		data = np.ma.masked_greater(data, 10 ** 19)
		data.units = nc_fid[self.variable].units
		return (data,lons)

	def return_mean(self,depth_idx=2):
		filename = os.path.join(datadir,self.type+'_'+self.variable+'_'+str(depth_idx)+'_mean.pkl')
		try:
		    with open(filename, 'rb') as file:
		    	data = pickle.load(file)
		    print('The mean was successfully loaded')
		    file.close()
		except:
			print('The mean has not been calculated. Calculating...')
			data = self.return_dataset(depth_idx=depth_idx)[0].mean(axis=0)
			with open(filename, 'wb') as file:
				pickle.dump(data, file)
			file.close()
			print('The Mean was successfully calculated and saved. ')
		return data

	def return_var(self,depth_idx=2):
		filename = os.path.join(datadir,self.type+'_'+self.variable+'_'+str(depth_idx)+'_var.pkl')
		try:
		    with open(filename, 'rb') as file:
		    	data = pickle.load(file)
		    print('The variance was successfully loaded')
		    file.close()
		except:
			print('The variance has not been calculated. Calculating...')
			data = np.nanvar(self.return_dataset(depth_idx=depth_idx)[0],axis=0)
			with open(filename, 'wb') as file:
				pickle.dump(data, file)
			file.close()
			print('The variance was successfully calculated and saved. ')
		return data

	def return_int_mean(self):
		filename = os.path.join(datadir,self.type+'_'+self.variable+'_int_mean.pkl')
		try:
		    with open(filename, 'rb') as file:
		    	data = pickle.load(file)
		    print('The mean was successfully loaded')
		    file.close()
		except:
			print('The mean has not been calculated. Calculating...')
			data = self.return_int()[0].mean(axis=0)
			with open(filename, 'wb') as file:
				pickle.dump(data, file)
			file.close()
			print('The Mean was successfully calculated and saved. ')
		return data

	def return_int_var(self):
		filename = os.path.join(datadir,self.type+'_'+self.variable+'_int_var.pkl')
		try:
		    with open(filename, 'rb') as file:
		    	data = pickle.load(file)
		    print('The variance was successfully loaded')
		    file.close()
		except:
			print('The variance has not been calculated. Calculating...')
			data = np.nanvar(self.return_int()[0],axis=0)
			with open(filename, 'wb') as file:
				pickle.dump(data, file)
			file.close()
			print('The variance was successfully calculated and saved. ')
		return data

class CM4O2(CM4):
	variable = 'o2'
	plot_file_handler = FilePathHandler(ROOT_DIR,CM4.type+'/'+variable+'/')

class CM4PO4(CM4):
	variable = 'po4'
	plot_file_handler = FilePathHandler(ROOT_DIR,CM4.type+'/'+variable+'/')

class CM4ThetaO(CM4):
	variable = 'thetao'
	plot_file_handler = FilePathHandler(ROOT_DIR,CM4.type+'/'+variable+'/')

class CM4Sal(CM4):
	variable = 'so'
	plot_file_handler = FilePathHandler(ROOT_DIR,CM4.type+'/'+variable+'/')

class CM4PH(CM4):
	variable = 'ph'
	plot_file_handler = FilePathHandler(ROOT_DIR,CM4.type+'/'+variable+'/')

class CM4CHL(CM4):
	variable = 'chl'
	plot_file_handler = FilePathHandler(ROOT_DIR,CM4.type+'/'+variable+'/')

class CM4DIC(CM4):
	variable = 'dissic'
	plot_file_handler = FilePathHandler(ROOT_DIR,CM4.type+'/'+variable+'/')


class SurfaceCM4(MappedBase):
	type = 'cm4'

	def return_units(self):
		filename = dict([(x,y) for x,y in CovCM4Global.get_filenames()])[self.variable]
		nc_fid = Dataset(filename[0])
		return nc_fid.variables[self.variable].units

	def return_dataset(self):
		master_list = CovCM4Global.get_filenames()
		file_list = dict([(x,y) for x,y in master_list])[self.variable]
		array_variable_list = []
		for file in np.sort(file_list):
			nc_fid = Dataset(file)
			array_variable_list.append(nc_fid[self.variable][:,:,:])
		lons = nc_fid["lon"][:]
		data = np.vstack(array_variable_list)
		data,lons = shiftgrid(180.5, data, lons, start=False)
		data = np.ma.masked_greater(data, 10 ** 19)
		data.units = nc_fid[self.variable].units
		return (data,lons)

	def return_mean(self):
		filename = os.path.join(datadir,self.type+'_'+self.variable+'_'+str(depth_idx)+'_mean.pkl')
		try:
		    with open(filename, 'rb') as file:
		    	data = pickle.load(file)
		    print('The mean was successfully loaded')
		    file.close()
		except:
			print('The mean has not been calculated. Calculating...')
			data = self.return_dataset()[0].mean(axis=0)
			with open(filename, 'wb') as file:
				pickle.dump(data, file)
			file.close()
			print('The Mean was successfully calculated and saved. ')
		return data

	def return_var(self):
		filename = os.path.join(datadir,self.type+'_'+self.variable+'_'+str(depth_idx)+'_var.pkl')
		try:
		    with open(filename, 'rb') as file:
		    	data = pickle.load(file)
		    print('The variance was successfully loaded')
		    file.close()
		except:
			print('The variance has not been calculated. Calculating...')
			data = np.nanvar(self.return_dataset()[0],axis=0)
			with open(filename, 'wb') as file:
				pickle.dump(data, file)
			file.close()
			print('The variance was successfully calculated and saved. ')
		return data

class CM4PC02(SurfaceCM4):
	variable = 'spco2'
	plot_file_handler = FilePathHandler(ROOT_DIR,SurfaceCM4.type+'/'+variable+'/')

