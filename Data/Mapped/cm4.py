from GeneralUtilities.Data.Mapped.mapped_base import MappedBase
from GeneralUtilities.Compute.list import LonList, LatList
import matplotlib.pyplot as plt
import datetime
import os
import numpy as np
from netCDF4 import Dataset
from TransitionMatrix.Utilities.Utilities import shiftgrid
import pickle
from GeneralUtilities.Data.Filepath.instance import FilePathHandler, make_folder_if_does_not_exist
from GeneralUtilities.Data.Mapped.__init__ import ROOT_DIR
from GeneralUtilities.Data.Filepath.instance import get_data_folder
import gsw
import OptimalArray.Utilities.CM4Mat as CM4Mat

data_directory = os.path.join(get_data_folder(),'Processed/CM4/')

def get_filenames():
	master_dict = {}
	for file in os.listdir(data_directory):
		if '.DS_Store' in file:
			continue
		if 'mean.pkl' in file:
			continue
		if 'var.pkl' in file:
			continue
		if 'int.pkl' in file:
			continue
		filename = os.path.join(data_directory,file)
		filename = filename
		var = file.split('_')[0]
		try:
			master_dict[var].append(filename)
		except KeyError:
			master_dict[var] = [filename]
	for var in master_dict.keys():
		master_dict[var] = np.sort(master_dict[var]).tolist()
	return list(master_dict.items())


class CM4(MappedBase):
	type = 'cm4'
	max_depth_lev = 25

	def return_filenames(self):
		return dict([(x,y) for x,y in get_filenames()])[self.variable]

	def return_dimensions(self):
		filename = self.return_filenames
		nc_fid = Dataset(filename[0])
		lats = nc_fid["lat"][:]
		lons = nc_fid["lon"][:]
		data = nc_fid[self.variable][0, 0,:,:]
		data,lons = shiftgrid(180.5, data, lons, start=False)
		return (LonList(lons),LatList(lats))

	def return_units(self):
		filename = self.return_filenames
		nc_fid = Dataset(filename[0])
		return nc_fid.variables[self.variable].units

	def return_dataset(self,depth_idx=2):
		master_list = get_filenames()
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

	def return_full_dataset(self):
		master_list = get_filenames()
		file_list = dict([(x,y) for x,y in master_list])[self.variable]
		array_variable_list = []
		for file in np.sort(file_list):
			print(file)
			nc_fid = Dataset(file)
			array_variable_list.append(nc_fid[self.variable][:,:self.max_depth_lev,:,:])
		lons = nc_fid["lon"][:]
		data = np.vstack(array_variable_list)
		data,lons = shiftgrid(180.5, data, lons, start=False)
		data = np.ma.masked_greater(data, 10 ** 19)
		data.units = nc_fid[self.variable].units
		return (data,lons)

	def return_full_mean(self):
		filename = os.path.join(data_directory,self.type+'_'+self.variable+'_full_mean.pkl')
		try:
		    with open(filename, 'rb') as file:
		    	data = pickle.load(file)
		    print('The mean was successfully loaded')
		    file.close()
		except:
			print('The mean has not been calculated. Calculating...')
			data,lons = self.return_full_dataset()
			data = data.mean(axis=0)
			with open(filename, 'wb') as file:
				pickle.dump(data, file)
			file.close()
			print('The Mean was successfully calculated and saved. ')
		return data

	def return_int(self):
		master_list = get_filenames()
		file_list = dict([(x,y) for x,y in master_list])[self.variable]
		array_variable_list = []
		for file in np.sort(file_list):
			nc_fid = Dataset(file)
			var_temp = nc_fid[self.variable][:,:self.max_depth_lev,:,:]
			if file ==np.sort(file_list)[0]:
				delta_z = np.diff(CM4Mat.CovCM4Global.get_depths().data)
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
		filename = os.path.join(data_directory,self.type+'_'+self.variable+'_'+str(depth_idx)+'_mean.pkl')
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
		filename = os.path.join(data_directory,self.type+'_'+self.variable+'_'+str(depth_idx)+'_var.pkl')
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
		filename = os.path.join(data_directory,self.type+'_'+self.variable+'_int_mean.pkl')
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
		filename = os.path.join(data_directory,self.type+'_'+self.variable+'_int_var.pkl')
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

	def return_depth_covariance(self):
		filename = os.path.join(data_directory,self.type+'_'+self.variable+'_vert_cov.pkl')
		try:
		    with open(filename, 'rb') as file:
		    	cov = pickle.load(file)
		    print('The vertical covariance was successfully loaded')
		    file.close()
		except:
			print('The variance has not been calculated. Calculating...')
			data,lons = self.return_full_dataset()
			cov = np.ones([self.max_depth_lev,self.max_depth_lev,data.shape[2],data.shape[3]])
			for i in range(data.shape[2]):
				print(i)
				for j in range(data.shape[3]):
					cov[:,:,i,j] = np.ma.cov(data[:,:,i,j].transpose()).filled()
			cov = np.ma.masked_greater(cov,10**10)
			with open(filename, 'wb') as file:
				pickle.dump(cov, file)
			file.close()
			print('The variance was successfully calculated and saved. ')
		return cov

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

class CM4Density(CM4):
	variable = 'density'
	plot_file_handler = FilePathHandler(ROOT_DIR,CM4.type+'/'+variable+'/')

	def return_filenames(self):
		return dict([(x,y) for x,y in get_filenames()])['thetao']

	def return_dimensions(self):
		filename = self.return_filenames()
		nc_fid = Dataset(filename[0])
		lats = nc_fid["lat"][:]
		lons = nc_fid["lon"][:]
		data = nc_fid['thetao'][0, 0,:,:]
		data,lons = shiftgrid(180.5, data, lons, start=False)
		return (LonList(lons),LatList(lats))

	def return_units(self):
		filename = dict([(x,y) for x,y in get_filenames()])['thetao']
		nc_fid = Dataset(filename[0])
		return nc_fid.variables['thetao'].units

	def return_dataset(self,depth_idx=2):
		depth = CM4Mat.CovCM4Global.get_depths()[depth_idx][1]
		master_list = get_filenames()
		data_dict = {}
		for var in ['thetao','so']:
			file_list = dict([(x,y) for x,y in master_list])[var]
			array_variable_list = []
			for file in np.sort(file_list):
				nc_fid = Dataset(file)
				array_variable_list.append(nc_fid[var][:, depth_idx,:,:])
			lons = nc_fid["lon"][:]
			data = np.vstack(array_variable_list)
			data,lons = shiftgrid(180.5, data, lons, start=False)
			data = np.ma.masked_greater(data, 10 ** 19)
			data_dict[var]=data
		data = gsw.density.rho(data_dict['so'],data_dict['thetao'],depth)
		return (data,lons)

	def return_full_dataset(self):
		master_list = get_filenames()
		file_dict = dict([(x,y) for x,y in master_list])
		file_list = zip(file_dict['thetao'],file_dict['so'])
		array_variable_list = []

		for file_t,file_s in file_list:
			nc_fid_t = Dataset(file_t)
			var_t = nc_fid_t['thetao'][:,:self.max_depth_lev,:,:]
			nc_fid_s = Dataset(file_s)
			var_s = nc_fid_s['so'][:,:self.max_depth_lev,:,:]
			
			z = [x[0] for x in CM4Mat.CovCM4Global.get_depths().data][:self.max_depth_lev]
			reshape_dims_z = np.ones(var_t.shape, dtype=int)
			for k,z in enumerate(z):
				reshape_dims_z[:,k,:,:] = reshape_dims_z[:,k,:,:]*z

			print('The data mask sum is ',var_s.mask.sum())
			var_temp = gsw.density.rho(var_s,var_t,reshape_dims_z)


			array_variable_list.append(var_temp)
		lons = nc_fid_t["lon"][:]
		data = np.vstack(array_variable_list)
		data,lons = shiftgrid(180.5, data, lons, start=False)
		data = np.ma.masked_greater(data, 10 ** 19)
		return (data,lons)

	def return_int(self):
		data,lons = self.return_full_dataset()
		delta_z = np.diff(CM4Mat.CovCM4Global.get_depths().data)
		delta_z = delta_z.flatten().tolist()[:self.max_depth_lev]		
		reshape_dims_delta_z = np.ones(data.shape, dtype=int)

		for k,delta in enumerate(delta_z):
			reshape_dims_delta_z[:,k,:,:] = reshape_dims_delta_z[:,k,:,:]*delta
		data = data*reshape_dims_delta_z
		data = data.sum(axis=1) #now we have computed the vertical integral
		return (data,lons)


class CM4HeatContent(CM4):
	variable = 'heat_content'
	plot_file_handler = FilePathHandler(ROOT_DIR,CM4.type+'/'+variable+'/')

	def return_int(self):
		master_list = get_filenames()
		file_list = dict([(x,y) for x,y in master_list])[self.variable]
		array_variable_list = []
		for file in np.sort(file_list):
			nc_fid = Dataset(file)
			var_temp = nc_fid[self.variable][:,:self.max_depth_lev,:,:]
			if file ==np.sort(file_list)[0]:
				delta_z = np.diff(CM4Mat.CovCM4Global.get_depths().data)
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



class SurfaceCM4(MappedBase):
	type = 'cm4'

	def return_units(self):
		filename = dict([(x,y) for x,y in get_filenames()])[self.variable]
		nc_fid = Dataset(filename[0])
		return nc_fid.variables[self.variable].units

	def return_dataset(self):
		master_list = get_filenames()
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
		filename = os.path.join(data_directory,self.type+'_'+self.variable+'_'+str(depth_idx)+'_mean.pkl')
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
		filename = os.path.join(data_directory,self.type+'_'+self.variable+'_'+str(depth_idx)+'_var.pkl')
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

