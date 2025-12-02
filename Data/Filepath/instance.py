import os
from datetime import date
import numpy as np 
from GeneralUtilities.__init__ import ROOT_DIR


def get_data_folder():
	return os.path.join(ROOT_DIR,'DataDir')

def get_compute_folder():
	return os.path.join(ROOT_DIR,'ComputeDir')

def make_folder_if_does_not_exist(folder_path):
	print('working on '+folder_path)
	if not os.path.exists(folder_path):
		os.makedirs(folder_path)
		print("The new directory is created!")		
	else:
		print("The directory already exists")	

class PathHandlerBase(object):
	def check_folder(self,folder):
		CHECK_FOLDER = os.path.isdir(folder)
		if not CHECK_FOLDER:
			os.makedirs(folder)
			print("created folder : ",folder)

	def file_return(self,base_folder,filename):
		self.check_folder(base_folder)
		return os.path.join(base_folder,filename)

	def tmp_file(self,filename):
		return self.file_return(self._tmp,filename)

	def out_file(self,filename):
		return self.file_return(self._out,filename)

	def store_file(self,filename):
		return self.file_return(self._store,filename)

	def data_file(self,filename):
		return self.file_return(self._data,filename)	

class FilePathHandler(PathHandlerBase):
	def __init__(self,init_root_dir,filename):
		init_root_dir_list = init_root_dir.split('/')
		idx = init_root_dir.split('/').index('Projects')
		dir_list = init_root_dir_list[idx+1:]
		try:
			dir_list.remove('Utilities')
		except:
			pass
		pipeline_base = os.path.join(get_data_folder(),'/'.join(dir_list))


		self._tmp = os.path.join(pipeline_base,filename,'tmp')
		self._out = os.path.join(pipeline_base,filename,'out')
		self._store = os.path.join(pipeline_base,filename,'store')

		self._data = '/'.join(init_root_dir.split('/')[:idx+2])+'/Data/'

		for dirs in [self._data,self._tmp,self._out,self._store]:
			self.check_folder(dirs)

class ComputePathHandler(PathHandlerBase):
	def __init__(self,init_root_dir,filename):
		init_root_dir_list = init_root_dir.split('/')
		idx = init_root_dir.split('/').index('Projects')
		dir_list = init_root_dir_list[idx+1:]
		try:
			dir_list.remove('Utilities')
		except:
			pass
		pipeline_base = os.path.join(get_compute_folder(),'/'.join(dir_list))


		self._tmp = os.path.join(pipeline_base,filename,'tmp')
		self._out = os.path.join(pipeline_base,filename,'out')
		self._store = os.path.join(pipeline_base,filename,'store')

		self._data = '/'.join(init_root_dir.split('/')[:idx+2])+'/Data/'

		for dirs in [self._data,self._tmp,self._out,self._store]:
			self.check_folder(dirs)				
