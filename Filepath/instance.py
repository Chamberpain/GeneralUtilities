import os
from datetime import date
import numpy as np 
from GeneralUtilities.__init__ import ROOT_DIR


def get_base_folder():
	return os.path.join(ROOT_DIR,'DataDir')

class FilePathHandler(object):
	def __init__(self,init_root_dir,filename):
		pipeline_base = init_root_dir.replace('Utilities','Pipeline')
		self._tmp = os.path.join(pipeline_base,filename,'tmp')
		self._out = os.path.join(pipeline_base,filename,'out')
		self._store = os.path.join(pipeline_base,filename,'store')

		folder_idx = init_root_dir.split('/').index('Projects')
		self._data = '/'.join(init_root_dir.split('/')[:folder_idx+2])+'/Data/'

		for dirs in [self._data,self._tmp,self._out,self._store]:
			self.check_folder(dirs)
						
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