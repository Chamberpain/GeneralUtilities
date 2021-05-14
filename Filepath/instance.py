import os
from datetime import date
import numpy as np 


def get_base_folder():
	file_path = '/Users/pchamberlain/Data/'
	if os.path.exists(file_path):
		pass
	else:
		file_path = '/home/pchamber/local_data/'
	return file_path

def does_file_exist(filename,mat):
	if os.path.isfile(filename+'.npy'):
		print(filename+'.npy')
		print('does exist, I will not remake the covariance')
		pass
	else:
		np.save(filename,mat)
		print(filename+'.npy')
		print('does not exist, I will make the covariance')

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
			os.makedirs(base_folder)
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


# def return_day_filepath(self):
# 	return date.today().strftime("%d-%m-%y")

# class KalmanNameUtilities(NameUtilities):
# 	project = 'kalman_smoother'

# 	def __init__(self,subproject=None):
# 		super().__init__(subproject = subproject)

# 	def return_day_filepath(self):
# 		today_string = super().return_day_filepath()
# 		print(today_string)