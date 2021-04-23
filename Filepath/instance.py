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
	def __init__(self,init_root_dir):
		pipeline_base = init_root_dir.replace('Utilities','Pipeline')
		self.tmp = os.path.join(pipeline_base,'tmp')
		self.out = os.path.join(pipeline_base,'out')
		self.store = os.path.join(pipeline_base,'store')

	def file_return(self,base_folder,filename):
		return os.path.join(base_folder,filename)

	def tmp_file(self,filename):
		self.file_return(self.tmp,filename)

	def out_file(self,filename):
		self.file_return(self.out,filename)

	def store_file(self,filename):
		self.file_return(self.store,filename)

# def return_day_filepath(self):
# 	return date.today().strftime("%d-%m-%y")

# class KalmanNameUtilities(NameUtilities):
# 	project = 'kalman_smoother'

# 	def __init__(self,subproject=None):
# 		super().__init__(subproject = subproject)

# 	def return_day_filepath(self):
# 		today_string = super().return_day_filepath()
# 		print(today_string)