import pickle

def load(file_path):
	"""
	Function loads data saved within the folder

	Parameters
	----------
	file_path: the file path of the data saved within the folder

	Returns
	-------
	unpickled data saved within the folder
	"""

	infile = open(file_path,'rb')
	dummy_file = pickle.load(infile)
	infile.close()
	return dummy_file

def save(file_path,data):
	"""
	Function saves data saved within the folder

	Parameters
	----------
	file_path: the file path of the data you wish to save within the folder
	data: the data you wish to save
	Returns
	-------
	none
	"""
	f = open(file_path,"wb")
	pickle.dump(data,f)
	f.close()

