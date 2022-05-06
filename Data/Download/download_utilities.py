import requests
import os
from GeneralUtilities.Data.Filepath.instance import make_folder_if_does_not_exist

def execute_download(url,filename,data_folder,verify=True):
	make_folder_if_does_not_exist(data_folder)
	filepath = os.path.join(data_folder,filename)
	if os.path.exists(filepath):
		print(filepath+' already exists, no download needed')
		pass
		return False
	else:
		print(filepath+' does not exist, downloading...')
		with open(filepath, "wb") as f:
			r = requests.get(url,verify=verify)
			f.write(r.content)
		return True