from GeneralUtilities.Data.Filepath.instance import get_data_folder
from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
import tarfile

data_folder = os.path.join(get_data_folder(),'Processed/AGVA/')

def download():
	for year in range(2004,2011):
		url = 'http://krill.ocean.washington.edu/agva/data/AGVA_v1.2_vel_%4d.tar.gz' % (year)
		filenamegz = url.split("/")[-1]
		filename = ".".join([filenamegz.split(".")[0],'nc'])
		if execute_download(url,filename,data_folder):
			os.rename(os.path.join(data_folder,filename),os.path.join(data_folder,filenamegz))
			file = tarfile.open(os.path.join(data_folder,filenamegz))
			file.extractall(data_folder)
			file.close()
			os.remove(os.path.join(data_folder,filenamegz))