from GeneralUtilities.Data.Filepath.instance import get_data_folder
from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
import gzip
import shutil

data_folder = os.path.join(get_data_folder(),'Processed/Roemmich-Gilson/')


def download():
	urls = ['https://sio-argo.ucsd.edu/gilson/argo_climatology/RG_ArgoClim_Temperature_2019.nc.gz',
	'https://sio-argo.ucsd.edu/gilson/argo_climatology/RG_ArgoClim_Salinity_2019.nc.gz']
	for url in urls:
		filenamegz = url.split("/")[-1]
		filename = ".".join(filenamegz.split(".")[:2])
		if execute_download(url,filename,data_folder,verify=False):
			os.rename(os.path.join(data_folder,filename),os.path.join(data_folder,filenamegz))
			with gzip.open(os.path.join(data_folder,filenamegz), 'rb') as f_in:
			    with open(os.path.join(data_folder,filename), 'wb') as f_out:
			        shutil.copyfileobj(f_in, f_out)	
			os.remove(os.path.join(data_folder,filenamegz))