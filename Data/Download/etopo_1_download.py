from GeneralUtilities.Data.Filepath.instance import get_data_folder
from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
import gzip
import shutil

data_folder = os.path.join(get_data_folder(),'Raw/ETopo1/')


def download():
	urls = ['https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gdal.grd.gz',
	'https://www.ngdc.noaa.gov/mgg/global/relief/ETOPO1/data/bedrock/grid_registered/netcdf/ETOPO1_Bed_g_gmt4.grd.gz']
	for url in urls:
		filenamegz = url.split("/")[-1]
		filename = ".".join(filenamegz.split(".")[:2])
		if execute_download(url,filename,data_folder):
			os.rename(os.path.join(data_folder,filename),os.path.join(data_folder,filenamegz))
			with gzip.open(os.path.join(data_folder,filenamegz), 'rb') as f_in:
			    with open(os.path.join(data_folder,filename), 'wb') as f_out:
			        shutil.copyfileobj(f_in, f_out)	
			os.remove(os.path.join(data_folder,filenamegz))