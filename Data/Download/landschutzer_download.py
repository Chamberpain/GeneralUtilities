from GeneralUtilities.Data.Filepath.instance import get_data_folder
from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
data_folder = os.path.join(get_data_folder(),'Processed/Landschutzer/')

def download():
	url = 'https://www.nodc.noaa.gov/archive/arc0105/0160558/3.3/data/0-data/spco2_1982-2015_MPI_SOM-FFN_v2016.nc'
	filename = url.split("/")[-1]
	execute_download(url,filename,data_folder)