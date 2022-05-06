from GeneralUtilities.Data.Filepath.instance import get_data_folder
from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
data_folder = os.path.join(get_data_folder(),'Raw/Modis/')

def download():
	for year in range(2003,2018):
		for month in range(1,13):
			url = 'https://jeodpp.jrc.ec.europa.eu/ftp/public/JRC-OpenData/GMIS/satellite/9km/GMIS_A_BBP_%02d_%4d.nc' % (month,year)
			filename = url.split("/")[-1]
			execute_download(url,filename,data_folder)