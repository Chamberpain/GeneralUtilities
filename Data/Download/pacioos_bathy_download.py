from GeneralUtilities.Data.Filepath.instance import get_data_folder
from GeneralUtilities.Data.Download.download_utilities import execute_download
import os

data_folder = os.path.join(get_data_folder(),'Raw/PacioosBathy/')


def download():
	url = 'https://pae-paha.pacioos.hawaii.edu/thredds/ncss/hmrg_bathytopo_1km_mhi?var=z&horizStride=1&addLatLon=true'
	filename = 'hmrg_bathytopo_1km_mhi.nc'
	execute_download(url,filename,data_folder)