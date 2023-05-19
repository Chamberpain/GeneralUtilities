from GeneralUtilities.Data.Filepath.instance import get_data_folder
from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
data_folder = os.path.join(get_data_folder(),'Raw/Aviso/')

def download():
        for year in range(2012, 2020):
                for month in range(1, 13):
                        url = f"https://coastwatch.noaa.gov/erddap/griddap/noaacwBLENDEDsshDaily.nc?sla%5B({year}-{month}-01):1:({year}-{month}-31T00:00:00Z)%5D%5B(-89.875):1:(89.875)%5D%5B(-179.875):1:(179.875)%5D"
                        filename = str(month) + ".nc"
                        folder = os.path.join(data_folder, str(year))
                        execute_download(url,filename,folder)
