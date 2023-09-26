from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
from OptimalArray.Utilities.MOM6Mat import CovMOM6GOM

def download():
	for year_idx in range(2008,2019):
			for month_idx in range(1,13):
				url = 'https://zenodo.org/record/7591445/files/MOM6_GoM_bgc_phys_%04d_%02d.nc' % (year_idx,month_idx) 
				filename = url.split("/")[-1]
				execute_download(url,filename,CovMOM6GOM.data_directory)