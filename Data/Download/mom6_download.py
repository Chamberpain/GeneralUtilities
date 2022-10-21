from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
from OptimalArray.Utilities.MOM6Mat import CovMOM6

def download():
	for month_idx in range(1,13):
		url = 'https://zenodo.org/record/6492027/files/MOM6_CCS_bgc_phys_2008_%02d.nc' % (month_idx) 
		filename = url.split("/")[-1]
		execute_download(url,filename,CovMOM6.data_directory)