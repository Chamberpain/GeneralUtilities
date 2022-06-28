from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
from OptimalArray.Utilities.CM4Mat import CovCM4


def download():
	start_date_list = [185001,187001,189001,191001,193001,195001,197001,199001,201001]
	end_date_list = [186912,188912,190912,192912,194912,196912,198912,200912,201412]
	var_list = ['chl','o2','ph','so','thetao']

	for start_date,end_date in zip(start_date_list,end_date_list):
		for var in var_list:
			url = 'http://esgdata.gfdl.noaa.gov/thredds/fileServer/gfdl_dataroot4/CMIP/NOAA-GFDL/GFDL-CM4/historical/r1i1p1f1/Omon/%s/gr/v20180701/%s_Omon_GFDL-CM4_historical_r1i1p1f1_gr_%6d-%6d.nc' % (var,var,start_date,end_date) 
			filename = url.split("/")[-1]
			execute_download(url,filename,CovCM4.data_directory)