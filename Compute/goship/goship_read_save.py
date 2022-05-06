import sys,os
sys.path.append(os.path.abspath("../"))
import soccom_proj_settings

import pandas as pd
import fnmatch
from netCDF4 import Dataset
import numpy as np
import time
from datetime import datetime

matches = []
frames = []
variable_list = ['pressure','salinity','temperature','oxygen','bottle_salinity','bottle_oxygen','nitrate','alkalinity','total_alkalinity'] #ph is temperature dependent, and i am too lazy to do the correction

def goship_file_reader(file_):
	data_type = ' '+match.split('/')[-2].split('_')[-1]
	nc_fid = Dataset(file_, 'r')
	df_token = pd.DataFrame()
	for variable in variable_list:
		try: 
			df_token[variable.capitalize()] = nc_fid[variable][:]
			df_token[variable+'_QC'] = nc_fid[variable+'_QC'][:]
			df_token.loc[df_token[variable+'_QC']!=2,variable.capitalize()]=np.nan #2 is the "acceptable" flag in goship files

		except IndexError: # this is the case where the requested variable is not in the goship data
			continue

	df_token['Lat']=nc_fid['latitude'][:][0]
	df_token['Lon']=nc_fid['longitude'][:][0]

	try:
		cast_number = ''.join(nc_fid['cast'][:])
	except TypeError:
		cast_number = str(nc_fid['cast'][0])

	try: 
		df_token['Date']=datetime.strptime(str(nc_fid['woce_date'][:][0]),'%Y%m%d')
	except ValueError:
		df_token['Date'] = np.nan #there are a few profile files without reasonable dates
	except IndexError:
		df_token['Date'] = pd.to_datetime(datetime.strptime(nc_fid.Created,'%Y%m%d %I:%S').date()) # not all the net cdf files are in the same format, this normalizes them
	try:
		df_token['Cruise']=nc_fid.WOCE_ID+data_type+' '+cast_number
		if nc_fid.WOCE_ID in ['','999','UNK','UNKNOWN','x','TEST']: #not all cruise metadata is as complete as we would like
			return pd.DataFrame()
	except AttributeError:
		df_token['Cruise']=nc_fid.WHP_SECTION_ID+data_type+' '+cast_number # not all net cdf files have the same naming convention, this references the correct one
		if nc_fid.WHP_SECTION_ID in ['','999','UNK','UNKNOWN','x','TEST']: #not all cruise metadata is as complete as we would like
			return pd.DataFrame()
	df_token = df_token.dropna(subset = ['Date','Pressure']) #without date or pressure, these measurements are useless...
	df_token['Cruise']=n
	nc_fid.close()
	return df_token



for root, dirnames, filenames in os.walk(soccom_proj_settings.goship_data_directory):
    for filename in fnmatch.filter(filenames, '*.nc'): #find all net cdf files in the specified directory
        matches.append(os.path.join(root, filename))

for n, match in enumerate(matches):
	if (n % 100) ==0:   # print every 100th value, so we know where we are and how long it will take 
		print 'file is ',match,', there are ',len(matches[:])-n,'goship profiles left'
	frames.append(goship_file_reader(match))
df_holder = pd.concat(frames) #merge all dataframes together
df_holder = df_holder[df_holder.Lat<-40]
df_holder = df_holder[df_holder.Lat>-90] # some null lat values are saved as -999. This eliminates them.
df_holder = df_holder[df_holder.Pressure>0] # some null Pressure values are saved as -999. This eliminates them.
df_holder['Type']='GOSHIP'

df_holder.loc[df_holder.Alkalinity.isnull(),'Alkalinity']=df_holder[df_holder.Alkalinity.isnull()].Total_alkalinity #need to verify that total_alkalinity and alkalinity are the same variable


df_holder = df_holder[['Cruise','Date','Temperature','Salinity','Oxygen','Nitrate','Alkalinity','Pressure','Lat','Lon','Type']]
df_holder = df_holder.sort_values(['Cruise','Date','Pressure'])    #sort in a reasonable manner
df_holder = df_holder.reset_index(drop=True)
df_holder.loc[df_holder.Oxygen<60,'Oxygen']=np.nan #oxygen = 0 is the nan case for a subset of these data.
# df_holder['Cruise'] = [x.strip() for x in df_holder.Cruise] 
df_holder.to_pickle(soccom_proj_settings.goship_file)