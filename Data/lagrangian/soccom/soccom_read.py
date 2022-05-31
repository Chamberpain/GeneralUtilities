import pandas as pd
import glob
import numpy as np
from netCDF4 import Dataset 
import sys,os
sys.path.append(os.path.abspath("../"))
import soccom_proj_settings
variable_list = ['Pressure','Temperature','Salinity','Oxygen','OxygenSat','Nitrate','pH25C','TALK_MLR','pHinsitu','DIC_MLR']   #initialize the column names

def soccom_file_reader(file_token):
    frame_token = []
    nc_token = Dataset(file_token)
    cruise = ''.join(nc_token['Cruise'][:].tolist()).strip()
    for date_header in range(nc_token.dimensions['N_PROF'].size):
        date = ''.join(nc_token['mon_day_yr'][date_header].tolist())
        Lon = nc_token['Lon'][date_header]
        if np.abs(Lon)>360:
            print 'Longitude is showing some funny values for float ',cruise
        Lat = nc_token['Lat'][date_header]
        if np.abs(Lat)>90:
            print 'Latitude is showing some funny values for float ',cruise
        PosQC = nc_token['Lat_QF'][date_header]
        dataframe_token = pd.DataFrame()
        for variable in variable_list: 
            dataframe_token[variable]=nc_token[variable][date_header]
            dataframe_token[variable+'_QF']=nc_token[variable+'_QF'][date_header]
        dataframe_token['Cruise']=cruise
        dataframe_token['Date']=pd.to_datetime(date)
        dataframe_token['Lat']=Lat
        dataframe_token['Lon']=Lon
        dataframe_token['PosQC']=PosQC
        frame_token.append(dataframe_token)
    holder = pd.concat(frame_token)
    return holder

def soccom_df(data_directory):
    frames = []             #initialize the variable that we will append dataframes to 
    file_names = glob.glob(data_directory)   # develop list of files in the soccom data directory
    for float_name in file_names:
        print float_name
        frames.append(soccom_file_reader(float_name))           #pandas is orders of magnitude faster at appending then concatenating
    df_holder = pd.concat(frames)
    df_holder = df_holder.dropna(subset=['Date','Lat','Lon'])
    df_holder.loc[df_holder.PosQC=='1','Lon']=np.nan
    df_holder.loc[df_holder.PosQC=='1','Lat']=np.nan
    df_holder.loc[df_holder.PosQC=='4','PosQC']=8
    df_holder.loc[df_holder.PosQC=='0','PosQC']=1
    for variable in variable_list:
        df_holder.loc[df_holder[variable+'_QF']!='0',variable]=np.nan
    df_holder['Type']='SOCCOM'
    df_holder = df_holder.dropna(subset=['Date','Lat','Lon'])
    df_holder = df_holder.sort_values(['Cruise','Date','Pressure'])    #sort in a reasonable manner
    df_holder = df_holder.reset_index(drop=True)
    df_holder = df_holder[['Date','Lat','Lon','Cruise','PosQC','Type']+variable_list]

    return df_holder

df_soccom = soccom_df(soccom_proj_settings.soccom_data_directory)
df_soccom['Cruise']=df_soccom['Cruise'].astype(str)
df_soccom.loc[df_soccom.Lon.values<0,['Lon']] = df_soccom[df_soccom.Lon<0].Lon.values+360
df_soccom.to_pickle(soccom_proj_settings.soccom_drifter_file)
