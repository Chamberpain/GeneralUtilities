from GeneralUtilities.Data.Filepath.instance import get_data_folder
from GeneralUtilities.Data.Download.download_utilities import execute_download
import os
data_folder = os.path.join(get_data_folder(),'Processed/glob_current/')


def execute_download():
    data_dict_list = [{'region':'hawaii.nc','lllon':-159,'urlon':-154,'lllat':16,'urlat':22},
    {'region':'puerto_rico.nc','lllon':-68.5,'urlon':-65,'lllat':16,'urlat':22.5},
    {'region':'monterey.nc','lllon':-126,'urlon':-121.5,'lllat':34,'urlat':39}]

    for data_dict in data_dict_list:
    	copernicusmarine.subset(
      dataset_id="cmems_obs-mob_glo_phy-cur_my_0.25deg_PT1H-i",
      variables=["ue", "uo", "utide", "ve", "vo", "vtide"],
      minimum_longitude=data_dict['lllon'],
      maximum_longitude=data_dict['urlon'],
      minimum_latitude=data_dict['lllat'],
      maximum_latitude=data_dict['urlat'],
      start_datetime="1993-01-01T23:00:00",
      end_datetime="2023-12-31T23:00:00",
      minimum_depth=0,
      maximum_depth=0,
      output_filename = data_dict['region'],
      output_directory = data_folder

    )
