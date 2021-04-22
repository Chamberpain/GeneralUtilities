from __future__ import print_function
import os
from netCDF4 import Dataset
import datetime
import geopy
import re
import numpy as np 
import geopy.distance
from data_save_utilities.lagrangian.drifter_base_class import BasePosition,Speed,BaseRead
from data_save_utilities.pickle_utilities import load,save



def data_adjust(data_instance, data_adjusted_instance):
	"""
	Function replaces adjusted data in the original data array

	Parameters
	----------
	data_instance: net cdf data instance
	adjusted: adjusted net cdf data instance

	Returns
	-------
	data processed in a usable form (string, array, float)
	"""
	masked_data_array = data_instance[:].data
	masked_adjusted_instance = data_adjusted_instance[:]
	masked_data_array[~masked_adjusted_instance.mask]=masked_adjusted_instance[~masked_adjusted_instance.mask]
	masked_data_array = np.ma.masked_equal(masked_data_array,data_instance._FillValue)
	return masked_data_array

def data_return(data_instance):
	"""
	Function reads in netcdf data instance and returns sensible data

	Parameters
	----------
	data_instance: net cdf data instance
	mask: specifies whether to mask the values or return fill value
	Returns
	-------
	data processed in a usable form (string, array, float)
	"""
	data_type = data_instance.datatype
	masked_array = data_instance[:]
	if data_type in ['S1','S']:
		data = ''.join([_.decode("utf-8") for _ in masked_array[~masked_array.mask].data.ravel()])
	elif data_type == 'float64':
		data = [_ for _ in masked_array[~masked_array.mask].data.ravel()]
		if len(data)==1:
			data = data[0]
	else: 
		print(data_type)
		print('I dont know what to do with this data')
		raise
	return data

def format_byte_list_to_string(byte_list):
	return ['None' if x is None else x.decode("utf-8") for x in byte_list]

class BaseDate(object):
	name = 'Date'

	def julian_time_parse(self,array,reference_date):
		time_instance = [reference_date +datetime.timedelta(days=_) for _ in array]
		return time_instance

	def parse_time(self,variable_instance,reference_date=None):
		"""
		Function parses a netcdf time related instance and returns a time instance

		Parameters
		----------
		variable instance of the time you wish decyfered
		reference date for juld 

		Returns
		-------
		date time instance
		"""

		if variable_instance.conventions == 'YYYYMMDDHHMISS':
			time_string = data_return(variable_instance)
			if time_string: #this tests for an empty string
				try:
					time_instance = datetime.datetime.strptime(time_string,'%Y%m%d%H%M%S')
				except ValueError:
					time_instance = None
			else:
				time_instance=None
		elif variable_instance.conventions == 'Relative julian days with decimal part (as parts of day)':
			if not reference_date:
				print('I need a reference date for Julian data')
				raise
			time_instance = julian_time_parse(data_return(variable_instance),reference_date)
		else:
			print(variable_instance.conventions)
			print('I dont know what to do with this')
			raise
		return time_instance

	def return_mask(self):
		return self._mask

	def is_problem(self):
		"""
		Returns boolean value to show that all the date tests have been passed.

		Current tests:
		All time differences are greater than 0
		All profiles happen before today 
		All profiles happen after the beginning of the argo program
		The maximum time difference between profiles cannot be greater than 9 months

		Returns
		-------
		Boolean value to show if this profile is a problem
		"""

		time_diff_list = [self._list[idx+1]-self._list[idx] for idx in range(len(self._list)-1)]
		seconds_diff_list = [_.days*24+_.seconds/3600. for _ in time_diff_list]
		test_1 = (np.array(seconds_diff_list)>0).all()	# make sure all time is going in the right direction
		if not test_1:
			print('time is not going in the right direction')
		test_2 = ((np.array(self._list))<datetime.datetime.today()).all() # make sure all profiles happen before today
		if not test_2:
			print('a profile happened in the future')
		test_3 = ((np.array(self._list))>datetime.datetime(1998,1,1)).all() # make sure all profiles are after beginning of argo program
		if not test_2:
			print('a profile happened before argo existed')
		test_4 = ((np.array([_.days for _ in time_diff_list])<270)).all() # make sure maximum time difference between profiles is less than 9 months
		if not test_2:
			print('maximum time between profiles is greater than 270 days')
		return ~(test_1&test_2&test_3&test_4)


class ArgoReader(BaseRead):
	"""
	class holds all argo meta, trajectory, profile data from argo meta class objects

	Parameters
	----------
	folder containing all data file opened with the netCDF4 Dataset library
	meta: boolean object that instructs file to open meta data netcdf file
	traj: boolean object that instructs file to open traj netcdf file and read
	prof: boolean object that instructs file to open profile netcdf file and read
	** note **
	these boolean objects are only to improve efficiency. By default all data will be loaded

	Returns
	-------
	Class object of netcdf file data.
	All lat/lon information formatted to -180 to 180
	All times are returned as datetime objects

	"""
	data_description = 'argo'
	def __init__(self,nc_folder,meta=True, traj= True, prof=True,tech=True):
		self.folder = nc_folder
		files = os.listdir(nc_folder)
		if meta:
			self.meta = self._file_reader(nc_folder,files,'.*meta.nc',self.MetaClass)
		if traj:
			self.traj = self._file_reader(nc_folder,files,'.*_traj.nc',self.TrajClass)
		if prof:
			self.prof = self._file_reader(nc_folder,files,'.*prof.nc',self.ProfClass)
		if tech:
			self.tech = self._file_reader(nc_folder,files,'.*tech.nc',self.TechClass)
		super(ArgoReader,self).__init__()

	def _file_reader(self,nc_folder,files,data_format,data_class):
		match_alg = re.compile(data_format)
		file_name = next(iter(list(filter(match_alg.match,files))+[None]))
		if file_name:
			print(os.path.join(nc_folder,file_name))
			return data_class(os.path.join(nc_folder,file_name))
		else:
			return None


	class BaseReadClass(object):
		pass
		""" base class that contains functions used by both ProfClass and TrajClass
			----------

		"""	
		class Position(BasePosition):
			""" class that is a holder of the position information

			"""     
			def __init__(self,nc_fid,mask):			
				self._list = [geopy.Point(_[0],_[1]) for _ in zip(nc_fid['LATITUDE'][mask],nc_fid['LONGITUDE'][mask])]
				if not self._list:
					print('the position information was empty, to prove it, take a look at')
					print(nc_fid['LATITUDE'][:])
					print(nc_fid['LONGITUDE'][:])

	class MetaClass():
		""" class to organize all of read in meta net cdf data
			----------
			file_path: the file path of the net cdf file you wish to read

		"""
		name = 'Meta'
		def __init__(self,file_path):
			print('I am opening Meta File')
			nc_fid = Dataset(file_path)
			self.id = data_return(nc_fid['PLATFORM_NUMBER'])
			self.date = self.Date(nc_fid)
			self.positioning_system = self._positioning_system_format(nc_fid.variables['POSITIONING_SYSTEM'])
			if bool(data_return(nc_fid.variables['LAUNCH_QC'])):
				launch_lat = data_return(nc_fid.variables['LAUNCH_LATITUDE'])
				launch_lon = data_return(nc_fid.variables['LAUNCH_LONGITUDE'])
				self.launch_loc = geopy.Point(launch_lat,launch_lon)
			nc_fid.close()


		def _positioning_system_format(self,data_instance):
			""" Performs necessary checks on data instance
				Parameters
				----------
				data_instance: net cdf positioning system data instance

				Returns
				-------
				positioning system string
			"""
			positioning_string = data_return(data_instance)
			if not positioning_string:
				return None
			if positioning_string in ['IRIDIUM','GTS','GPSIRIDIUM','IRIDIUMGPS','GPSIRIDIUMRAFOS']:
				positioning_string = 'GPS'
			assert positioning_string in ['ARGOS','GPS']
			return positioning_string

		class Date(BaseDate):
			def __init__(self,nc_fid):
				self.launch_date = self.parse_time(nc_fid.variables['LAUNCH_DATE'])
				if bool(data_return(nc_fid.variables['START_DATE_QC'])):
					self.start_date = self.parse_time(nc_fid.variables['START_DATE'])

	class TrajClass(BaseReadClass):
		""" class to organize all of read in trajectory net cdf data
			----------
			file_path: the file path of the net cdf file you wish to read

		"""	
		name = 'Traj'
		def __init__(self,file_path):
			print('I am opening Traj File')
			nc_fid = Dataset(file_path)
			self.date = self.Date(nc_fid)
			mask = self.date.return_mask()
			self.position_accuracy = nc_fid['POSITION_ACCURACY'][mask]
			self.pos = self.Position(nc_fid,mask)
			self.speed = Speed(self.date,self.pos,speed_limit=5)
			nc_fid.close()

		class Date(BaseDate):
			def __init__(self,nc_fid):
				try:
					date_data = data_adjust(nc_fid['JULD'],nc_fid['JULD_ADJUSTED'])
				except IndexError:
					date_data = nc_fid['JULD'][:]
				mask = (np.array([_ in ['1','2'] for _ in format_byte_list_to_string(nc_fid['POSITION_QC'][:].tolist())]))&\
				(~nc_fid['LONGITUDE'][:].mask)&\
				(~nc_fid['LATITUDE'][:].mask)&\
				(~nc_fid['POSITION_ACCURACY'][:].mask)&\
				(~date_data.mask)
				self._list = self.julian_time_parse(date_data[mask].data, self.parse_time(nc_fid['REFERENCE_DATE_TIME']))
				self._mask = mask


		# def assign_depth(self,depth):
		# 	try:
		# 		self.depth = [depth.return_z(_) for _ in self.pos]

	class ProfClass(BaseReadClass):
		""" class to organize all of read in profile net cdf data
			----------
			file_path: the file path of the net cdf file you wish to read

		"""		
		name = 'Profile'
		def __init__(self,file_path):
			print('I am opening Prof File')
			nc_fid = Dataset(file_path)
			self.date = self.Date(nc_fid)
			mask = self.date.return_mask()
			self.pos = self.Position(nc_fid,mask)
			self.speed = Speed(self.date,self.pos,speed_limit=5)
			nc_fid.close()

		class Date(BaseDate):
			def __init__(self,nc_fid):
				date_data = nc_fid['JULD_LOCATION'][:]
				mask = (np.array([_ in ['1','2'] for _ in format_byte_list_to_string(nc_fid['POSITION_QC'][:].tolist())]))&\
				(~nc_fid['LONGITUDE'][:].mask)&\
				(~nc_fid['LATITUDE'][:].mask)&\
				(~date_data.mask)
				
				self._list = self.julian_time_parse(date_data[mask].data, self.parse_time(nc_fid['REFERENCE_DATE_TIME']))
				self._mask = mask

	class TechClass():
		name = 'Tech'
		def __init__(self,file_path):
			print('I am opening Tech File')
			nc_fid = Dataset(file_path)
			self.drift_depth = self.DriftDepth(nc_fid)


		class DriftDepth():
			def __init__(self,nc_fid):


				variable_name_list = [''.join(format_byte_list_to_string(i)).replace(" ","") for i in nc_fid['TECHNICAL_PARAMETER_NAME'][:].data]
				
				indices = [] 
				variable_list = ['PRESSURE_InternalVacuumProfileStart_mbar','PRES_ParkMaximum_dBAR','PRES_ParkMinimum_dBAR']
				indices += [i for i, x in enumerate(variable_name_list) if x in variable_list]
				# if not indices:
				# 	print(np.unique(variable_name_list))

				holder = [''.join(format_byte_list_to_string(nc_fid['TECHNICAL_PARAMETER_VALUE'][:].data[_,:])).replace(" ","") for _ in indices]
				self._list = []
				for dummy in holder:
					self._list.append(int(''.join([_ for _ in dummy if _.isdigit()])))
			


			def is_problem(self):
				print(self._list)
				if not self._list:
					return False
				test_1 = (np.array(self._list)>500).all()	# allow a window of 500 on either side of drift depth
				if not test_1:
					print('drift depth less than 500 mb')
				test_2 = (np.array(self._list)<1500).all() # 
				if not test_2:
					print('drift depth greater than 1500 mb')
				return ~(test_1&test_2)

def aggregate_argo_list(num=-1,list_save=False,recompile=False):
	"""
	function that returns dictionary of argo read classes

	Parameters
	----------
	num: the length of the dictionary you want (default is all files)
	list_save: whether you want to save the list (default is false)

	Returns
	-------
	dictionary of argo read classes.
	"""
	__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))


	def compile_matches():
		base_folder = '/Users/pchamberlain/Data/Raw/ARGO/'
		matches = []
		for root, dirnames, filenames in os.walk(base_folder):
			meta_match = re.compile('.*meta.nc') # all folders will have a meta file
			if filter(meta_match.match,filenames):
				matches.append(root)
		save(__location__+"/argo_matches.pkl",matches)
		return matches

	# depth = Depth()
	if recompile:
		compile_matches()
	else:
		try:
			matches = load(__location__+"/argo_matches.pkl")
		except IOError:
			print('there was an error loading '+__location__+"argo_matches.pkl")
			matches = compile_matches()

	num = 0
	while num<len(matches):
		print('I am opening number match ',num)
		yield ArgoReader(matches[num])
		num +=1 


# def time_parser(date,ref_date):
#     try:
#         x = jdcal.jd2gcal(date,ref_date)
#         year = x[0]
#         month = x[1]
#         day = x[2]
#         return [datetime.date(year,month,day)]
#     except:
#         return [np.nan]


# def list_multiplier(list_,col_num):
#     list_ = [[x]*col_num for x in list_]
#     list_ = [item for sublist in list_ for item in sublist]
#     return list_


# def argo_file_reader(file_):
#     nc_fid = Dataset(file_, 'r')
#     pos_type = ''.join([x for x in nc_fid.variables['POSITIONING_SYSTEM'][:][0].tolist() if x is not None])
#     cruise = ''.join([x for x in nc_fid.variables['PLATFORM_NUMBER'][:][0].tolist() if x is not None])
#     mode = nc_fid.variables['DATA_MODE'][:]
#     frames = []
#     for variable in np.unique(mode).data:
#         truth_list = mode==variable
#         if variable is 'R':
#             lat = nc_fid.variables['LATITUDE'][truth_list].tolist()
#             lon = nc_fid.variables['LONGITUDE'][truth_list].tolist()
#             pos_qc = nc_fid.variables['POSITION_QC'][truth_list].tolist()
#             date = nc_fid.variables['JULD'][truth_list].tolist()
#             date_qc = nc_fid.variables['JULD_QC'][truth_list]
#             sal = nc_fid.variables['PSAL'][truth_list,:].flatten().tolist()
#             pressure = nc_fid.variables['PRES'][truth_list,:].flatten().tolist()
#             temp = nc_fid.variables['TEMP'][truth_list,:].flatten().tolist()

#             sal_error = nc_fid.variables['PSAL_ERROR'][truth_list,:].flatten().tolist()
#             pressure_error = nc_fid.variables['PRES_ERROR'][truth_list,:].flatten().tolist()
#             temp_error = nc_fid.variables['TEMP_ERROR'][truth_list,:].flatten().tolist()
#             cycle = nc_fid.variables['CYCLE_NUMBER'][truth_list].tolist()

#             sal_qc = nc_fid.variables['PSAL_QC'][truth_list,:].flatten().tolist()
#             pressure_qc = nc_fid.variables['PRES_QC'][truth_list,:].flatten().tolist()
#             temp_qc = nc_fid.variables['TEMP_QC'][truth_list,:].flatten().tolist()

#         else:
#             lat = nc_fid.variables['LATITUDE'][truth_list].tolist()
#             lon = nc_fid.variables['LONGITUDE'][truth_list].tolist()
#             pos_qc = nc_fid.variables['POSITION_QC'][truth_list].tolist() 
#             try:
#                 date = nc_fid.variables['JULD_ADJUSTED'][truth_list].tolist()
#                 date_qc = nc_fid.variables['JULD_ADJUSTED_QC'][truth_list]
#             except KeyError:
#                 date = nc_fid.variables['JULD'][truth_list].tolist()
#                 date_qc = nc_fid.variables['JULD_QC'][truth_list]
#             sal = nc_fid.variables['PSAL_ADJUSTED'][truth_list,:].flatten().tolist()
#             pressure = nc_fid.variables['PRES_ADJUSTED'][truth_list,:].flatten().tolist()
#             temp = nc_fid.variables['TEMP_ADJUSTED'][truth_list,:].flatten().tolist()

#             sal_error = nc_fid.variables['PSAL_ADJUSTED_ERROR'][truth_list,:].flatten().tolist()
#             pressure_error = nc_fid.variables['PRES_ADJUSTED_ERROR'][truth_list,:].flatten().tolist()
#             temp_error = nc_fid.variables['TEMP_ADJUSTED_ERROR'][truth_list,:].flatten().tolist()
#             try:
#                 cycle = nc_fid.variables['CYCLE_NUMBER_ADJUSTED'][truth_list].tolist()
#             except KeyError:
#                 cycle = nc_fid.variables['CYCLE_NUMBER'][truth_list].tolist()
#             sal_qc = nc_fid.variables['PSAL_ADJUSTED_QC'][truth_list,:].flatten().tolist()
#             pressure_qc = nc_fid.variables['PRES_ADJUSTED_QC'][truth_list,:].flatten().tolist()
#             temp_qc = nc_fid.variables['TEMP_ADJUSTED_QC'][truth_list,:].flatten().tolist()
#         dummy, col_num = nc_fid.variables['PSAL'][truth_list,:].shape

#         lat = list_multiplier(lat,col_num)
#         lon = list_multiplier(lon,col_num)
#         pos_qc =  list_multiplier(pos_qc,col_num)
#         cycle = list_multiplier(cycle,col_num)

#         ref_date = ''.join(nc_fid.variables['REFERENCE_DATE_TIME'][:].tolist())
#         ref_date = sum(jdcal.gcal2jd(ref_date[:4],ref_date[4:6],ref_date[6:8]))
#         date = list_multiplier([item for sublist in [time_parser(x,ref_date) for x in date] for item in sublist],col_num)
#         date_qc = list_multiplier(date_qc,col_num)

#         frames.append(pd.DataFrame({'Cycle':cycle,'Date':date,'DateQC':date_qc,'Lon':lon,'Lat':lat,
#             'Pressure':pressure,'Temperature':temp,'Salinity':sal,'PosQC':pos_qc,'PressureQC':pressure_qc,'Pressureerror':pressure_error,
#             'SalQC':sal_qc,'Salerror':sal_error,'TempQC':temp_qc,'Temperror':temp_error}))

#     df_holder = pd.concat(frames)

#     df_holder['Type']=pos_type
#     df_holder['Cruise']=cruise

#     df_holder = df_holder.dropna(subset = ['Date'])
#     df_holder = df_holder[(df_holder.PressureQC.isin(['1','2']))&(df_holder.SalQC.isin(['1','2']))&(df_holder.TempQC.isin(['1','2']))&
#     (df_holder.PosQC.isin(['1','2','8']))&(df_holder.DateQC.isin(['1','2']))]

#     if (df_holder.Salerror.values>0.5).any()|(df_holder.Pressureerror.values>20).any()|(df_holder.Temperror.values>0.5).any():
#         # floats that are outside these parameters need to be eliminated because it is unclear what future data processing problems they may cause.
#         print 'I found a float outside of error parameters'
#         f.write(str(cruise[0])+' is rejected because it was found outside error parameters \n')
#         nc_fid.close()
#         return pd.DataFrame()
#     nc_fid.close()
#     return df_holder

# def argo_df(data_directory):
#     global f
#     f = open('argo_df_changelog.txt','w')
#     frames = []
#     matches = []
#     for root, dirnames, filenames in os.walk(data_directory): #this should be a class generator for paralization (accoording to gui)
#         for filename in fnmatch.filter(filenames, '*prof.nc'): #yield function should be used
#             matches.append(os.path.join(root, filename))
#     for n, match in enumerate(matches):
#         print 'file is ',match,', there are ',len(matches[:])-n,'floats left'
#         t = time.time()
#         frames.append(argo_file_reader(match))
#         if debug:
#             print 'Building and merging datasets took ', time.time()-t 
#     df_holder = pd.concat(frames)
#     df_holder.Date = pd.to_datetime(df_holder.Date)
#     df_holder = df_holder.dropna(subset = ['Pressure'])
#     f.close()
#     return df_holder

# def traj_file_reader(file_):
#     nc_fid = Dataset(file_, 'r')
#     lat = nc_fid.variables['LATITUDE'][:].tolist()
#     lon = nc_fid.variables['LONGITUDE'][:].tolist()
#     num = nc_fid.variables['CYCLE_NUMBER'][:].tolist()
#     pos_type = ''.join([x for x in nc_fid.variables['POSITIONING_SYSTEM'][:].tolist() if x is not None])
#     pos_acc = nc_fid.variables['POSITION_ACCURACY'][:]
#     pos_qc = nc_fid['POSITION_QC'][:]
#     pos_acc = np.ma.array(pos_acc,mask=(pos_qc!='1'))

#     date = nc_fid.variables['JULD'][:]
#     date_qc = nc_fid.variables['JULD_QC'][:]
#     date = np.ma.array(date,mask=(date_qc!='1'))

#     cruise = ''.join([a for a in nc_fid['PLATFORM_NUMBER'][:].tolist() if a])
#     ref_date = ''.join(nc_fid.variables['REFERENCE_DATE_TIME'][:].tolist())
#     ref_date = datetime.datetime.strptime(ref_date,'%Y%m%d000000').date()
#     date_list = []
#     for n,day in enumerate(date):
#         try:
#             date_list.append(ref_date+datetime.timedelta(days=day))
#         except TypeError:
#             date_list.append(np.nan)
#     # date = [item for sublist in [time_parser(x,ref_date) for x in date] for item in sublist]
#     df_holder = pd.DataFrame({'Cruise':cruise,'Date':date_list,'Lon':lon,'Lat':lat,'Position Type':pos_type,'Position Accuracy':pos_acc,'Position QC':pos_qc})
	

#     df_holder = df_holder.dropna(subset=['Date'])
#     df_holder = df_holder[((df_holder['Position Type']=='ARGOS')&(df_holder['Position Accuracy'].isin(['1','2','3'])))|(df_holder['Position Type']=='GPS')]
#     nc_fid.close()

#     if df_holder.empty:
#         return df_holder
#     df_holder = df_holder.dropna(subset=['Lat','Lon'])
#     df_holder = df_holder.drop_duplicates(subset=['Date'])
#     df_holder.Date = pd.to_datetime(df_holder.Date)

#     if  df_holder.Date.min()<pd.to_datetime(datetime.date(1996,1,1)):
#         print 'huge problem'
#         print df_holder.Date.min()
#         plt.figure()
#         plt.subplot(2,1,1)
#         df_holder.set_index('Date').Lat.plot()
#         plt.title('Latitude')
#         plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom='off',      # ticks along the bottom edge are off
#     top='off',         # ticks along the top edge are off
#     labelbottom='off') # labels along the bottom edge are off
#         plt.xlabel('')
#         plt.subplot(2,1,2)
#         df_holder.set_index('Date').Lon.plot()
#         plt.title('Longitude')
#         plt.savefig(cruise)
#         plt.close()
#         return pd.DataFrame()
#     if (np.diff(df_holder.Lat)/df_holder.Date.diff().dt.days.values[1:]>0.7).any()|((np.diff(oceans.wrap_lon360(df_holder.Lon))/df_holder.Date.diff().dt.days.values[1:]>0.7)&(np.diff(df_holder.Lon)/df_holder.Date.diff().dt.days.values[1:]>0.7)).any():
#         'the damn float is going too fast'
#         print df_holder[(df_holder['Lat'].diff().abs()/df_holder.Date.diff().dt.days>.7)|((df_holder['Lon'].diff().abs()/df_holder.Date.diff().dt.days>.7)&(df_holder['Lon'].apply(oceans.wrap_lon180).diff().abs()/df_holder.Date.diff().dt.days>.7))]
#         print df_holder[df_holder['Lon'].apply(oceans.wrap_lon180).diff().abs()/df_holder.Date.diff().dt.days>.7]
#         print df_holder[((df_holder['Lon'].diff().abs()/df_holder.Date.diff().dt.days>.7)&(df_holder['Lon'].apply(oceans.wrap_lon180).diff().abs()/df_holder.Date.diff().dt.days>.7))]
#         # if not df_holder[(df_holder['Lat'].diff().abs()/df_holder.Date.diff().dt.days>.3)].empty:
#         #     df_holder['Lat'] = pd.rolling_median(df_holder['Lat'], window=2, center=True).fillna(method='bfill').fillna(method='ffill')
#         #     print 'I smoothed the lat'
#         # if not df_holder[((df_holder['Lon'].diff().abs()/df_holder.Date.diff().dt.days>.3)&(df_holder['Lon'].apply(oceans.wrap_lon180).diff().abs()/df_holder.Date.diff().dt.days>.3))].empty:
#         #     df_holder['Lon'] = pd.rolling_median(df_holder['Lon'], window=2, center=True).fillna(method='bfill').fillna(method='ffill')
#         #     print 'I smoothed the lon'
#         plt.figure()
#         plt.subplot(2,1,1)
#         df_holder.set_index('Date').Lat.plot()
#         plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom='off',      # ticks along the bottom edge are off
#     top='off',         # ticks along the top edge are off
#     labelbottom='off') # labels along the bottom edge are off
#         plt.xlabel('')
#         plt.title('Latitude')
#         plt.subplot(2,1,2)
#         df_holder.set_index('Date').Lon.plot()
#         plt.title('Longitude')
#         plt.savefig(cruise)
#         plt.close()
#         return pd.DataFrame()
#     return df_holder

# def traj_df(data_directory):
#     frames = []
#     matches = []
#     float_type = ['Argo']
#     for root, dirnames, filenames in os.walk(root_dir):
#         meta_match = re.compile('.*meta.nc')
#         if filter(meta_match.match,filenames):
#             matches.append(root)
#     for n, match in enumerate(matches):
#         try:
#             print 'file is ',match,', there are ',len(matches[:])-n,'floats left'
#             t = time.time()
#             frames.append(traj_file_reader(match))
#             print 'Building and merging datasets took ', time.time()-t 
#         except IndexError:
#             continue
#     df_holder = pd.concat(frames)
#     df_holder.Date = pd.to_datetime(df_holder.Date)
#     return df_holder