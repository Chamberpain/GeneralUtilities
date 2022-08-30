from GeneralUtilities.Compute.list import TimeList
import datetime


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
	try:
		data_type = data_instance.datatype
	except AttributeError:
		data_type = data_instance.dtype
	masked_array = data_instance[:]
	if data_type in ['S1','S']:
		data = ''.join([_.decode("utf-8",errors='ignore') for _ in masked_array[~masked_array.mask].data.ravel()])
	elif data_type == 'float64':
		data = [_ for _ in masked_array[~masked_array.mask].data.ravel()]
		if len(data)==1:
			data = data[0]
	else: 
		print(data_type)
		print('I dont know what to do with this data')
		raise
	return data

def julian_time_parse(array,reference_date):
	time_instance = [reference_date +datetime.timedelta(days=_) for _ in array]
	return time_instance

def parse_time(variable_instance,reference_date=None):
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
