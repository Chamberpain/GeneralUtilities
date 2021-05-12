from __future__ import print_function

def find_nearest(items, pivot,test=True):
	nearest = min(items, key=lambda x: abs(x - pivot))
	item_range = max(items)-min(items)
	if test:
		assert (nearest-pivot)<0.1*item_range # only allow 10% extrapolation
	return	nearest

def flat_list(non_flat_list):
	flat_list = [item for sublist in non_flat_list for item in sublist]
	return flat_list

def time_parser(juld_list,ref_date = datetime.datetime(1950,1,1,1,1)):
	return [ref_date + datetime.timedelta(days=x) for x in juld_list]