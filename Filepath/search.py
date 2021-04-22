def find_files(dir_,ext_):
	import os
	import fnmatch
	matches = []
	for root, dirnames, filenames in os.walk(dir_):
	    for filename in fnmatch.filter(filenames,ext_):
	        matches.append(os.path.join(root, filename))
	return matches