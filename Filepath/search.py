import os
import fnmatch
import time

def find_files(data_directory,fmt,function=str):
	frames = []
	matches = []
	for root, dirnames, filenames in os.walk(data_directory):
		for filename in fnmatch.filter(filenames,fmt):
			matches.append(os.path.join(root, filename))
	for n, match in enumerate(matches):
		print('file is ',match,', there are ',len(matches[:])-n,'floats left')
		t = time.time()
		frames.append(function(match))
		print('Building and merging datasets took ', time.time()-t)
	return frames