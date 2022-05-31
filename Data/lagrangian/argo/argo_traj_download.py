import ftplib
import pickle
import re
import os
import ftplib

"File to download all traj files from the ifremer server. This is a silly way to do things as they supposedly offer Rsync, but I only found this out after I coded this up."

root_folder = '/Users/paulchamberlain/Data/Traj/'
reg_express = re.compile('Rtraj.nc')
ftp_link = 'ftp.ifremer.fr'
ftp = ftplib.FTP(ftp_link)
ftp.login()
with open('backup.txt','rb') as fp:
	filenames = pickle.load(fp)
k = len(filenames)
for n,file_name in enumerate(filenames):
	print 'working on ',file_name
	print 'downloaded ',n,' files. ',(k-n),' remaining'
	save_name = root_folder+file_name.split('/')[-1] # use path as filename prefix, with underscores
	if file_name.split('/')[-1] in os.listdir(root_folder):
		print 'I am skipping because I found a duplicate'
		continue
	else:
		try:
			ftp.retrbinary(('RETR '+file_name), open(save_name, 'wb').write)
		except:
			raise