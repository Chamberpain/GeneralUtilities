import os
from GeneralUtilities.Data.Filepath.instance import get_data_folder
import subprocess

argo_link = "https://www.seanoe.org/data/00311/42182/data/102270.tar.gz"
argo_folder = os.path.join(get_data_folder(),'Raw/Argo/')

def download():
	# download tar file into argo folder and untar
	subprocess.call(["wget", argo_link, "-P" ,argo_folder])
	tarred_file = os.path.join(argo_folder, '102270.tar.gz')
	subprocess.call(["tar", "-xvzf", tarred_file, '-C', argo_folder])
	
	# untar nested tars
	dac_subfolder = os.path.join(argo_folder, "202305-ArgoData/dac")
	for f in os.listdir(dac_subfolder):
		if '.tar.gz' not in str(f):
			continue
		complete_path = os.path.join(dac_subfolder, f)
		subprocess.call(['tar', '-xvzf', complete_path, '-C', argo_folder])

	# delete temp folder and tar files
	to_del = os.path.join(argo_folder, '202305-ArgoData')
	subprocess.call('rm', '-rf', to_del)
	subprocess.call(["rm", tarred_file])

download()
