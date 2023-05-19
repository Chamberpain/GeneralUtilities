import os 
from GeneralUtilities.Data.Filepath.instance import get_data_folder
import subprocess 

argo_link = "https://www.seanoe.org/data/00311/42182/data/102270.tar.gz"
argo_folder = os.path.join(get_data_folder(),'Raw/Argo/')

def download():
	subprocess.call(["wget", argo_link, "-P" ,argo_folder])
	tarred_file = os.path.join(argo_folder, '102270.tar.gz')
	subprocess.call(["tar", "-xvzf", tarred_file])
	
