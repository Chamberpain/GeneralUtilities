import os 
from GeneralUtilities.Data.Filepath.instance import get_data_folder
import subprocess 

argo_folder = os.path.join(get_data_folder(),'Raw/Argo/')

def download():
	subprocess.call(["rsync", "-avzh", "--delete","vdmzrs.ifremer.fr::argo/",argo_folder])