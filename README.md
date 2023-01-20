# GeneralUtilities

It is necessary to create a symbolic link from the GeneralUtilities base folder to the folder that will contain your data. This symbolic linked folder will be called DataDir.

Data used in these projects can be automatically downloaded from the GeneralUtilities/Data/Download folder.
master_download.py will download run all the download scripts and store them in the DataDir.

If you use conda to manage your packages, it is recommended to create a new environment when using this code. Details on how to do that can be found here (https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands). Then the following will need to be installed:


conda install ipython
conda install -c conda-forge cartopy
conda install requests
conda install scipy
conda install -c conda-forge geopy
conda install netCDF4
conda install -c conda-forge gsw
