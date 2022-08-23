# GeneralUtilities

It is necessary to create a symbolic link from the GeneralUtilities base folder to the folder that will contain your data. This symbolic linked folder will be called DataDir.

The data download scripts can now automatically download the necessary locations. For the argo code, you will need to add either data or a symbolic link to DataDir/Raw/Argo.
Then in 
GeneralUtilities/Data/lagrangian/argo/argo_read.py
Run the function full_argo_list()
And you will compile a database of the argo files. These databases will be called all_dict_argo and all_dict_argo_bgc and will be in the DataDir/Raw/Argo folder.



