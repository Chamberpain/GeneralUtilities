from GeneralUtilities.Data.Download.agva_download import download as agva_download
from GeneralUtilities.Data.Download.cm4_download import download as cm4_download
from GeneralUtilities.Data.Download.etopo_1_download import download as etopo_1_download
from GeneralUtilities.Data.Download.landschutzer_download import download as landschutzer_download
from GeneralUtilities.Data.Download.modis_download import download as modis_download
from GeneralUtilities.Data.Download.pacioos_bathy_download import download as pacioos_bathy_download
from GeneralUtilities.Data.Download.roemmich_gilson_download import download as roemmich_gilson_download


for d_func in [agva_download,cm4_download,etopo_1_download,landschutzer_download,modis_download,pacioos_bathy_download,roemmich_gilson_download]:
	d_func()
