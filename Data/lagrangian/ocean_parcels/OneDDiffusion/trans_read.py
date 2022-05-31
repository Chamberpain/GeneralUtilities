from __future__ import print_function
from data_save_utilities.lagrangian.ocean_parcels.base_oceanparcels_read import OPReader

class OneDDiffusion(OPReader):
    file_path = '/Users/pchamberlain/Projects/data_save_utilities/lagrangian/ocean_parcels/OneDDiffusion/'
    data_description = 'OP_OneDDiffusion'
    def __init__(self,x,y,t,k):
        super(OneDDiffusion,self).__init__(x,y,t,k)
