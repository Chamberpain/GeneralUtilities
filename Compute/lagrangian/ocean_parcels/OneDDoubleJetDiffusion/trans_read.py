from __future__ import print_function
from data_save_utilities.lagrangian.ocean_parcels.base_oceanparcels_read import OPReader

class OneDDoubleJetDiffusion(OPReader):
    file_path = '/Users/pchamberlain/Projects/data_save_utilities/lagrangian/ocean_parcels/OneDDoubleJetDiffusion/'
    data_description = 'OP_OneDDoubleJetDiffusion'
    def __init__(self,x,y,t,k):
        super(OneDDiffusion,self).__init__(x,y,t,k)
