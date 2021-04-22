from __future__ import print_function
from data_save_utilities.lagrangian.ocean_parcels.base_oceanparcels_read import OPReader

class OneDJet(OPReader):
    file_path = '/Users/pchamberlain/Projects/data_save_utilities/lagrangian/ocean_parcels/OneDJet/'
    data_description = 'OP_OneDJet'
    def __init__(self,x,y,t,k):
        super(OneDJet,self).__init__(x,y,t,k)
