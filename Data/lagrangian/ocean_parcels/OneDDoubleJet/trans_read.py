from __future__ import print_function
from data_save_utilities.lagrangian.ocean_parcels.base_oceanparcels_read import OPReader

class OneDDoubleJet(OPReader):
    file_path = '/Users/pchamberlain/Projects/data_save_utilities/lagrangian/ocean_parcels/OneDDoubleJet/'
    data_description = 'OP_OneDDoubleJet'
    def __init__(self,x,y,t,k):
        super(OneDDoubleJet,self).__init__(x,y,t,k)
