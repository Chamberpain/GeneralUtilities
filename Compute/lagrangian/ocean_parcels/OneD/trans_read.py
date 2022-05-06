from __future__ import print_function
from data_save_utilities.lagrangian.ocean_parcels.base_oceanparcels_read import OPReader

class OneD(OPReader):
    file_path = '/Users/pchamberlain/Projects/data_save_utilities/lagrangian/ocean_parcels/OneD/'
    data_description = 'OP_OneD'
    def __init__(self,x,y,t,k):
        super(OneD,self).__init__(x,y,t,k)
