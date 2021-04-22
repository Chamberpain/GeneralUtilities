import numpy as np
from data_save_utilities.lagrangian.ocean_parcels.idealized_base import IdealizedBase

class OneDDiffusion(IdealizedBase):
    file_path = '/Users/pchamberlain/Projects/data_save_utilities/lagrangian/ocean_parcels/OneDDiffusion/'
    def __init__(self,K_bar = 0.1,alpha = 1):
        super(OneDDiffusion,self).__init__(K_bar = K_bar,alpha = alpha)

    def construct_diffusivity_and_velocity(self):
        beta = np.zeros(self.y_K.shape)         # Placeholder for fraction term in K(y) formula
        Kh_meridional =self.K_bar*np.ones(self.y_K.shape)
        Kh_meridional = np.concatenate((np.array([0]), Kh_meridional, np.array([0])))
        Kh_meridional[:3] = 0
        Kh_meridional[-3:] = 0
        _,Kh_meridional = np.meshgrid(Kh_meridional,Kh_meridional)
        U = self.speed*np.ones([self.Ny,self.Ny])
        return {'U': U,
                'V': np.zeros([self.Ny,self.Ny]),
                'Kh_zonal': self.K_bar*np.ones([self.Ny,self.Ny]),
                'Kh_meridional':  Kh_meridional}





