import numpy as np
from data_save_utilities.lagrangian.ocean_parcels.idealized_base import IdealizedBase

class OneDDoubleJetDiffusion(IdealizedBase):
    file_path = '/Users/pchamberlain/Projects/data_save_utilities/lagrangian/ocean_parcels/OneDDoubleJetDiffusion/'
    def __init__(self,K_bar = 0.1,alpha = 1,additive_diffusion=10**-2):
        self.additive_diffusion = additive_diffusion
        super(OneDDoubleJetDiffusion,self).__init__(K_bar = K_bar,alpha = alpha)
    def construct_diffusivity_and_velocity(self):
        beta = np.zeros(self.y_K.shape)         # Placeholder for fraction term in K(y) formula

        for yi in range(len(self.y_K)):
            if self.y_K[yi] < self.L/2:
                beta[yi] = self.y_K[yi]*np.power(self.L - 2*self.y_K[yi], 1/self.alpha)
            elif self.y_K[yi] >= self.L/2:
                beta[yi] = (self.L - self.y_K[yi])*np.power(2*self.y_K[yi] - self.L, 1/self.alpha)

        Kh_constant = 0.1*(2*(1+self.alpha)*(1+2*self.alpha))/(self.alpha**2*np.power(self.L, 1+1/self.alpha))
        Kh_meridional = Kh_constant*beta+self.additive_diffusion
        Kh_meridional = np.concatenate((np.array([0]), Kh_meridional, np.array([0])))
        Kh_meridional[:3] = 0
        Kh_meridional[-3:] = 0
        _,Kh_meridional = np.meshgrid(Kh_meridional,Kh_meridional)
        U = self.speed/Kh_constant*Kh_meridional
        return {'U': U,
                'V': np.zeros([self.Ny,self.Ny]),
                'Kh_zonal': self.K_bar*np.ones([self.Ny,self.Ny]),
                'Kh_meridional':  Kh_meridional}