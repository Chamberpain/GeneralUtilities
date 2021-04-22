from data_save_utilities.lagrangian.ocean_parcels.OneD.generator import OneD
from data_save_utilities.lagrangian.ocean_parcels.OneDDiffusion.generator import OneDDiffusion
from data_save_utilities.lagrangian.ocean_parcels.OneDJet.generator import OneDJet
from data_save_utilities.lagrangian.ocean_parcels.OneDDoubleJet.generator import OneDDoubleJet
from data_save_utilities.lagrangian.ocean_parcels.OneDDoubleJetDiffusion.generator import OneDDoubleJetDiffusion

recompile = True
for generator in [OneD,OneDDiffusion,OneDJet,OneDDoubleJet,OneDDoubleJetDiffusion]:
	dummy = generator()
	if recompile:
		dummy.recompile()
	dummy.plot_particles()
	dummy.plot_initial_conditions()