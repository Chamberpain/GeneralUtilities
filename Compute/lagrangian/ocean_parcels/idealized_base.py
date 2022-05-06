from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
from parcels import FieldSet, ParticleSet, JITParticle, plotTrajectoriesFile
from parcels import AdvectionDiffusionM1
from datetime import timedelta
import xarray as xr

class IdealizedBase(object):
	L = 1. # Basin scale
	Ny = 103 # Number of grid cells in y_direction (101 +2, one level above and one below, where fields are set to zero)
	y = np.linspace(-0.01, 1.01, Ny)
	y_K = np.linspace(0., 1., (Ny-2))     # y-coordinates used for setting diffusivity
	seconds = 1
	speed = 10.
	def __init__(self,K_bar = 0.25,alpha = 1.):
		self.K_bar = K_bar # Average diffusivity
		self.alpha = alpha # Profile steepness
		self.data = self.construct_diffusivity_and_velocity()

	def recompile(self):

		def periodicBC(particle, fieldset, time):
			if particle.lon < fieldset.halo_west:
				particle.lon += fieldset.halo_east - fieldset.halo_west
			elif particle.lon > fieldset.halo_east:
				particle.lon -= fieldset.halo_east - fieldset.halo_west

		fieldset = self.construct_fieldset()
		filename = self.file_path+"output.nc"

		dt = 0.0001
		pset = ParticleSet.from_list(fieldset,
									 pclass=JITParticle,
									 lon=np.random.uniform(size=1000,low=0.0,high=1.0),
									 lat=np.random.uniform(size=1000,low=0.1,high=0.9),
									 time=np.zeros(1000),
									 lonlatdepth_dtype=np.float64)
		output_file = pset.ParticleFile(name=filename,
										 outputdt=timedelta(seconds=dt))
		pset.execute(AdvectionDiffusionM1 + pset.Kernel(periodicBC),
					  runtime=timedelta(seconds=self.seconds),
					  dt=timedelta(seconds=dt),
					  output_file=output_file,
					  verbose_progress=True)
		output_file.close()

	def construct_fieldset(self):
		dims = {'lon': np.linspace(-0.01, 1.01, self.Ny, dtype=np.float32),
				'lat': np.linspace(-0.01, 1.01, self.Ny, dtype=np.float32)}
		fieldset = FieldSet.from_data(self.data, dims, mesh='flat', allow_time_extrapolation=True)
		fieldset.add_constant('dres', 0.00005)
		fieldset.add_constant('halo_west', fieldset.U.grid.lon[0])
		fieldset.add_constant('halo_east', fieldset.U.grid.lon[-1])
		fieldset.add_periodic_halo(zonal=True)
		return fieldset

	def plot_initial_conditions(self):
		y = np.linspace(-0.01, 1.01, 103)
		plt.subplot(1,2,1)

		plt.plot(self.data['Kh_meridional'][:,0],self.y)
		plt.ylabel("y")
		plt.xlabel(r"$K_{meridional}$")
		plt.subplot(1,2,2)
		plt.plot(self.data['U'][:,0],self.y)
		plt.ylabel("y")
		plt.xlabel(r"$Velocity$")
		plt.savefig(self.file_path+'initial_conditions')
		plt.close()

		
	def plot_particles(self):
		filename = self.file_path+"output.nc"
		M1_out = xr.open_dataset(filename)
		fig, ax = plt.subplots(1, 2)
		fig.set_figwidth(12)

		for data, ai, dim, ystart, ylim in zip([M1_out.lat, M1_out.lon], ax, ('y', 'x'), (0.75, 0), [(0, 1), (0, 1)]):
			ai.plot(np.arange(0, self.seconds+0.0002, 0.0001), data.T, alpha=0.2)
			ai.set_xlabel("t")
			ai.set_ylabel(dim)
			ai.set_xlim(0, self.seconds)
			ai.set_ylim(ylim)
		fig.suptitle("`AdvectionDiffusionM1` Simulation: Particle trajectories in the x- and y-directions against time")
		plt.savefig(self.file_path+'particle_plot_1')
		plt.close()
		fig, ax = plt.subplots(1, 2)
		fig.set_figwidth(12)
		mask = M1_out.lat.data[:,0]<0.5
		for data, ai, dim, ystart, ylim in zip([M1_out.lat[mask,:], M1_out.lon[mask,:]], ax, ('y', 'x'), (0.75, 0), [(0, 1), (0, 1)]):
			ai.plot(np.arange(0, (self.seconds+0.0002), 0.0001), data.T, alpha=0.2)
			ai.set_xlabel("t")
			ai.set_ylabel(dim)
			ai.set_xlim(0, self.seconds)
			ai.set_ylim(ylim)
		fig.suptitle("`AdvectionDiffusionM1` Simulation: Particle trajectories in the x- and y-directions against time")
		plt.savefig(self.file_path+'particle_plot_2')
		plt.close()
		plotTrajectoriesFile(filename, mode='hist2d', bins=[30, 20],show_plt=False)
		plt.savefig(self.file_path+'particle_hist')
		plt.close()