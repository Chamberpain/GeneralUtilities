from mpl_toolkits.basemap import pyproj
import numpy as np
from matplotlib.patches import Polygon
from mpl_toolkits.basemap import Basemap
from data_save_utilities.depth.depth_utilities import ETopo1Depth as Depth
import matplotlib.pyplot as plt


class Basemap(Basemap):

	@classmethod
	def auto_map(cls,urlat,lllat,urlon,lllon,lon_0,aspect=False,spacing=10,bar=False,depth=False,resolution='l'):

		m = cls(projection='cea',llcrnrlat=lllat,urcrnrlat=urlat,\
				llcrnrlon=lllon,urcrnrlon=urlon,resolution=resolution,lon_0=lon_0,\
				fix_aspect=aspect)
			# m.drawmapboundary(fill_color='darkgray')
		m.fillcontinents(color='dimgray',zorder=10)
		m.drawmeridians(np.arange(0,360,2*spacing),labels=[True,False,False,False])
		m.drawparallels(np.arange(-90,90,spacing),labels=[False,False,False,True])
		if depth:
			depth = Depth()
			depth.z[depth.z>0] = 0 
			depth.z[depth.z<-6000] = -6000 
			XX,YY = np.meshgrid(depth.x[::5],depth.y[::5])
			XX,YY = m(XX,YY)
			m.contourf(XX,YY,depth.z[::5,::5]/1000.,vmin=-6,vmax=0)
		if bar:
			m.colorbar(label='Depth (km)')

		return m

	def ellipse(self, x0, y0, a, b, n, phi=0,ax=None,line=True, **kwargs):
		ax = kwargs.pop('ax', None) or self._check_ax()
		g = pyproj.Geod(a=self.rmajor, b=self.rminor)
		# Gets forward and back azimuths, plus distances between initial
		# points (x0, y0)
		azf, azb, dist = g.inv([x0, x0], [y0, y0], [x0+a, x0], [y0, y0+b])
		tsid = dist[0] * dist[1] # a * b
		seg = []
		AZ = np.linspace(azf[0], 360. + azf[0], n)
		for i, az in enumerate(AZ):
			if np.allclose(0., y0) and (np.allclose(90., az) or
				np.allclose(270., az)):
				continue
			azr = 2. * np.pi / 360. * (az+phi+90)
			A = dist[0] * np.sin(azr)
			B = dist[1] * np.cos(azr)
			r = tsid / (B**2. + A**2.)**0.5
			lon, lat, azb = g.fwd(x0, y0, az, r)
			x, y = self(lon, lat)
			# Add segment if it is in the map projection region.
			if x < 1e20 and y < 1e20:
				seg.append((x, y))
		if line:
			seg = [self(x,y,inverse=True) for x,y in seg]
			return seg
		else: 
			poly = Polygon(seg, **kwargs)
			ax.add_patch(poly)
		# Set axes limits to fit map region.
			self.set_axes_limits(ax=ax)
			return poly
