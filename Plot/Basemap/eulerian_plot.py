from __future__ import print_function
import os
from mpl_toolkits.basemap import Basemap,shiftgrid
import numpy as np
import pandas as pd
import pyproj
from matplotlib.patches import Polygon
import datetime
import matplotlib.pyplot as plt
import scipy.io
import pickle
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
from scipy.ndimage.filters import gaussian_filter
from data_save_utilities.depth.depth_utilities import ETopo1Depth as Depth

pacific_coord =[160,290,-75,-30] 
atlantic_coord=[-60,20,-75,-30]
big_weddell_coord=[-80,40,-77,-50]
weddell_coord=[-60,20,-77,-60]
ross_coord=[170,250,-77,-61]
lynne_coord =[-70,50,-80,-32]

map_features=[['Weddell Sea',-45,-70],['Ross Sea',195,-70], ['Ross Polynya',180, -80],['Maud Rise Polynya',0,-70],['Cosmonaut Polynya',45,-70],['Darnley Polynya',75,-70],['Weddel Polyna',-60,-70], ['Mertz Polyna',150,-70]]
datadir = os.getenv("HOME")+'/Data/Processed/eulerian_plot/basemap/'
sound_speed = 1.467 #km/s
slow = 1/sound_speed


#!!!!!!!!!!!!!   These variables should be accessed from a data file !!!!!!!!!!!!!!!
SOSO_weddel_coord = {'W1a':(-63.9538,-00.0035),'W1b':(-63.953,-00.006),'W1c':(-63.953,-00.028),'W1d':(-63.967,-00.052)
,'W1e':(-63.994,-00.044),'W1f':(-63.994,-00.044),'W2a':(-66.509,-00.0338),'W2b':(-66.511,-00.032)
,'W2c':(-66.511,-00.032),'W2d':(-66.511,-00.032),'W2e':(-66.511,-00.025),'W2f':(-66.515,-00.011)
,'W3a':(-64.4915,9.8255),'W3b':(-64.4915,9.8255),'W4a':(-66.618,-27.105),'W4b':(-66.615,-27.118)
,'W4c':(-66.612,-27.122),'W4d':(-66.608,-27.121),'W5a':(-65.619,-36.392),'W5b':(-65.617,-36.421)
,'W5c':(-65.621,-36.422),'W6a':(-63.703,-50.870),'W6b':(-63.712,-50.842),'W6c':(-63.718,-50.832)
,'W6d':(-63.702,-50.827),'W7a':(-69.0112,-34.0025),'W8a':(-65.5752,-37.1221),'W8b':(-65.5752,-37.1221)
,'W9a':(-69.061,-17.4315),'W9b':(-69.059,-17.384),'W9c':(-69.058,-17.389),'W10a':(-64.389,-45.873)
,'W10b':(-64.398,-45.866),'W10c':(-64.382,-45.869),'W11a':(-68.995,-6.945),'W11b1':(-69.005,-6.982)
,'W11b2':(-69.005,-6.982),'W11c':(-69.006,-6.983),'W12a':(-65.968,-12.252),'W13a':(-70.893,-28.891)
,'W14a':(-68.483,-44.111)}
SOSO_drake_coord = {'US1-1':(-59.914,-109.9105),'US1-2':(-55.996,-109.9102),'US1-3':(-58.1010,-97.96)
,'US1-4':(-59.9655,-86.1536),'US1-5':(-55.9630,-85.9863),'US1-6':(-58.3343,-74.4142),'US2-1':(-58.1690,-102.169)
,'US2-2':(-59.939,-78.003),'UK1-1':(-59.977,-65.990),'UK1-2':(-58.014,-61.886),'UK1-3':(-56.023,-57.784)
,'UK1-4':(-58.051,-53.645)}


def nabla(array):
	A,B = np.gradient(array)
	return A/(111/6.)+B/(111*np.cos(np.deg2rad(60))/6.)



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
	# def __init__(self,region=None,padding=2,dataframe=None,date=None,resolution='l',timedelta=0,cruise=None,soso_df=None,fix_aspect=True):
	# 	if region is None:
	# 		print 'Region must be specified, program is now exiting.'
	# 		raise NameError('No Region')
	# 	if (region is 'float')&(dataframe is None):
	# 		print 'Map cannot be made for float region if no float data is provided'
	# 		raise NameError('No Float Data')
	# 	if dataframe is not None:
	# 		df = dataframe
	# 	else:
	# 		df = pd.read_csv(os.path.join(datadir,'interpoalted_argo_all.csv'),index_col = 'Unnamed: 0',usecols=['Unnamed: 0','Cruise','Date','Lat','Lon','Type','PosQC'])
	# 	df['Lon180'] = df.Lon[:]
	# 	df.loc[df.Lon180.values>180,['Lon180']] = df[df.Lon180>180].Lon180-360
	# 	df = df.dropna(subset = ['Lat','Lon'])
	# 	df.Date = pd.to_datetime(df.Date)

	# 	if resolution is 'l':
	# 		self.subsample=4
	# 	elif resolution is 'm':
	# 		self.subsample=2
	# 	elif resolution is 'h':
	# 		self.subsample=1

	# 	# if date is not None:
	# 	self.date = date
	# 	#     df = df[(df.Date>=min(self.date,self.date+datetime.timedelta(days=timedelta)))&(df.Date<=max(self.date,self.date+datetime.timedelta(days=timedelta)))]
	# 	# else:
	# 	#     self.date=date
	# 	#     self.soso_df=soso_df
	# 	if cruise is not None:
	# 		self.cruise=cruise
	# 		df = df[df.Cruise==cruise]
	# 	self.dataframe = df



	# 	if region=='float':
	# 		lllat=self.dataframe.Lat.min()-padding
	# 		if lllat<-90:
	# 			lllat = -90
	# 		urlat=self.dataframe.Lat.max()+padding
	# 		if urlat>90:
	# 			urlat = 90
	# 		if abs(self.dataframe['Lon'].max()-self.dataframe['Lon'].min())<abs(self.dataframe['Lon180'].max()-self.dataframe['Lon180'].min()):
	# 			print 'going to 360 coord'
	# 			lllon=self.dataframe.Lon.min()-padding
	# 			urlon=self.dataframe.Lon.max()+padding
	# 			if urlon > 360: 
	# 				urlon = 360
	# 			if lllon<0:
	# 				lllon=0
	# 			self.lon_0=180

	# 		else: 
	# 			print 'going to 180 coord'
	# 			lllon=self.dataframe.Lon180.min()-padding
	# 			urlon=self.dataframe.Lon180.max()+padding
	# 			if urlon > 180: 
	# 				urlon = 180
	# 			if lllon<-180:
	# 				lllon=-180
	# 			self.lon_0=0
	# 			self.dataframe.loc[:,['Lon']]=self.dataframe.Lon180

	# 		Basemap.__init__(self,projection='mill',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution='c',lon_0=self.lon_0,fix_aspect=fix_aspect)
	# 	if region == 'polar':
	# 		Basemap.__init__(self,projection='spstere',boundinglat=-40,lon_0=180,resolution='f')
	# 		self.lon_0=0
	# 	elif region == 'atlantic':
	# 		lllat=atlantic_coord[2]
	# 		urlat=atlantic_coord[3]
	# 		lllon=atlantic_coord[0]
	# 		urlon=atlantic_coord[1]
	# 		self.lon_0=0
	# 		self.dataframe=self.dataframe[(self.dataframe.Lat>lllat)&(self.dataframe.Lat<urlat)&(self.dataframe.Lon180>lllon)&(self.dataframe.Lon180<urlon)]
	# 		Basemap.__init__(self,projection='mill',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution='c',lon_0=self.lon_0)
	# 		self.dataframe.loc[:,['Lon']]=self.dataframe.Lon180

	# 	elif region == 'pacific':
	# 		lllat=pacific_coord[2]
	# 		urlat=pacific_coord[3]
	# 		lllon=pacific_coord[0]
	# 		urlon=pacific_coord[1]
	# 		self.dataframe=self.dataframe[(self.dataframe.Lat>lllat)&(self.dataframe.Lat<urlat)&(self.dataframe.Lon>lllon)&(self.dataframe.Lon<urlon)]
	# 		self.lon_0=180
	# 		Basemap.__init__(self,projection='mill',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution='c',lon_0=self.lon_0)    

	# 	elif region == 'weddell':
	# 		lllat=weddell_coord[2]
	# 		urlat=weddell_coord[3]
	# 		lllon=weddell_coord[0]
	# 		urlon=weddell_coord[1]
	# 		self.lon_0=0
	# 		self.dataframe=self.dataframe[(self.dataframe.Lat>lllat)&(self.dataframe.Lat<urlat)&(self.dataframe.Lon180>lllon)&(self.dataframe.Lon180<urlon)]
	# 		Basemap.__init__(self,projection='mill',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution='c',lon_0=self.lon_0,fix_aspect=fix_aspect)
	# 		self.dataframe.loc[:,['Lon']]=self.dataframe.Lon180        

	# 	elif region == 'big_weddell':
	# 		lllat=big_weddell_coord[2]
	# 		urlat=big_weddell_coord[3]
	# 		lllon=big_weddell_coord[0]
	# 		urlon=big_weddell_coord[1]
	# 		self.lon_0=0
	# 		self.dataframe=self.dataframe[(self.dataframe.Lat>lllat)&(self.dataframe.Lat<urlat)&(self.dataframe.Lon180>lllon)&(self.dataframe.Lon180<urlon)]
	# 		Basemap.__init__(self,projection='mill',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution='c',lon_0=self.lon_0)
	# 		self.dataframe.loc[:,['Lon']]=self.dataframe.Lon180        


	# 	elif region == 'ross':
	# 		lllat=ross_coord[2]
	# 		urlat=ross_coord[3]
	# 		lllon=ross_coord[0]
	# 		urlon=ross_coord[1]
	# 		self.lon_0=180
	# 		self.dataframe=self.dataframe[(self.dataframe.Lat>lllat)&(self.dataframe.Lat<urlat)&(self.dataframe.Lon>lllon)&(self.dataframe.Lon<urlon)]
	# 		Basemap.__init__(self,projection='mill',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution='c',lon_0=self.lon_0)
	# 		self.dataframe.loc[:,['Lon']]=self.dataframe

	# 	elif region == 'lynne':
	# 		lllat=lynne_coord[2]
	# 		urlat=lynne_coord[3]
	# 		lllon=lynne_coord[0]
	# 		urlon=lynne_coord[1]
	# 		self.lon_0=0
	# 		self.dataframe=self.dataframe[(self.dataframe.Lat>lllat)&(self.dataframe.Lat<urlat)&(self.dataframe.Lon180>lllon)&(self.dataframe.Lon180<urlon)]
	# 		Basemap.__init__(self,projection='mill',llcrnrlat=lllat,urcrnrlat=urlat,llcrnrlon=lllon,urcrnrlon=urlon,resolution='f',lon_0=self.lon_0,fix_aspect=False)
	# 		self.dataframe.loc[:,['Lon']]=self.dataframe.Lon180        


	# 	self.drawcoastlines(linewidth=1.5)
	# 	self.fillcontinents(color=(0.57,0.57,0.57) ,zorder=8)


	def drifter_plot(self,markersz=2,lineplot=False,dataframe=None,drift_type=None,color='bo'):
		if drift_type is not None:
			df = self.dataframe[self.dataframe.PosQC==drift_type]
			label_name = drift_type            
		#     if drift_type == 'SOCCOM':
		#         color = 'ro'
		#     elif drift_type == 'Argo':
		#         color = 'yo'
		else:
			df = self.dataframe
			label_name = 'All Drifters'
		#     color = 'bo'
		if dataframe is not None:
			df = dataframe

		if lineplot is True:
			for name in df.Cruise.unique():
				mask = df.Cruise==name
				df.loc[mask,'Diff'] = df[mask]['Date'].diff().dt.days>30
				df.loc[mask,'Diff'] = df[mask].Diff.apply(lambda x: 1 if x else 0).cumsum()
				for g in df[mask].groupby('Diff').groups:
					frame = df[mask].groupby('Diff').get_group(g)
					xcoord = list(frame.Lon.values)
					ycoord = list(frame.Lat.values)
					xpos,ypos = self(xcoord,ycoord)
					self.plot(xpos,ypos,linewidth=0.5,color=color,zorder=1)               
		else:
			
			try:
				xcoord = list(df[df.PosQC=='Kalman Interpolation'].Lon.values)
				ycoord = list(df[df.PosQC=='Kalman Interpolation'].Lat.values)
				xpos,ypos = self(xcoord,ycoord)
				self.plot(xpos,ypos,'c*',markersize=markersz,label='Kalman Interpolated',markeredgewidth=0.05)
				xcoord = list(df[df.PosQC=='Acoustic Tracked'].Lon.values)
				ycoord = list(df[df.PosQC=='Acoustic Tracked'].Lat.values)
				xpos,ypos = self(xcoord,ycoord)
				self.plot(xpos,ypos,'r^',markersize=markersz+1,label='RAFOS Enabled',markeredgewidth=0.1)
			except TypeError:
				pass
			xcoord = list(df[df.PosQC==1].Lon.values)
			ycoord = list(df[df.PosQC==1].Lat.values)
			xpos,ypos = self(xcoord,ycoord)
			self.plot(xpos,ypos,'bo',markersize=markersz+2,label='Satellite Tracked',markeredgewidth=0.1)
			xcoord = list(df[df.PosQC==8].Lon.values)
			ycoord = list(df[df.PosQC==8].Lat.values)
			xpos,ypos = self(xcoord,ycoord)
			self.plot(xpos,ypos,'go',markersize=markersz,label='Linearly Interpolated',markeredgewidth=0.1)
			if self.date:
				xcoord = list(df[df.Date==self.date].Lon.values)
				ycoord = list(df[df.Date==self.date].Lat.values)
				xpos,ypos = self(xcoord,ycoord)
				self.plot(xpos,ypos,'c*',markersize=markersz+12,label='Current Position',markeredgewidth=0.5,zorder=10)

	def toa_to_dist(self,toa):
		dist = toa/slow
		return dist

	def linespace(self,line_space=10,fontsz=10):
		parallels = np.arange(-90,0,float(line_space)/2)
		self.drawparallels(parallels,labels=[1,0,0,0],fontsize=fontsz)
		meridians = np.arange(-360.,360.,float(line_space))
		self.drawmeridians(meridians,labels=[0,0,0,1],fontsize=fontsz)

 
	def annotate(self):
		df_annotate = self.dataframe.drop_duplicates(subset=['Cruise'])
		namelist = df_annotate.Cruise.values
		xcoord = list(df_annotate.Lon.values)
		ycoord = list(df_annotate.Lat.values)
		xpos,ypos = self(xcoord,ycoord)
		x2 = 10*np.sin(np.linspace(-np.pi/2,np.pi/2,len(namelist)))
		y2 = 10*np.cos(np.linspace(-np.pi/2,np.pi/2,len(namelist)))
		for k in range(len(namelist)):
			plt.annotate(namelist[k], xy=(xpos[k], ypos[k]),  xycoords='data',xytext=(x2[k], y2[k]), textcoords='offset points',color='midnightblue',fontsize=10,fontweight='bold')


	def shoot(self, lon, lat, azimuth, maxdist=None):
		"""Shooter Function
		Original javascript on http://williams.best.vwh.net/gccalc.htm
		Translated to python by Thomas Lecocq
		"""
		glat1 = lat*np.pi / 180.
		glon1 = lon*np.pi / 180.
		s = maxdist / 1.852
		faz = azimuth * np.pi / 180.
		EPS= 0.00000000005
		if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
			alert("Only N-S courses are meaningful, starting at a pole!")
		a=6378.13/1.852
		f=1/298.257223563
		r = 1 - f
		tu = r * np.tan(glat1)
		sf = np.sin(faz)
		cf = np.cos(faz)
		if (cf==0):
			b=0.
		else:
			b=2. * np.arctan2 (tu, cf)
		cu = 1. / np.sqrt(1 + tu * tu)
		su = tu * cu
		sa = cu * sf
		c2a = 1 - sa * sa
		x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
		x = (x - 2.) / x
		c = 1. - x
		c = (x * x / 4. + 1.) / c
		d = (0.375 * x * x - 1.) * x
		tu = s / (r * a * c)
		y = tu
		c = y + 1
		while (np.abs (y - c) > EPS):
			sy = np.sin(y)
			cy = np.cos(y)
			cz = np.cos(b + y)
			e = 2. * cz * cz - 1.
			c = y
			x = e * cy
			y = e + e - 1.
			y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
				  d / 4. - cz) * sy * d + tu
		b = cu * cy * cf - su * sy
		c = r * np.sqrt(sa * sa + b * b)
		d = su * cy + cu * sy * cf
		glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
		c = cu * cy - su * sy * cf
		x = np.arctan2(sy * sf, c)
		c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
		d = ((e * cy * c + cz) * sy * c + y) * sa
		glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
		baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
		glon2 *= 180./np.pi
		glat2 *= 180./np.pi
		baz *= 180./np.pi
	 
		return (glon2, glat2, baz)
	 
	def equi(self,centerlon, centerlat, radius, *args, **kwargs):
		glon1 = centerlon
		glat1 = centerlat
		X = []
		Y = []
		for azimuth in range(-180, 180):
			glon2, glat2, baz = self.shoot(glon1, glat1, azimuth, radius)
			X.append(glon2)
			Y.append(glat2)
		X.append(X[0])
		Y.append(Y[0])
		#~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
		if self.lon_0 == 180:
			X = np.array(X)+180

		X,Y = self(X,Y)
		self.plot(X,Y,**kwargs)


	def ellipse(self, x0, y0, a, b, n, phi=0,ax=None, **kwargs):
		"""
		Draws a polygon centered at ``x0, y0``. The polygon approximates an
		ellipse on the surface of the Earth with semi-major-axis ``a`` and 
		semi-minor axis ``b`` degrees longitude and latitude, made up of 
		``n`` vertices.

		For a description of the properties of ellipsis, please refer to [1].

		The polygon is based upon code written do plot Tissot's indicatrix
		found on the matplotlib mailing list at [2].

		Extra keyword ``ax`` can be used to override the default axis instance.

		Other \**kwargs passed on to matplotlib.patches.Polygon

		RETURNS
			poly : a maptplotlib.patches.Polygon object.

		REFERENCES
			[1] : http://en.wikipedia.org/wiki/Ellipse


		"""
		# if (self.lon_0==0)&(x0>180):
		#     x0 = x0-360
		# if (self.lon_0==180)&(x0>180):
		#     x0 = x0-720
			# print 'I have subtracted'
		print(self.lon_0)
		a = a/abs(np.cos(np.deg2rad(y0)))
		ax = kwargs.pop('ax', None) or self._check_ax()
		g = pyproj.Geod(a=self.rmajor, b=self.rminor)
		# Gets forward and back azimuths, plus distances between initial
		# points (x0, y0)
		azf, azb, dist = g.inv([x0, x0], [y0, y0], [x0+a, x0], [y0, y0+b])
		tsid = dist[0] * dist[1] # a * b

		# Initializes list of segments, calculates \del azimuth, and goes on 
		# for every vertex
		seg = []
		AZ = np.linspace(azf[0], 360. + azf[0], n)
		for i, az in enumerate(AZ):
			# Skips segments along equator (Geod can't handle equatorial arcs).
			# az =+ phi
			if np.allclose(0., y0) and (np.allclose(90., az) or
				np.allclose(270., az)):
				continue
			# In polar coordinates, with the origin at the center of the 
			# ellipse and with the angular coordinate ``az`` measured from the
			# major axis, the ellipse's equation  is [1]:
			#
			#                           a * b
			# r(az) = ------------------------------------------
			#         ((b * cos(az))**2 + (a * sin(az))**2)**0.5
			#
			# Azymuth angle in radial coordinates and corrected for reference
			# angle.
			azr = 2. * np.pi / 360. * (az+phi)
			A = dist[0] * np.sin(azr)
			B = dist[1] * np.cos(azr)
			r = tsid / (B**2. + A**2.)**0.5
			lon, lat, azb = g.fwd(x0, y0, az, r)
			x, y = self(lon, lat)
			# Add segment if it is in the map projection region.
			if x < 1e20 and y < 1e20:
				seg.append((x, y))
		segx = np.array(zip(*seg)[0])
		if not (all(item >= 0 for item in segx) or all(item < 0 for item in segx)):
			return
		segy = zip(*seg)[1]
		#     segx = [abs(x) for x in segx]
		seg = zip(segx,segy)
		poly = Polygon(seg, **kwargs)
		ax.add_patch(poly)
		# Set axes limits to fit map region.
		self.set_axes_limits(ax=ax)
		return poly


			# az +=phi
			# # Skips segments along equator (Geod can't handle equatorial arcs).
			# if np.allclose(0., y0) and (np.allclose(90., az) or
			#     np.allclose(270., az)):
			#     continue
			# # In polar coordinates, with the origin at the center of the 
			# # ellipse and with the angular coordinate ``az`` measured from the
			# # major axis, the ellipse's equation  is [1]:
			# #                           a * b
			# # r(az) = ------------------------------------------
			# #         ((b * cos(az))**2 + (a * sin(az))**2)**0.5
			# # Azymuth angle in radial coordinates and corrected for reference
			# # angle.
			# azr = 2. * np.pi / 360. * (az + 90.)
			# A = dist[0] * np.sin(azr)
			# B = dist[1] * np.cos(azr)
			# r = tsid / (B**2. + A**2.)**0.5
			# lon, lat, azb = g.fwd(x0, y0, az, r)
			# x, y = self(lon, lat)


	def ice_map(self,mean_month=False,col='crimson',line=8):
		di = pd.read_csv(os.path.join(datadir,'ALLICE.csv'),index_col = 'Unnamed: 0')
		di.Date = pd.to_datetime(di.Date)
		if mean_month is not False:
			date = np.sort(di[di.Date<=datetime.date(1986,mean_month,1)]['Date'].unique())[-1]
		else:    
			date = np.sort(di[di.Date<self.date]['Date'].unique())[-1]
			print(date)
		if self.lon_0==0:
			di = di[(di.Date == date)&(di.Lon>self.lonmin)&(di.Lon<self.lonmax)]
		elif self.lon_0==180:
			di.loc[di.Lon.values<0,['Lon']] = di[di.Lon<0].Lon+360
			di = di[(di.Date == date)&(di.Lon>self.lonmin)&(di.Lon<self.lonmax)]        
		for shape in di.Shape.unique():
			di_working = di[di.Shape==shape]
			try: self.plot(di_working.Lon.values,di_working.Lat.values,linewidth=line,color=col,latlon=True)
			except IndexError: 
				continue

	def streamline_plot(self):
		mat = scipy.io.loadmat(os.path.join(datadir,'agva_mean_streamfunction_1000dbar.mat'))
		XC = mat['long'][:][0]
		XC = np.array(XC.tolist())
		XC[XC>290] = XC[XC>290]-360
		YC = mat['latg'][:][0]
		streamline=mat['mean_streamfunction_1000dbar']


		streamline_mask, XC = shiftgrid(180, streamline, XC, start=False)
		XC,YC = np.meshgrid(XC,YC)
		XX,YY = self(XC,YC)

		levels = np.arange(np.nanmin(streamline_mask),np.nanmax(streamline_mask),15)
		self.contour(XX,YY,streamline_mask,linewidths=2,colors='orange',animated=True,levels = levels, label='Streamlines',alpha=0.5)

	def f_over_h(self, type='streamline',subsample=10):
		subsamplex = subsample
		subsampley = subsample
		mat = scipy.io.loadmat(os.path.join(datadir,'grid.mat'))
		depth = mat['Depth'].T

		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		if self.lon_0==0:
			depth, XC = shiftgrid(180, depth, XC, start=False)
		depth = gaussian_filter(depth,7)
		depth = np.ma.array(depth,mask=depth<500)
		XC,YC = np.meshgrid(XC,YC)
		omega = 7.2921 * 10**-5
		f = 2*omega*np.sin(np.deg2rad(YC))
		out = f/depth
		XX,YY = self(XC,YC)

		if type=='streamline':
			# mean = out.mean()
			# std = out.std()
			# outmin = mean-std/2
			# outmax = mean+std/2
			# print outmin 
			# print outmax
			# out = np.ma.array(out,mask=(out>outmax)|(out<outmin))
			levels = np.sort([(out.min()-out.max())/(float(x)+1.) + out.max() for x in range(20)])
			# print levels
			cmap='inferno'
			CS = self.contour(XX,YY,out,levels=levels,colors = plt.get_cmap('Accent')(range(20)))
			# plt.clabel(CS,fmt='%.6e', fontsize=9, inline=1)

		elif type=='quiver':
			U,V = np.gradient(out)
			U = -U
			V = V    #because of the definition of the streamfunction
			self.quiver(XX[::subsamplex,::subsampley],YY[::subsamplex,::subsampley],U[::subsamplex,::subsampley],V[::subsamplex,::subsampley])

		elif type=='gradient':
			out = nabla(out)
			levels = np.sort([(out.min()-out.max())/(float(x)+1.) + out.max() for x in range(20)])
			# print levels
			cmap='inferno'
			CS = self.contour(XX,YY,out,levels=levels,colors = plt.get_cmap('Accent')(range(20)))

	def buoyancy_plot(self):
		file = os.path.join(datadir,'averageBF.npy')
		mat = scipy.io.loadmat(os.path.join(datadir,'grid.mat'))
		b = np.load(file)
		b_mask = np.ma.array(b,mask=((b==0.)|(b==1.)))
		b_mask[b_mask<-200]=-200
		b_mask[b_mask>200]=200
		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		if self.lon_0==0:
			b_mask, XC = shiftgrid(180, b_mask, XC, start=False)
		XC,YC = np.meshgrid(XC,YC)
		XX,YY = self(XC,YC)
		self.pcolormesh(XX[::self.subsample,::self.subsample],YY[::self.subsample,::self.subsample],b_mask[::self.subsample,::self.subsample],cmap=plt.cm.PuOr,alpha=0.4)
		plt.colorbar(label='Buoyancy Flux ($W\ m^{-2}$)',fraction=0.046, pad=0.04)

	def plot_features(self,fontsz=16):
		for feature in map_features: 
			XX,YY = self(feature[1],feature[2])
			plt.annotate(feature[0], xy=(XX,YY),  xycoords='data',xytext=(XX,YY),color='teal',fontsize=fontsz,fontweight='bold') # ,arrowprops=dict(arrowstyle="fancy", color='b'))

	def co2flux_plot(self):
		df = pd.read_csv(os.path.join(datadir,'sumflux_2006c.txt'), skiprows = 51,sep=r"\s*")
		x = np.arange(df.LON.min(),df.LON.max()+5,5)
		y = np.arange(df.LAT.min(),df.LAT.max()+4,4)
		XC,YC = np.meshgrid(x,y)
		CO2 = np.zeros([len(y),len(x)])
		di = df.iterrows()
		for i in range(len(df)):
			row = next(di)[1]
			CO2[(row['LON']==XC)&(row['LAT']==YC)] = row['FLUXGMSW06']
		if self.lon_0==0:
			CO2, x = shiftgrid(180, CO2, x, start=False)
		CO2_mask = np.ma.masked_equal(CO2,0.)
		XC,YC = np.meshgrid(x,y)
		XX,YY = self(XC,YC)
		self.pcolormesh(XX,YY,CO2_mask,cmap=plt.cm.PRGn)
		plt.colorbar(label='CO2 Flux $gm C/m^2/yr$',fraction=0.046, pad=0.04)

	def temp_plot(self,level='bottom',colormap=plt.cm.Greys_r):
		if level is 'top':
			level=0
		elif level is 'bottom':
			level = 18
		temp = np.load(os.path.join(datadir,'averageT_'+str(level)+'.npy'))
		mat = scipy.io.loadmat(os.path.join(datadir,'grid.mat'))
		temp_mask = np.ma.masked_equal(temp,0.)
		temp_mask[temp_mask>2]=2
		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		if self.lon_0==0:
			temp_mask, XC = shiftgrid(180, temp_mask, XC, start=False)
		XC,YC = np.meshgrid(XC,YC)
		XX,YY = self(XC,YC)
		clevs = np.linspace(np.nanmin(temp_mask),np.nanmax(temp_mask),40)
		self.contourf(XX[::self.subsample,::self.subsample],YY[::self.subsample,::self.subsample],temp_mask[::self.subsample,::self.subsample],cmap=colormap,levels=clevs,animated=True,alpha=.75)
		plt.colorbar(label='Temperature $(^\circ C)$',fraction=0.046, pad=0.04,ticks=[-1.5, -.5, .5,1.5])

	def ssh_plot(self):
		ssh = np.load(os.path.join(datadir,'sshave.npy'))
		mat = scipy.io.loadmat(os.path.join(datadir,'grid.mat'))
		ssh_mask = np.ma.masked_equal(ssh,0.)
		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		if self.lon_0==0:
			ssh_mask, XC = shiftgrid(180, ssh_mask, XC, start=False)
		XC,YC = np.meshgrid(XC,YC)
		XX,YY = self(XC,YC)
		clevs = np.arange(-4,4,0.08)
		cs = self.contour(XX[::self.subsample,::self.subsample],YY[::self.subsample,::self.subsample],ssh_mask[::self.subsample,::self.subsample],clevs,linewidths=0.5,colors='k',animated=True)
		self.contourf(XX[::self.subsample,::self.subsample],YY[::self.subsample,::self.subsample],ssh_mask[::self.subsample,::self.subsample],clevs,cmap=plt.cm.RdBu_r,animated=True)
		plt.clabel(cs,inline=1,fontsize=6)


	def eaverage(self,level='bottom',day=0):
		subsample = self.subsample
		if level is 'top':
			level=0
		elif level is 'bottom':
			level = 23
		with open(os.path.join(datadir,'low_res_sose_mat.pickle'),'rb') as handle: 
			mat = pickle.load(handle)
		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		with open(os.path.join(datadir,'Eaverage_'+str(level)+'.pickle'),'rb') as handle: 
			plot_data = pickle.load(handle)    
		if self.lon_0==0:
			plot_data, XC = shiftgrid(180, plot_data, XC, start=False)
		XC,YC = np.meshgrid(XC,YC)
		XX,YY = self(XC,YC)
		self.contourf(XX[::subsample,::subsample],YY[::subsample,::subsample],plot_data[::subsample,::subsample],cmap=plt.cm.Greys,norm = LogNorm(),animated=True)
		plt.colorbar(label='J/kg',fraction=0.046, pad=0.04)



	def eke(self,level='bottom',day=0):
		subsample = self.subsample
		if level is 'top':
			level=0
		elif level is 'bottom':
			level = 23
		with open(os.path.join(datadir,'low_res_sose_mat.pickle'),'rb') as handle: 
			mat = pickle.load(handle)
		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		with open(os.path.join(datadir,'EKE_'+str(level)+'.pickle'),'rb') as handle: 
			data = pickle.load(handle)    
		plot_data = data[day,:,:]
		if self.lon_0==0:
			plot_data, XC = shiftgrid(180, plot_data, XC, start=False)
		XC,YC = np.meshgrid(XC,YC)
		XX,YY = self(XC,YC)
		self.contourf(XX[::subsample,::subsample],YY[::subsample,::subsample],plot_data[::subsample,::subsample],cmap=plt.cm.Purples,norm = LogNorm(),animated=True,alpha=0.6)
		plt.colorbar(label='Mean Eddy Kinetic Energy (J/kg)',fraction=0.046, pad=0.04,shrink=0.5)



	def etotal(self,level='bottom',day=0):
		subsample = self.subsample
		if level is 'top':
			level=0
		elif level is 'bottom':
			level = 23
		with open(os.path.join(datadir,'low_res_sose_mat.pickle'),'rb') as handle: 
			mat = pickle.load(handle)
		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		with open(os.path.join(datadir,'Etotal_'+str(level)+'.pickle'),'rb') as handle: 
			data = pickle.load(handle)    
		plot_data = data[day,:,:]
		if self.lon_0==0:
			plot_data, XC = shiftgrid(180, plot_data, XC, start=False)
		XC,YC = np.meshgrid(XC,YC)
		XX,YY = self(XC,YC)
		self.contourf(XX[::subsample,::subsample],YY[::subsample,::subsample],plot_data[::subsample,::subsample],cmap=plt.cm.Greys,norm = LogNorm(),animated=True)
		plt.colorbar(label='J/kg',fraction=0.046, pad=0.04)

	def bathy(self,color=plt.cm.Oranges,contour=False):
		cdict = {'red':  ((0.0, 0.63, 0.63),
				(0.2, 0.66, 0.66),
				   (0.4, 0.70, 0.70),
				   (0.6, .78, .78),
				   (0.8, 0.92, 0.92),
				   (1.0, 0.98, 0.98)),

		 'green':((0.0, 0.63, 0.63),
				(0.2, 0.66, 0.66),
				   (0.4, 0.70, 0.70),
				   (0.6, .78, .78),
				   (0.8, 0.92, 0.92),
				   (1.0, 0.98, 0.98)),

		 'blue': ((0.0, 0.63, 0.63),
				(0.2, 0.66, 0.66),
				   (0.4, 0.70, 0.70),
				   (0.6, .78, .78),
				   (0.8, 0.92, 0.92),
				   (1.0, 0.98, 0.98)),
		}
		bathy = LinearSegmentedColormap('bathy', cdict)
		subsample = self.subsample
		with open(os.path.join(datadir,'low_res_sose_mat.pickle'),'rb') as handle: 
			mat = pickle.load(handle)
		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		plot_data = mat['Depth'].T/1000.
		if self.lon_0==0:
			plot_data, XC = shiftgrid(180, plot_data, XC, start=False)
		XC,YC = np.meshgrid(XC,YC)
		XX,YY = self(XC,YC)
		levels = [0,1,2,3,4,5]
		if not contour:
			self.contourf(XX[::subsample,::subsample],YY[::subsample,::subsample],plot_data[::subsample,::subsample],levels,cmap=bathy,animated=True,vmax=5,vmin=0)
			cbar = self.colorbar()
			cbar.set_label('Depth (km)')
		else:
			self.contour(XX[::subsample,::subsample],YY[::subsample,::subsample],plot_data[::subsample,::subsample],cmap=plt.get_cmap('Greys'),alpha=.8,vmax=5,vmin=-2)


	def quiver_plot(self,level='bottom'):
		subsample = self.subsample*4+5
		if level is 'top':
			level=0
		elif level is 'bottom':
			level = 23
		U = np.load(os.path.join(datadir,'averageU_'+str(level)+'.npy'))
		V = np.load(os.path.join(datadir,'averageV_'+str(level)+'.npy'))
		mat = scipy.io.loadmat(os.path.join(datadir,'grid.mat'))
		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		if self.lon_0==0:
			U, dummy = shiftgrid(180, U, XC, start=False)
			V, XC = shiftgrid(180, V, XC, start=False)
		XC,YC = np.meshgrid(XC,YC)
		XX,YY = self(XC,YC)
		V_mask = np.ma.masked_equal(V,0.)
		U_mask = np.ma.masked_equal(U,0.)
		self.quiver(XX[::subsample,::subsample],YY[::subsample,::subsample],U_mask[::subsample,::subsample],V_mask[::subsample,::subsample],linewidths=0.5)

	def SOSO(self,markersz = 15,region=None):
		if region == 'weddell':
			x = [a[1] for a in SOSO_weddel_coord.values()]
			y = [a[0] for a in SOSO_weddel_coord.values()]
			xpos, ypos = self(x,y)
			self.plot(xpos,ypos,'*',markersize = markersz, color = 'gold',label = 'Weddell Sound Sources')
		if region == 'dimes':
			x1 = [a[1] for a in SOSO_drake_coord.values()]
			y1 = [a[0] for a in SOSO_drake_coord.values()]
			xpos, ypos = self(x1,y1)
			self.plot(xpos,ypos,'s',markersize = markersz/2, color = 'gold',label = 'Dimes Sound Sources') #,label = 'Sound Sources')

	def LOP_plot(self):
		soso_df = self.soso_df[(self.soso_df.Date==self.date)&(self.soso_df.Cruise==self.cruise)]
		for it in soso_df.iterrows():
			if it[1]['PosQC'] in SOSO_weddel_coord: 
				# print SOSO_weddel_coord[it[1]['PosQC']][0]
				# print SOSO_weddel_coord[it[1]['PosQC']][1]
				self.equi(SOSO_weddel_coord[it[1]['PosQC']][1],SOSO_weddel_coord[it[1]['PosQC']][0], self.toa_to_dist(it[1]['Observation']),lw=2.,label='source '+it[1]['PosQC'])
			if it[1]['PosQC'] in SOSO_drake_coord: 
				self.equi(self, SOSO_drake_coord[it[1]['PosQC']][1],SOSO_drake_coord[it[1]['PosQC']][0], self.toa_to_dist(it[1]['Observation']),lw=2.,label='source '+it[1]['PosQC'])

	def covariance_plot(self):
		w = self.dataframe[self.dataframe.Date==self.date].eig1.values[0]
		v = self.dataframe[self.dataframe.Date==self.date].eig2.values[0]
		print(w)
		print(v)

		xcoord = self.dataframe[self.dataframe.Date==self.date].Lon.values[0]
		ycoord = self.dataframe[self.dataframe.Date==self.date].Lat.values[0]
		angle = np.degrees(np.arctan(v[1,np.argmax(w)]/v[0,np.argmax(w)]))
		print(angle)
		poly = self.ellipse(xcoord, ycoord, 2*max(w)*np.sqrt(5.991),2*min(w)*np.sqrt(5.991),400, phi=angle, facecolor='c', zorder=10,
			alpha=0.3)
		
	def region_plot(self,region=[]):
		regions=[(-60,20,-76,-60,'weddell'),(140,280,-76,-60,'ross')]
		if region:
			regions = [x for x in regions if x[4] in region]
		for region in regions:
			lons = np.array([region[1],region[1]]+np.linspace(region[1],region[0],200).tolist()+[region[0],region[0],region[1]])
			lons[lons>180] = lons[lons>180]-360
			lats = [region[2],region[3]]+([region[3]]*200)+[region[3],region[2],region[2]]
			self.plot(lons.tolist(),lats,'k--',latlon=True,zorder=11,linewidth=7)

	def a12_stations(self,markersz=6):
		file = os.path.join(datadir,'a12.header')
		df = pd.read_csv(file,sep=r"\s+",skiprows=(0,1),header=None,usecols=(2,3),names=('Lat','Lon'))
		xcoord = list(df.Lon.values)
		ycoord = list(df.Lat.values)
		xpos,ypos = self(xcoord,ycoord)
		self.scatter(xpos,ypos,color='g',zorder=3)

	def polarstern_stations(self):
		file = os.path.join(datadir,'polarstern_ant302_locations_latlon.txt')
		df = pd.read_csv(file,sep=r"\s+",skiprows=(0,1),header=None,usecols=(0,1),names=('Lat','Lon'))
		xcoord = list(df.Lon.values)
		ycoord = list(df.Lat.values)
		xpos,ypos = self(xcoord,ycoord)
		self.scatter(xpos,ypos,color='r',zorder=4)

	def orsi_fronts(self):
		file_list = ['saf.asc','pf.asc','stf.asc','sbdy.asc','saccf.asc']
		for name in file_list:
			file = os.path.join(datadir,name)
			token = open(file,'r')
			coord = [i.strip().split() for i in token.readlines()]
			xcoord,ycoord = zip(*[(float(i[0]),float(i[1])) for i in coord if abs(float(i[0]))<179])
			xpos,ypos = self(xcoord,ycoord)
			self.plot(xpos,ypos,color='k',zorder=2,linewidth=1)
	def show(self):
		plt.show()
