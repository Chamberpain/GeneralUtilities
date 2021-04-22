import sys,os
sys.path.append(os.path.abspath("./data"))
from mpl_toolkits.basemap import shiftgrid
from eulerian_plot import SBasemap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
import scipy.io
from scipy.ndimage.filters import gaussian_filter
import scipy.optimize as sopt
from mpl_toolkits.mplot3d import axes3d
import oceans
from LatLon import LatLon
from eulerian_plot import big_weddell_coord
import seaborn

debug = False
datadir = os.getenv("HOME")+'/iCloud/Data/Processed/eulerian_plot/basemap/'




class f_over_h_contour(object): 

	def __init__(self,smooth_level=5):
		self.m = SBasemap(region='big_weddell')
		mat = scipy.io.loadmat(os.path.join(datadir,'grid.mat'))
		depth = mat['Depth'].T
		XC = mat['XC'][:,0]
		YC = mat['YC'][0,:]
		# if m.lon_0==0:
		depth, XC = shiftgrid(180, depth, XC, start=False)
		depth = gaussian_filter(depth,smooth_level)

		if debug:
			print XC.min()
		if debug:
			print XC.max()
		XD,YD = np.meshgrid(XC,YC)
		depth = depth[(XD>big_weddell_coord[0])&(XD<big_weddell_coord[1])&(YD<big_weddell_coord[3])]
		rows =  YC[YC<big_weddell_coord[3]].shape[0]
		columns = XC[(XC>big_weddell_coord[0])&(XC<big_weddell_coord[1])].shape[0]
		depth = depth.reshape([rows,columns])
		self.XC = XD[(XD>big_weddell_coord[0])&(XD<big_weddell_coord[1])&(YD<big_weddell_coord[3])].reshape([rows,columns])
		self.YC = YD[(XD>big_weddell_coord[0])&(XD<big_weddell_coord[1])&(YD<big_weddell_coord[3])].reshape([rows,columns])
		self.Xm,self.Ym = self.m(self.XC,self.YC)
		omega = 7.2921 * 10**-5
		f = 2*omega*np.sin(np.deg2rad(self.YC))
		self.out = f/depth
		self.A,self.B = np.gradient(self.out)


	def plot_f_h_data(self):
		plt.figure()
		plt.plot(x2,y2,'ro')
		plt.plot(lon,lat)
		plt.plot(x_g,y_g)
		plt.plot(x_s, y_s,'ro')
		plt.plot(x_g[0],y_g[0],'ko')
		plt.show()

	def _rect_inter_inner(self,x1,x2):
	    n1=x1.shape[0]-1
	    n2=x2.shape[0]-1
	    X1=np.c_[x1[:-1],x1[1:]]
	    X2=np.c_[x2[:-1],x2[1:]]
	    S1=np.tile(X1.min(axis=1),(n2,1)).T
	    S2=np.tile(X2.max(axis=1),(n1,1))
	    S3=np.tile(X1.max(axis=1),(n2,1)).T
	    S4=np.tile(X2.min(axis=1),(n1,1))
	    return S1,S2,S3,S4

	def _rectangle_intersection_(self,x1,y1,x2,y2):
	    S1,S2,S3,S4=self._rect_inter_inner(x1,x2)
	    S5,S6,S7,S8=self._rect_inter_inner(y1,y2)

	    C1=np.less_equal(S1,S2)
	    C2=np.greater_equal(S3,S4)
	    C3=np.less_equal(S5,S6)
	    C4=np.greater_equal(S7,S8)

	    ii,jj=np.nonzero(C1 & C2 & C3 & C4)
	    return ii,jj

	def intersection(self,x1,y1,x2,y2):
	    """
	INTERSECTIONS Intersections of curves.
	   Computes the (x,y) locations where two curves intersect.  The curves
	   can be broken with NaNs or have vertical segments.
	usage:
	x,y=intersection(x1,y1,x2,y2)
	    Example:
	    a, b = 1, 2
	    phi = np.linspace(3, 10, 100)
	    x1 = a*phi - b*np.sin(phi)
	    y1 = a - b*np.cos(phi)
	    x2=phi
	    y2=np.sin(phi)+2
	    x,y=intersection(x1,y1,x2,y2)
	    plt.plot(x1,y1,c='r')
	    plt.plot(x2,y2,c='g')
	    plt.plot(x,y,'*k')
	    plt.show()
	    """
	    ii,jj=self._rectangle_intersection_(x1,y1,x2,y2)
	    n=len(ii)

	    dxy1=np.diff(np.c_[x1,y1],axis=0)
	    dxy2=np.diff(np.c_[x2,y2],axis=0)

	    T=np.zeros((4,n))
	    AA=np.zeros((4,4,n))
	    AA[0:2,2,:]=-1
	    AA[2:4,3,:]=-1
	    AA[0::2,0,:]=dxy1[ii,:].T
	    AA[1::2,1,:]=dxy2[jj,:].T

	    BB=np.zeros((4,n))
	    BB[0,:]=-x1[ii].ravel()
	    BB[1,:]=-x2[jj].ravel()
	    BB[2,:]=-y1[ii].ravel()
	    BB[3,:]=-y2[jj].ravel()

	    for i in range(n):
	        try:
	            T[:,i]=np.linalg.solve(AA[:,:,i],BB[:,i])
	        except:
	            T[:,i]=np.NaN


	    in_range= (T[0,:] >=0) & (T[1,:] >=0) & (T[0,:] <=1) & (T[1,:] <=1)

	    xy0=T[2:,in_range]
	    xy0=xy0.T
	    return xy0[:,0],xy0[:,1]

	def find_nearest(self,array,value):
		assert (type(value) in [float,int,np.float64]),'find nearest failed because value was not right format'
		diff_array = np.abs(array-value)
		assert (value >= array.min()), 'find nearest failed because value was less than array, this was '+str(value)
		assert (value <= array.max()), 'find nearest failed because value was greater than array'
		assert diff_array.min()<=np.diff(array).max()/2, 'find nearest failed because array difference was out of bounds '+str(diff_array.min())+' '+str(np.diff(array).max()/2)
		idx = diff_array.argmin()
		return idx

	def find_nearest_array(self,array_lat,array_lon,lat_,lon_):
		assert type(lat_) in [float,int,np.float64]
		assert type(lon_) in [float,int,np.float64]
		assert lon_<=180
		idx = (np.abs(np.array(array_lat)-lat_)+np.abs(np.array(array_lon)-lon_)).argmin()
		return idx

	def on_curve(self,array_lat,array_lon,lat_,lon_):
		assert type(lat_) in ([float,int,np.float64]),'find nearest array failed because lat was not right format'
		assert type(lon_) in ([float,int,np.float64]),'find nearest array failed because lon was not right format'
		assert (lon_<=180),'find nearest array failed because lon was out of bounds'
		tolerance = 0.5
		truth = min(np.abs(np.array(array_lat)-lat_)+np.abs(np.array(array_lon)-lon_))<tolerance
		if debug:
			if truth:
				'on curve'
			else:
				'not on curve'
		return truth

	def distance_calculator(self,array_lat,array_lon):
		distance_list = [LatLon(i[0],i[1]) for i in zip(array_lat,array_lon)]
		calc_list = []
		for i in range(len(distance_list)-1):
			try:
				calc_list.append((distance_list[i]-distance_list[i+1]).magnitude)
			except ValueError:
				raise
		return np.sum(calc_list)


	def PV(self,x):
		""" accepts as input a LatLon object and returns the PV at that location"""
		assert (x.lon.decimal_degree<=180),'PV failed because lon was out of bounds'
		column = self.find_nearest(self.XC[0,:],x.lon.decimal_degree)
		row = self.find_nearest(self.YC[:,0],x.lat.decimal_degree)
		return self.out[row,column]

	def dPV(self,x):
		""" accepts as input a Lat,Lon tuple and returns the differental PV in both directions"""
		assert x[1]<=180,'dPV failed because lon was out of bounds'
		column = self.find_nearest(self.XC[0,:],x[1])
		row = self.find_nearest(self.YC[:,0],x[0])
		return np.array([self.A[row,column], self.B[row,column]])

	def f_over_h_calculator(self,pos1,pos2):
		""" Takes as input 2 positions in LatLon formatand calculates the shortest 
		f/h following contour between the two"""
		#this establishes a wrap 180 coordinate system
		assert (pos1.lon.decimal_degree<=180),'pos1 lon was out of bounds'
		assert (pos2.lon.decimal_degree<=180),'pos2 lon was out of bounds'
		plt.close()	#we do this to make the algorithm run efficiently, otherwise it bogs down with all the points
		pv_end = self.PV(pos2)	# get the potential vorticity at the end
		pv_start = self.PV(pos1)  # get the potential vorticity at the beginning
		CS = self.m.contour(self.Xm,self.Ym,self.out,pv_start) # get the contour of points which follow the pv_start geometry
		guesses = [np.array([pos2.lat.decimal_degree,pos2.lon.decimal_degree ])]	#guesses is iteratively added to, but the ending pv is a good first guess
		direction = np.sign(pv_start-pv_end)	#determine direction of the gradient
		for dummy in range(300):	#iterate 300 times
			x = guesses[-1]			#take the last guess in the list 
			s = direction*self.dPV(x)	#this is the gradient of PV at the location
			alpha = .05/ np.sqrt(s[0]**2+s[1]**2) #this is the total slope 
			next_guess = x + alpha * s #iterate in the direction of the slope
			if next_guess[0]<self.YC[:,0].min():
				break
			if next_guess[0]>self.YC[:,0].max():
				break
			if next_guess[1]<self.XC[0,:].min():
				break
			if next_guess[1]>self.XC[0,:].max():
				break			
			if np.isnan(next_guess).any():
				break
			else:
				guesses.append(next_guess)
		y_g,x_g = zip(*guesses)
		y_s,x_s = (pos1.lat.decimal_degree,pos1.lon.decimal_degree)
		p = [element.vertices.tolist() for element in CS.collections[0].get_paths()[:]]
		_didnt_work = 1
		_along_f_over_h = np.nan
		_across_f_over_h = np.nan
		_along_pos = np.nan
		_across_pos = np.nan
		for item in p:
			x,y = zip(*item)
			lon, lat = self.m(x,y,inverse=True)
			if debug:
				print 'there was an element in the vertices'
			if self.on_curve(np.array(lat),np.array(lon),y_s,x_s):
				if debug:
					print 'there was an element on the curve'


				idx2 = self.find_nearest_array(lat,lon,y_s,x_s)

				start_diff = LatLon(y_s,x_s) - LatLon(lat[idx2],lon[idx2])

				pos_list = [LatLon(i[0],i[1])+start_diff for i in zip(lat,lon)]
				lat,lon = zip(*[(i.lat.decimal_degree,i.lon.decimal_degree) for i in pos_list])

				x2,y2 = self.intersection(np.array(x_g),np.array(y_g),np.array(lon),np.array(lat))				
				if x2.any():						
					x2 = x2[0]; y2= y2[0]
					idx1 = self.find_nearest_array(lat,lon,y2,x2)
					idx3 = self.find_nearest_array(y_g,x_g,y2,x2)
					start_idx = min(idx1,idx2)
					end_idx = max(idx1,idx2)
					sign_idx = np.sign(idx1-idx2)
					if sign_idx ==0: 
						sign_idx = 1
					if not sign_idx:
						start_idx -= 1
					else:
						end_idx += 1
					_along_pos = zip(lat[start_idx:end_idx][::sign_idx],lon[start_idx:end_idx][::sign_idx])
					_across_pos = zip(y_g[:idx3],x_g[:idx3])

					_along_f_over_h = self.distance_calculator(lat[start_idx:end_idx],lon[start_idx:end_idx])
					_across_f_over_h = self.distance_calculator(y_g[:idx3],x_g[:idx3])
					_didnt_work = 0
					if debug:
						print 'found one that was on curve with an intersection'
						plt.close()
						plt.plot(lon,lat)
						plt.plot(x_g,y_g)
						plt.plot(x_s,y_s,'ro',markersize=20)
						x_e,y_e = (pos2.lon.decimal_degree,pos2.lat.decimal_degree)
						plt.plot(x_e,y_e,'y*',markersize=20)
						plt.plot(x2,y2,'b*',markersize=20)
						print min(lon)
						print max(lon)
						print min(lat)
						print max(lat)
						print x_s
						print y_s
						plt.show()
				else:
					if debug:
						print 'there was no intersection'
		return _didnt_work,_along_f_over_h,_across_f_over_h,_along_pos,_across_pos


	def test1(self):
		df = pd.read_pickle(soccom_proj_settings.interpolated_drifter_file)
		df['Lon'] = oceans.wrap_lon180(df.Lon)
		df = df.dropna(subset = ['Lat','Lon'])
		plt.figure(figsize=(7,6))
		self.m.drifter_plot(lineplot=True,color='r',drift_type=8,dataframe=df)
		self.m.f_over_h()
		self.m.show()

	def test2(self):
		df = pd.read_pickle(soccom_proj_settings.interpolated_drifter_file)
		df = df[df.Lat<-60]
		df = df[(df.Lon<20)|(df.Lon>300)]
		for pv in df.PV[::500]:
			CS = self.m.contour(self.Xm,self.Ym,self.out,pv)
			p = [element.vertices.tolist() for element in CS.collections[0].get_paths()[:]]  
			plt.figure()
			for item in p:
				x,y = zip(*item)
				lon, lat = self.m(x,y,inverse=True)
				plt.plot(lon,lat) 
			plt.show()

	def test3(self):
		interpolated_drifter_file = os.getenv("HOME")+'/iCloud/Data/Processed/observing_ice_covered/interpolated_drifter.pickle'
		dataf = pd.read_pickle(interpolated_drifter_file)
		dataf = dataf[dataf.PosQC==1]
		dataf['Lon'] = oceans.wrap_lon180(dataf.Lon)
		dataf = dataf[dataf.Lat<-60]
		dataf = dataf[(dataf.Lon>-60)&(dataf.Lon<20)]
		along_f_over_h = []
		across_f_over_h = []
		pv_start = []
		pv_end = []
		start = []
		end = []
		didnt_work = []
		for name in dataf.Cruise.unique():
			mask = dataf.Cruise==name
			dataf.loc[mask,'Diff'] = dataf[mask]['Date'].diff().dt.days>30
			dataf.loc[mask,'Diff'] = dataf[mask].Diff.apply(lambda x: 1 if x else 0).cumsum()
			for g in dataf[mask].groupby('Diff').groups:
				frame = dataf[mask].groupby('Diff').get_group(g)
				didnt_work_temp, along_f_over_h_temp, across_f_over_h_temp,dummy1,dummy2 = self.f_over_h_calculator(frame)
				along_f_over_h.append(along_f_over_h_temp)
				across_f_over_h.append(across_f_over_h_temp)
				didnt_work.append(didnt_work_temp)
		plt.figure()
		seaborn.jointplot(np.array(along_f_over_h),np.array(across_f_over_h),kind='reg').set_axis_labels('along f/h','across f/h')
		plt.show()