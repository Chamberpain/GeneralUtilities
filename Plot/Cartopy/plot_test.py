# import sys,os
# sys.path.append(os.path.abspath("../scripts"))
# sys.path.append(os.path.abspath("/Users/paulchamberlain/git/chamberpain-working/scripts"))
# sys.path.append(os.path.abspath("/Users/paulchamberlain/Data/"))
# import soccom_proj_settings
from eulerian_plot import SBasemap
import matplotlib.pyplot as plt
import pandas as pd
import datetime
file = '/Users/paulchamberlain/Data/ARGO/interpolated_argo_all.pickle'
df = pd.read_pickle(file)
df = df[df.PosQC==1]


df = df[df.Cruise==5904153]
print df 

for region in ['polar','atlantic','pacific','weddell','ross','float']:
# for region in ['float']:

	print region
	plt.figure()
	a = SBasemap(region=region,gridlines=False,dataframe=df)
	a.streamline_plot()
	a.drifter_plot(dataframe=df)
	plt.title(region)
	print 'I made it to the end of the loop'
# for floater in df.Cruise.unique()[::50]:
# 	print floater
# 	plt.figure()
# 	frame = df[df.Cruise==floater]


# 	print frame.Lon.min()
# 	print frame.Lon.max()
# 	a.streamline_plot()		
# a.bathy()
# a.delta_bathy()
# a.streamline_plot()
# a.f_over_h()
# a.plot()
# a.ice_plot(date)
# a.bathy(subsample=50)
# plt.figure()
# frame = pd.DataFrame({'Lat':[-50,-52,0,-89],'Lon':[0.5,0.1,359,358],'Cruise':['Cow','Cow','Sheep','Sheep']})
# a.drifter_plot(dataframe=frame,lineplot=True)
# frame = pd.DataFrame({'Lat':[-50,-55,-51,-52,0,-89],'Lon':[0.5,90,180,181,359,358],'Cruise':['Cow','Cow','Sheep','Sheep','Horse','Horse']})
# a = SBasemap(dataframe=frame,region='float')
# a.drifter_plot(dataframe=frame,lineplot=True)
# a.bathy(subsample=30)
# plt.figure()
# frame = pd.DataFrame({'Lat':[-50,-52,0,-89],'Lon':[0.5,0.1,359,358],'Cruise':['Cow','Cow','Sheep','Sheep']})
# a = SBasemap(dataframe=frame,region='float')
# a.drifter_plot(dataframe=frame,lineplot=True)
# a.bathy(subsample=30)
# plt.figure()
# a = SBasemap(region='atlantic')
# a.bathy(subsample=30)
# plt.figure()
# a = SBasemap(region='pacific')
# a.bathy(subsample=30)
# plt.figure()
# a = SBasemap(region='weddell')
# a.bathy(subsample=30)
# plt.figure()
# a = SBasemap(region='ross')
# a.bathy(subsample=30)
plt.show()



# import pandas as pd
# import matplotlib.pyplot as plt
# import numpy as np
# import datetime
# This testing routine should be made into a unit test, will do this when I get a little time

# def ellipse_test(m):
# 	for lon in np.linspace(-180,180.,20):
# 		print 'lon is ',lon
# 		print 'm(lon) is ',m(lon,-60)
# 		y0 = -60
# 		x0 = lon
# 		# if x0 > 0:
# 		# x0 = x0+360
# 		# x0,y0 = m(x0,y0)
# 		a = 3
# 		b = 3
# 		n = 100
# 		m.ellipse(x0,y0,a,b,n,phi=np.random.randn()*360)
# 	return m

# # def circle_test(m)


# df = pd.read_pickle(soccom_proj_settings.interpolated_drifter_file)
# df['Lon180'] = df.Lon
# df.loc[df.Lon180.values>180,['Lon180']] = df[df.Lon180>180].Lon180-360
# df = df.dropna(subset = ['Lat','Lon'])
# df.Date = pd.to_datetime(df.Date)
# for date in [[2013,1,1]]:
# 	# for frame in [df[(df.Lon180<-170)|(df.Lon180>170)],df[(df.Lon<10)|(df.Lon>350)]]:
# 	# 	plt.figure()
# 	# 	m = SBasemap(dataframe=frame,date=date,timedelta = 100)
# 	# 	# m.buoyancy_plot()
# 	# 	# m.co2flux_plot()
# 	# 	# m = ellipse_test(m)
# 	# 	# m.ice_map()
# 	# 	# m.streamline_plot()
# 	# 	plt.title('Dataframe tests')
# 	# 	m.plot_features()
# 	# 	# m.eke()
# 	# 	m.bathy()
# 	# 	m.drifter_plot(drift_type='Argo',markersz=20) 
# 	# 	# m.etotal()
# 	# 	# m.eaverage()
# 	for region in ['polar','atlantic','pacific']:
# 		print region
# 		plt.figure()
# 		m = SBasemap(region=region,date=date,timedelta = -100)
# 		ellipse_test(m)
# 		# for lon in np.linspace(0,360,12):
#   #  			m.equi(lon, -60,500)
# 		# m.plot_features(fontsz=8)
# 		# m.buoyancy_plot()
# 		# m.co2flux_plot()
# 		# m = ellipse_test(m)
# 		# m.bathy()
# 		# m.drifter_plot()
# 		# m.annotate()
# 		# m.SOCCOM()
# 		# m.eke()
# 		# m.streamline_plot()
# 		# m.ice_map()
# 		# m.plot_features()
# 		# m.etotal()
# 		# m.eaverage()
# 		plt.title(region)
# 	plt.show()