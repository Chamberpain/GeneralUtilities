import sys,os
sys.path.append(os.path.abspath("../../"))
import soccom_proj_settings
from eulerian_plot import SBasemap
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import datetime
# This testing routine should be made into a unit test, will do this when I get a little time

def ellipse_test(m):
	for lon in np.linspace(-180,180.,20):
		print 'lon is ',lon
		print 'm(lon) is ',m(lon,-60)
		y0 = -60
		x0 = lon
		# if x0 > 0:
		# x0 = x0+360
		# x0,y0 = m(x0,y0)
		a = 3
		b = 3
		n = 100
		m.ellipse(x0,y0,a,b,n,phi=np.random.randn()*360)
	return m

# def circle_test(m)
a12_floats = ['5904475','5904474','5904476','5904469','5904477','5904467','5904478','5904471','5904472','5904468','5904397','5904473']
df = pd.read_pickle(soccom_proj_settings.interpolated_drifter_file)
df = df[df.Cruise.isin(df[df.Lat<-50].Cruise.unique().tolist()+a12_floats)]
df['Lon180'] = df.Lon
df.loc[df.Lon180.values>180,['Lon180']] = df[df.Lon180>180].Lon180-360
df = df.dropna(subset = ['Lat','Lon'])
plt.figure(figsize=(7,6))
m = SBasemap(region='lynne',dataframe=df,date=datetime.date(2015,2,20),resolution='h')
m.drifter_plot(lineplot=True,color='y',drift_type=1)
m.drifter_plot(lineplot=True,color='r',drift_type=8)
m.linespace(line_space=10)
m.bathy(color=plt.cm.Greys)
m.ice_map(line=3,col='orange')
df = df[df.Cruise.isin(a12_floats)]
print df.Cruise.unique()
markersz=2
for floater in df.Cruise.unique():
	xcoord = df[df.Cruise==floater].Lon.values.tolist()
	ycoord = df[df.Cruise==floater].Lat.values.tolist()
	xpos,ypos = m(xcoord,ycoord)
	m.plot(xpos,ypos,'bo',markersize=markersz,markeredgewidth=0.1,zorder=10)
	m.plot([xpos[i] for i in [0,-1]],[ypos[i] for i in [0,-1]],'bo',markersize=markersz+4,markeredgewidth=0.1,zorder=10)
m.a12_stations()
m.polarstern_stations()
m.orsi_fronts()
plt.savefig('float_map',format='eps')
plt.show()

# df.Date = pd.to_datetime(df.Date)
# for date in [[2013,1,1]]:
	# for frame in [df[(df.Lon180<-170)|(df.Lon180>170)],df[(df.Lon<10)|(df.Lon>350)]]:
	# 	plt.figure()
	# 	m = SBasemap(dataframe=frame,date=date,timedelta = 100)
	# 	# m.buoyancy_plot()
	# 	# m.co2flux_plot()
	# 	# m = ellipse_test(m)
	# 	# m.ice_map()
	# 	# m.streamline_plot()
	# 	plt.title('Dataframe tests')
	# 	m.plot_features()
	# 	# m.eke()
	# 	m.bathy()
	# 	m.drifter_plot(drift_type='Argo',markersz=20) 
	# 	# m.etotal()
	# 	# m.eaverage()
	# for region in ['polar','atlantic','pacific']:
	# 	print region
	# 	plt.figure()
	# 	m = SBasemap(region=region,date=date,timedelta = -100)
		# ellipse_test(m)
		# for lon in np.linspace(0,360,12):
  #  			m.equi(lon, -60,500)
		# m.plot_features(fontsz=8)
		# m.buoyancy_plot()
		# m.co2flux_plot()
		# m = ellipse_test(m)
		# m.bathy()
		# m.drifter_plot()
		# m.annotate()
		# m.SOCCOM()
		# m.eke()
		# m.streamline_plot()
		# m.ice_map()
		# m.plot_features()
		# m.etotal()
		# m.eaverage()
	# 	plt.title(region)
	# plt.show()