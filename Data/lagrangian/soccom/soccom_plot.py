import sys,os
sys.path.append(os.path.abspath("../"))
import soccom_proj_settings
import pandas as pd
import numpy as np
from scipy.interpolate import griddata
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import cmocean
import matplotlib.colors as colors
import matplotlib.cm as cmx
import subprocess


variable_list = {'Temperature':(cmocean.cm.thermal,'Temperature ($^\circ$C)'),
'Salinity':(cmocean.cm.haline,'Salinity (PSS)'),'Oxygen':(cmocean.cm.oxy,'Oxygen (micromole/liter)'),
'OxygenSat':(cmocean.cm.oxy_r,'Oxygen Saturation (%)'),'Nitrate':(cmocean.cm.matter,'Nitrate (micromole/kg)'),
'pH25C':(cmocean.cm.curl,'PH at 25$^\circ$C'),'TALK_MLR':(cmocean.cm.algae,'Alkalinity (micromole/kg)'),
'pHinsitu':(cmocean.cm.balance,'Insitu pH'),'DIC_MLR':(cmocean.cm.delta,'Dissolved Inorganic Carbon (micromole/kg)')} 

def property_plot(df_):
	date_list = (df_.Date.drop_duplicates()-df_.Date.min()).dt.days.values.tolist()
	plt.gca().xaxis.set_major_formatter(mdates.DateFormatter('%m/%Y'))
	if max(date_list)>360:
		plt.gca().xaxis.set_major_locator(mdates.MonthLocator(interval=4))
	else:
		plt.gca().xaxis.set_major_locator(mdates.MonthLocator())
	x = df_.Date.unique().shape[0]
	y = df_[variable].values.shape[0]/x
	tokenA = pd.DataFrame(df_[variable].values.reshape(x,y).T,columns=df_.Date.unique(),index=df_[df_.Date==df_.Date.min()].Pressure.values)
	token = tokenA.reset_index().dropna(subset=['index']).set_index('index') #this is because some pressure values in the record are nan.
	token[token.index<1600]	
	pressure_list = token.index.values.tolist()
	X,Y=np.meshgrid(df_.Date.unique(),pressure_list)
	plt.pcolormesh(X,Y,token.interpolate(method='akima',limit=15, limit_direction='both').values)#,shading='gouraud',cmap=variable_list[variable][0],alpha=0.9) 
	plt.ylim([0,1600])
	cbar = plt.colorbar(label=variable_list[variable][1])

	CS = plt.contour(X,Y,token.interpolate().values,linewidth=3,colors='k',alpha=.9)
	# plt.clabel(CS, inline=1, fontsize=7,fmt='%1.1f')
	plt.gcf().autofmt_xdate()
	plt.gca().invert_yaxis()
	plt.ylabel('Pressure (db)')
	plt.xlabel('Date')
	plt.title('Float '+str(cruise))
	plt.savefig(str(cruise)+'_'+variable+'_section')
	plt.close()


def line_plot(df_):
	NCURVES = len(df_.Date.unique())
	values = range(NCURVES)
	jet = cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	scalarMap.set_array([])
	for n,date in enumerate(df_.Date.unique()):
		colorVal = scalarMap.to_rgba(values[n])
		df_holder = df_[df.Date==date].dropna(subset=[variable])
		try:
			df_holder.set_index(variable)['Pressure'].plot(color=colorVal)
		except TypeError:	#this is the case when the entire profile has been flagged as eroneous and needs to be passed
			pass
	plt.ylim([0,1600])
	var_max = df_[variable].max()
	var_min = df_[variable].min()
	var_over = (var_max-var_min)*0.05
	plt.xlim([(var_min-var_over),(var_max+var_over)])
	plt.colorbar(scalarMap,label='Profile Number')
	plt.gca().invert_yaxis()
	plt.title('Float '+str(cruise))
	plt.ylabel('Pressure (db)')
	plt.xlabel(variable_list[variable][1])
	plt.savefig(str(cruise)+'_'+variable+'_lineplot')
	plt.close()

def ts_plot(df_):
	NCURVES = len(df_.Date.unique())
	values = range(NCURVES)
	jet = cm = plt.get_cmap('jet') 
	cNorm  = colors.Normalize(vmin=0, vmax=values[-1])
	scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
	scalarMap.set_array([])
	for n,date in enumerate(df_.Date.unique()):
		colorVal = scalarMap.to_rgba(values[n])
		df_holder = df_[df.Date==date].dropna(subset=['Salinity'])
		df_holder.set_index('Salinity')['Temperature'].plot(color=colorVal)
	t_max = df_['Temperature'].max()
	t_min = df_['Temperature'].min()
	t_over = (t_max-t_min)*0.05
	s_max = df_['Salinity'].max()
	s_min = df_['Salinity'].min()
	s_over = (s_max-s_min)*0.05
	plt.ylim([(t_min-t_over),(t_max+t_over)])
	plt.xlim([(s_min-s_over),(s_max+s_over)])
	plt.colorbar(scalarMap,label='Profile Number')
	plt.title('Float '+str(cruise))
	plt.ylabel('Tempeature ($^\circ$C)')
	plt.xlabel('Salinity (PSS)')
	plt.savefig(str(cruise)+'_tsplot')	
	plt.close()

df = pd.read_pickle(soccom_proj_settings.soccom_drifter_file)
df = df[df.Cruise.isin(['5904469','5904472'])]
for cruise in df.Cruise.unique():
	print cruise
	for variable in variable_list:
		print variable
		df_token = df[(df.Cruise==cruise)][['Pressure','Date','PosQC']+[variable]]
		property_plot(df_token)
		line_plot(df_token)
	ts_plot(df[df.Cruise==cruise])
subprocess.call(['./convert_eps.sh'])