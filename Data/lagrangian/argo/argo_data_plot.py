from argo_read import all_argo_holder,argo_reader,prof_class,meta_class,traj_class,tech_class
import numpy as np 
import matplotlib.pyplot as plt
import copy 

def plot_gridded_speed_and_velocity_mean_and_variance():
	all_argo = all_argo_holder()	
	import sys
	import os
	sys.path.append(os.path.join(os.getcwd(),'../../'))
	from transition_matrix.makeplots.plot_utils import basemap_setup
	speed,dx,dy = all_argo.get_full_speed_list()
	data = np.array(speed)
	data[data==None]=np.nan
	data = data.astype(float)
	number = all_argo.data_to_gridded(~np.isnan(data),np.count_nonzero)
	plt.figure()
	XX,YY,m = basemap_setup(all_argo.lat_bins,all_argo.lon_bins,'Argo')
	number = np.ma.masked_less(number,20)
	m.pcolormesh(XX,YY,number)
	plt.colorbar()
	speed_mean = all_argo.data_to_gridded(speed,np.nanmean)
	plt.figure()
	XX,YY,m = basemap_setup(all_argo.lat_bins,all_argo.lon_bins,'Argo')
	m.pcolormesh(XX,YY,np.ma.masked_array(speed_mean,number.mask))
	plt.colorbar()
	dx_mean = all_argo.data_to_gridded(dx,np.nanmean)
	dy_mean = all_argo.data_to_gridded(dy,np.nanmean)
	plt.figure()
	XX,YY,m = basemap_setup(all_argo.lat_bins,all_argo.lon_bins,'Argo')
	m.quiver(XX,YY,np.ma.masked_array(dx_mean,number.mask),dy_mean)
	speed_var = all_argo.data_to_gridded(speed,np.nanvar)
	plt.figure()
	XX,YY,m = basemap_setup(all_argo.lat_bins,all_argo.lon_bins,'Argo')
	m.pcolormesh(XX,YY,np.ma.masked_array(speed_var,number.mask))
	plt.show()


def plot_gridded_speed_and_velocity_mean_and_variance():
	import sys
	import os
	sys.path.append(os.path.join(os.getcwd(),'../../'))
	from transition_matrix.makeplots.plot_utils import basemap_setup

	all_argo = all_argo_holder()
	argos_list,gps_list = all_argo.get_pos_list()	
	argos_class_list = [all_argo.class_list.get(key) for key in argos_list]
	gps_class_list = [all_argo.class_list.get(key) for key in gps_list]

	def get_speed_from_class_list(name_list):
		all_argo.class_list = name_list
		speed,dx,dy = all_argo.get_full_speed_list()
		data = np.array(speed)
		data[data==None]=np.nan
		data = data.astype(float)
		number = all_argo.data_to_gridded(~np.isnan(data),np.count_nonzero)
		speed_mean = all_argo.data_to_gridded(speed,np.nanmean)
		dx_mean = all_argo.data_to_gridded(dx,np.nanmean)
		dy_mean = all_argo.data_to_gridded(dy,np.nanmean)
		return (number,speed_mean,dx_mean,dy_mean)

	argos_number,argos_speed,argos_dx,argos_dy = get_speed_from_class_list(dict(zip(argos_list,argos_class_list)))
	gps_number,gps_speed,gps_dx,gps_dy = get_speed_from_class_list(dict(zip(gps_list,gps_class_list)))

	def basemap_pcolormesh(data):
		plt.figure()
		XX,YY,m = basemap_setup(all_argo.lat_bins,all_argo.lon_bins,'Argo')
		m.pcolormesh(XX,YY,data)
		plt.colorbar()

	def basemap_quiver(data_x,data_y):
		plt.figure()
		XX,YY,m = basemap_setup(all_argo.lat_bins,all_argo.lon_bins,'Argo')
		m.quiver(XX,YY,data_x,data_y)

	basemap_pcolormesh(argos_number)
	plt.title('Argos Number')
	basemap_pcolormesh(gps_number)
	plt.title('GPS Number')
	basemap_pcolormesh(argos_speed-gps_speed)
	plt.title('Argos Mean Speed - GPS Mean Speed')
	basemap_quiver(argos_dx-gps_dx,argos_dy-gps_dy)
	plt.title('Argos Mean Vel - GPS Mean Vel')
	plt.show()