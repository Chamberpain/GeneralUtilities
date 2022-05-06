import numpy as np
import geopy

degree_dist = 111.32 # km/degree
seconds_per_day = 24*3600. 
seconds_per_sideral_day = 23*3600 + 56*60 + 4.1
r_earth = 6378 # km
omega = (2*np.pi)/seconds_per_sideral_day
g = 9.81 # m/s^2


def calculate_f(point):
	assert issubclass(geopy.Point,point.__class__)
	return 2 * omega * np.sin(np.deg2rad(point.latitude))

def calculate_barotropic_rossby_def_rad(point,D=4): #D is in km
	return np.sqrt(g*D)/abs(calculate_f(point))