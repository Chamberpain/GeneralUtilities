import numpy as np
import geopy

degree_dist = 111.32 # km/degree
seconds_per_day = 24*3600. 
seconds_per_sideral_day = 23*3600 + 56*60 + 4.1
r_earth = 6378 # km
omega = np.pi/seconds_per_sideral_day

def calculate_f(point):
	assert issubclass(geopy.Point,point.__class__)
	return 2 * omega * np.sin(np.deg2rad(point.latitude))