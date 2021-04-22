def convert_lon_to_180(masked_array):
	masked_array[masked_array>180]=masked_array[masked_array>180]-360
	masked_array[masked_array<-180]=masked_array[masked_array<-180]+360
	return masked_array
