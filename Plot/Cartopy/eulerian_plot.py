import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib.pyplot as plt
import numpy as np
from cartopy import geodesic



def ellipse(geod,lon, lat, b, a, n_samples=360,phi=0):
    """
    Return the coordinates of a geodetic ellipse of a given
    x radius of a and y radius of b around a lon/lat point.

    Radius is in meters in the geodetic's coordinate system.

    """

    tsid = a*b
    radius = []
    phi_r = phi* np.pi / 180
    for az in (np.linspace(360, 0, n_samples)* np.pi / 180):
        A = a * np.sin(az)*np.cos(phi_r)+b*np.cos(az)*np.sin(phi_r)
        B = b * np.cos(az)*np.cos(phi_r)-b*np.sin(az)*np.sin(phi_r)
        r = tsid / (B**2. + A**2.)**0.5
        radius.append(r)

    lons, lats, back_azim = geod.fwd(np.repeat(lon, n_samples),
                                     np.repeat(lat, n_samples),
                                     np.linspace(360, 0, n_samples),
                                     np.array(radius),
                                     radians=False,
                                     )
    return (lons,lats)


class BaseCartopy():
    def __init__(self,lat_grid,lon_grid,ax=False,projection=ccrs.PlateCarree()):
        assert max(lat_grid)<=90
        assert min(lat_grid)>=-90
        assert max(lon_grid)<=180
        assert min(lon_grid)>=-180

        self.lat_grid = lat_grid
        self.lon_grid = lon_grid

        if not ax:
            fig = plt.figure()
            self.ax = fig.add_subplot(1, 1, 1, projection=projection)
        else:
            self.ax = ax

    def meshgrid_xx_yy(self):
        return np.meshgrid(self.lon_grid,self.lat_grid)

    def finish_map(self):
        self.ax.add_feature(cfeature.LAND,zorder=10)
        self.ax.add_feature(cfeature.COASTLINE,zorder=10)
        self.ax.set_aspect('auto')
        gl = self.ax.gridlines(draw_labels=True)
        gl.xlabels_top = False
        gl.ylabels_right = False

    def get_map(self):
        XX,YY = self.meshgrid_xx_yy()
        return (XX,YY,self.ax)



class GlobalCartopy(BaseCartopy):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)      
        print('I am plotting global region')
        llcrnrlon=-180.
        llcrnrlat=-80.
        urcrnrlon=180.
        urcrnrlat=80
        self.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat], crs=ccrs.PlateCarree())
        self.finish_map()

class HypernavCartopy(BaseCartopy):
    def __init__(self,nc,float_pos_dict,*args,pad=1,**kwargs):
        super().__init__(*args,**kwargs)     
        try:
            urlat = nc['lat'][:].max()
            lllat = nc['lat'][:].min()
            urlon = nc['lon'][:].max()
            lllon = nc['lon'][:].min()
        except TypeError:
            urlat = max(nc.lats)
            lllat = min(nc.lats)
            urlon = max(nc.lons)
            lllon = min(nc.lons)        
        llcrnrlon=(lllon-pad)
        llcrnrlat=(lllat-pad)
        urcrnrlon=(urlon+pad)
        urcrnrlat=(urlat+pad)
        self.ax.scatter(float_pos_dict['lon'],float_pos_dict['lat'],c='pink',linewidths=5,marker='x',s=80,zorder=10)
        self.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat], crs=ccrs.PlateCarree())
        self.finish_map()

class PointCartopy(BaseCartopy):
    def __init__(self,geopy_point,*args,pad=1,**kwargs):
        super().__init__(*args,**kwargs)          
        llcrnrlon=(geopy_point.longitude-pad)
        llcrnrlat=(geopy_point.latitude-pad)
        urcrnrlon=(geopy_point.longitude+pad)
        urcrnrlat=(geopy_point.latitude+pad)
        self.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat], crs=ccrs.PlateCarree())
        self.finish_map()


