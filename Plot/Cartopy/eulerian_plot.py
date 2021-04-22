# import sys,os
# sys.path.append(os.path.abspath("./data"))
# from mpl_toolkits.basemap import Basemap,shiftgrid
# import numpy as np
# import pandas as pd
# import pyproj
# from matplotlib.patches import Polygon
import datetime
import osr # for some bizarre reasion, this needs to be imported relatively early or program will crash
import pyproj
import shapefile
import scipy.io
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy
import cartopy.util
import cartopy.crs as ccrs
from cartopy.io.shapereader import Reader
from cartopy.feature import ShapelyFeature
from oceans.colormaps import cm as ocm
import oceans
import os
import ftplib
import zipfile
import glob
import numpy as np 
import itertools
import regionmask
from mpl_toolkits.basemap import shiftgrid
# import scipy.io
# import pickle
# from matplotlib.colors import LogNorm

region_dict = {'pacific':((160,290,-75,-30),180),'atlantic':((-60,20,-75,-30),0),
                'weddell':((-60,20,-77,-60),0),'ross':((170,250,-77,-60),180)}
map_features=[['Weddell Sea',-45,-70],['Ross Sea',195,-70], ['Ross Polynya',180, -80],['Maud Rise Polynya',0,-70],['Cosmonaut Polynya',45,-70],['Darnley Polynya',75,-70],['Weddel Polyna',-60,-70], ['Mertz Polyna',150,-70]]
degree_dist = 111
# datadir = os.sep.join([os.path.dirname(os.path.abspath(__file__)), 'data'])



#these things should be moved to a data directory
# sound_speed = 1.467 #km/s
# slow = 1/sound_speed


# SOSO_weddel_coord = {'W1a':(-63.9538,-00.0035),'W1b':(-63.953,-00.006),'W1c':(-63.953,-00.028),'W1d':(-63.967,-00.052)
# ,'W1e':(-63.994,-00.044),'W1f':(-63.994,-00.044),'W2a':(-66.509,-00.0338),'W2b':(-66.511,-00.032)
# ,'W2c':(-66.511,-00.032),'W2d':(-66.511,-00.032),'W2e':(-66.511,-00.025),'W2f':(-66.515,-00.011)
# ,'W3a':(-64.4915,9.8255),'W3b':(-64.4915,9.8255),'W4a':(-66.618,-27.105),'W4b':(-66.615,-27.118)
# ,'W4c':(-66.612,-27.122),'W4d':(-66.608,-27.121),'W5a':(-65.619,-36.392),'W5b':(-65.617,-36.421)
# ,'W5c':(-65.621,-36.422),'W6a':(-63.703,-50.870),'W6b':(-63.712,-50.842),'W6c':(-63.718,-50.832)
# ,'W6d':(-63.702,-50.827),'W7a':(-69.0112,-34.0025),'W8a':(-65.5752,-37.1221),'W8b':(-65.5752,-37.1221)
# ,'W9a':(-69.061,-17.4315),'W9b':(-69.059,-17.384),'W9c':(-69.058,-17.389),'W10a':(-64.389,-45.873)
# ,'W10b':(-64.398,-45.866),'W10c':(-64.382,-45.869),'W11a':(-68.995,-6.945),'W11b1':(-69.005,-6.982)
# ,'W11b2':(-69.005,-6.982),'W11c':(-69.006,-6.983),'W12a':(-65.968,-12.252),'W13a':(-70.893,-28.891)
# ,'W14a':(-68.483,-44.111)}
# SOSO_drake_coord = {'US1-1':(-59.914,-109.9105),'US1-2':(-55.996,-109.9102),'US1-3':(-58.1010,-97.96)
# ,'US1-4':(-59.9655,-86.1536),'US1-5':(-55.9630,-85.9863),'US1-6':(-58.3343,-74.4142),'US2-1':(-58.1690,-102.169)
# ,'US2-2':(-59.939,-78.003),'UK1-1':(-59.977,-65.990),'UK1-2':(-58.014,-61.886),'UK1-3':(-56.023,-57.784)
# ,'UK1-4':(-58.051,-53.645)}


def download_seaice(date,local_root):
    month_list = ['01_Jan/','02_Feb/','03_Mar/','04_Apr/','05_May/','06_Jun/','07_Jul/','08_Aug/','09_Sep/','10_Oct/','11_Nov/','12_Dec/']
    month_string = month_list[date.month-1]

    # first need to download the zipped file into the directory
    file_path = 'extent_S_'+str(date.year)+'{:02}_polygon_v2.1.zip'.format(date.month)
    ftp = ftplib.FTP('sidads.colorado.edu')
    ftp.login()
    ftp.cwd('/DATASETS/NOAA/G02135/south/monthly/shapefiles/shp_extent/'+month_string)
    with open(os.path.join(local_root,file_path), 'w') as fobj:
        ftp.retrbinary('RETR %s' % file_path, fobj.write)
    ftp.quit()

    #second need to unzip the file
    zip_ref = zipfile.ZipFile(os.path.join(local_root,file_path), 'r')
    zip_ref.extractall(local_root)
    zip_ref.close()
    os.remove(os.path.join(local_root,file_path)) #delete the original zip file

    #third need to rename the zipped file because the prj4 libraries seem to have a problem with . characters
    for filename in glob.glob(os.path.join(local_root,file_path.split('.')[0]+'*')):
        print 'renameing %s' % filename
        split = filename.split('.')
        os.rename(filename, split[0]+'-'+split[1]+'.'+split[2])




class SBasemap():
    def __init__(self, dataframe=None,region=None,padding=2,row=1,column=1,plot_num=1,gridlines=True): #,date=None,resolution='l',timedelta=0,cruise=None,soso_df=None):


####################  dataframe setup section ##################
        self.region=region
        if region is None:
            raise NameError('Need to specify region')

        if (dataframe is None)&(region is 'float'):
            raise NameError('Dataframe must be specified if region is float')
        elif (region is 'float'):
            lllat=dataframe.Lat.min()-padding
            if lllat<-90:
                lllat = -90
            urlat=dataframe.Lat.max()+padding
            if urlat>90:
                urlat = 90

            if dataframe.Lon.apply(oceans.wrap_lon180).diff().max()[0]>dataframe.Lon.apply(oceans.wrap_lon360).diff().max()[0]:
                print 'going to 360 coord'
                print padding
                lllon=dataframe.Lon.min()-padding
                urlon=dataframe.Lon.max()+padding
                if urlon > 360: 
                    urlon = 360
                if lllon<0:
                    lllon= 0
                self.lon_mid = 180
            else: 
                print 'going to 180 coord'
                lllon=dataframe.Lon.apply(oceans.wrap_lon180).min()[0]-padding
                urlon=dataframe.Lon.apply(oceans.wrap_lon180).max()[0]+padding
                if urlon > 180: 
                    urlon = 180
                if lllon<-180:
                    lllon=-180
                self.lon_mid=0
            self.ax = plt.subplot(row,column,plot_num,projection=cartopy.crs.Miller(central_longitude=self.lon_mid))

            self.ax.set_extent([lllon, urlon, lllat, urlat])        
            i = list(itertools.product([lllon,urlon],[lllat,urlat]))
            i[2], i[3] = i[3], i[2]  
            self.mask = regionmask.Regions_cls(region,[0],[region],[region],[i])


        elif region == 'polar':
            self.lon_mid=-60
            self.ax = plt.subplot(row,column,plot_num,projection=cartopy.crs.SouthPolarStereo(central_longitude=self.lon_mid))
            self.ax.set_extent([-180-self.lon_mid, 180-self.lon_mid, -90, -45], ccrs.PlateCarree())
            theta = np.linspace(0, 2*np.pi, 100)
            center, radius = [0.5, 0.5], 0.5
            verts = np.vstack([np.sin(theta), np.cos(theta)]).T
            circle = mpath.Path(verts * radius + center)
            self.ax.set_boundary(circle,transform=self.ax.transAxes)

        else:
            (lllon,urlon,lllat,urlat),self.lon_mid= region_dict[region]
            self.ax = plt.subplot(row,column,plot_num,projection=cartopy.crs.Miller(central_longitude=self.lon_mid))
            self.ax.set_extent([lllon, urlon, lllat, urlat])   

        if gridlines:
            self.ax.gridlines()

        self.ax.add_feature(cartopy.feature.LAND)
        self.ax.add_feature(cartopy.feature.COASTLINE)
#### set the region mask
        if region in region_dict:
            i = list(itertools.product(region_dict[region][0][:2],region_dict[region][0][2:]))
            i[2], i[3] = i[3], i[2]  
            self.mask = regionmask.Regions_cls(region,[0],[region],[region],[i])

    def drifter_plot(self,dataframe,markersz=2,lineplot=False,color='b',marker='o',label=None):
        print self.ax.get_extent()
        if self.region != 'polar':
            lllon,urlon,dummy,dummy = self.ax.get_extent()
            lllon += self.lon_mid
            urlon += self.lon_mid
            print lllon
            print urlon
            if self.lon_mid == 180:
                dataframe = dataframe[(dataframe.Lon>lllon)&(dataframe.Lon<urlon)]
            elif self.lon_mid == 0:
                dataframe = dataframe[(dataframe.Lon.apply(oceans.wrap_lon180)>lllon)&(dataframe.Lon.apply(oceans.wrap_lon180)<urlon)]
        if dataframe.empty:
            return
        if lineplot is True:
            for name in dataframe.Cruise.unique():
                df_ = dataframe[dataframe.Cruise==name].drop_duplicates(subset=['Lat','Lon'])

                # if (df_.Lon.diff().abs>60)&(df_.Lon.diff().abs<350):
                #     index_list = df_[(df_.Lon.diff().abs()>60)&(df_.Lon.diff().abs()<350)].index.values.tolist() 
                #     frame_list = [df_[df_.index.isin(range(i,j))] for i, j in zip([0]+index_list, index_list+[df_.index.max()])]
                #     for frame in frame_list:
                #         xcoord = frame.Lon.values
                #         ycoord = frame.Lat.values
                #         self.ax.plot(xcoord,ycoord,linewidth=0.5,alpha = 0.6,transform=ccrs.PlateCarree())
                # else:
                xcoord = df_.Lon.values
                ycoord = df_.Lat.values
                self.ax.plot(xcoord,ycoord,linewidth=0.5,alpha = 0.6,transform=ccrs.PlateCarree())
        else:            
                xcoord = list(dataframe.Lon.values)
                print xcoord
                ycoord = list(dataframe.Lat.values)
                self.ax.scatter(xcoord,ycoord,color=color,marker=marker,s=markersz,transform=ccrs.PlateCarree())             

##########################        Utililties for processing data ###############################




    # def parse_dataset(self,data,XC,YC):
    # #this function parses the data we wish to plot to our specific longitude and latitude extent
    #     if self.region is 'polar':
    #         return (data,XC)
    #     else:
    #         lllon,urlon,lllat,urlat = self.ax.get_extent()
    #         print lllon
    #         print urlon
    #         print urlat
    #         print lllat
    #         if self.lon_mid == 0
    #             lllon_index = min(range(len(XC)), key=lambda i: abs(XC[i]-lllon))
    #             urlon_index = min(range(len(XC)), key=lambda i: abs(XC[i]-urlon))
    #             lllat_index = min(range(len(YC)), key=lambda i: abs(YC[i]-lllat))
    #             urlat_index = min(range(len(YC)), key=lambda i: abs(YC[i]-urlat))
    #             XC = XC[lllon_index:urlon_index]
    #             YC = YC[lllat_index:urlat_index]
    #             data = data[lllat_index:ur_lat_index,lllon_index:urlon_index]
    #         elif self.lon_mid == 180
    #             XC180 = oceans.wrap180(XC)
    #             lllon_index = min(range(len(XC)), key=lambda i: abs(XC[i]-lllon))
    #             urlon_index = min(range(len(XC)), key=lambda i: abs(XC[i]-urlon))
    #             lllat_index = min(range(len(YC)), key=lambda i: abs(YC[i]-lllat))
    #             urlat_index = min(range(len(YC)), key=lambda i: abs(YC[i]-urlat))
    #             XC = XC[lllon_index:urlon_index]
    #             YC = YC[lllat_index:urlat_index]
    #             data = data[lllat_index:ur_lat_index,lllon_index:urlon_index]

    #         return (data,XC)


###########################.    routines for plotting data ##################################


    def ice_plot(self,date,color='b'):
        local_root = '/Users/paulchamberlain/Data/SeaIce/'
        if type(date) is datetime.date:
            file = 'extent_S_'+str(date.year)+'{:02}_polygon_v2-1.prj'.format(date.month)
            local_path = os.path.join(local_root,file)
        else:
            print 'date is not in the datetime format'
            raise NameError('date must be in datetime format')
        try:
            prj_file = open(os.path.join(local_root,file),'r')
        except:
            print 'there was an error opening the proj4 file, downloading file'
            download_seaice(date,local_root)
            prj_file = open(os.path.join(local_root,file),'r')
        prj_txt = prj_file.read()
        srs = osr.SpatialReference()
        srs.ImportFromESRI([prj_txt])

        shpProj = pyproj.Proj(srs.ExportToProj4())
        sf = shapefile.Reader(os.path.join(local_root,file))
        for n,shape in enumerate(sf.shapes()):
            print n
            shpX, shpY = zip(*shape.points)
            Lon,Lat = shpProj(shpX,shpY,inverse=True)
            index_list = np.where((abs((np.diff(Lon)))>15)|(abs((np.diff(Lat)))>10))[0]
# This is a total hack. The downloaded data is not parsed correctly and is drawn as one continuous line. 
# This detects large changes in longitude and parses the polygons appropriately.
            index_list = (index_list+1).tolist()
            Lon = [list(Lon)[i:j] for i, j in zip([0]+index_list, index_list+[None])]
            Lat = [list(Lat)[i:j] for i, j in zip([0]+index_list, index_list+[None])]
            for k,(x,y) in enumerate(zip(Lon,Lat)):
                x[-1]=x[-1]+0.00001  
                y[-1]=y[-1]+0.00001
#catopy seems to be really stupid and if 2 lines are plotted at the same point it breaks. This will hopefully be fixed. 
#the fact that the last point in these polygons is repeated could be used to break up the polygons so they are plotted properly
                self.ax.plot(x,y,linewidth=2,color=color,transform=ccrs.Geodetic())             


    # def toa_to_dist(self,toa):
    #     dist = toa/slow
    #     return dist

 
    # def annotate(self):
    #     df_annotate = df.drop_duplicates(subset=['Cruise'])
    #     namelist = df_annotate.Cruise.values
    #     xcoord = list(df_annotate.Lon.values)
    #     ycoord = list(df_annotate.Lat.values)
    #     xpos,ypos = self(xcoord,ycoord)
    #     x2 = 10*np.sin(np.linspace(-np.pi/2,np.pi/2,len(namelist)))
    #     y2 = 10*np.cos(np.linspace(-np.pi/2,np.pi/2,len(namelist)))
    #     for k in range(len(namelist)):
    #         plt.annotate(namelist[k], xy=(xpos[k], ypos[k]),  xycoords='data',xytext=(x2[k], y2[k]), textcoords='offset points',color='midnightblue',fontsize=10,fontweight='bold')


    # def shoot(self, lon, lat, azimuth, maxdist=None):
    #     """Shooter Function
    #     Original javascript on http://williams.best.vwh.net/gccalc.htm
    #     Translated to python by Thomas Lecocq
    #     """
    #     glat1 = lat*np.pi / 180.
    #     glon1 = lon*np.pi / 180.
    #     s = maxdist / 1.852
    #     faz = azimuth * np.pi / 180.
    #     EPS= 0.00000000005
    #     if ((np.abs(np.cos(glat1))<EPS) and not (np.abs(np.sin(faz))<EPS)):
    #         alert("Only N-S courses are meaningful, starting at a pole!")
    #     a=6378.13/1.852
    #     f=1/298.257223563
    #     r = 1 - f
    #     tu = r * np.tan(glat1)
    #     sf = np.sin(faz)
    #     cf = np.cos(faz)
    #     if (cf==0):
    #         b=0.
    #     else:
    #         b=2. * np.arctan2 (tu, cf)
    #     cu = 1. / np.sqrt(1 + tu * tu)
    #     su = tu * cu
    #     sa = cu * sf
    #     c2a = 1 - sa * sa
    #     x = 1. + np.sqrt(1. + c2a * (1. / (r * r) - 1.))
    #     x = (x - 2.) / x
    #     c = 1. - x
    #     c = (x * x / 4. + 1.) / c
    #     d = (0.375 * x * x - 1.) * x
    #     tu = s / (r * a * c)
    #     y = tu
    #     c = y + 1
    #     while (np.abs (y - c) > EPS):
    #         sy = np.sin(y)
    #         cy = np.cos(y)
    #         cz = np.cos(b + y)
    #         e = 2. * cz * cz - 1.
    #         c = y
    #         x = e * cy
    #         y = e + e - 1.
    #         y = (((sy * sy * 4. - 3.) * y * cz * d / 6. + x) *
    #               d / 4. - cz) * sy * d + tu
    #     b = cu * cy * cf - su * sy
    #     c = r * np.sqrt(sa * sa + b * b)
    #     d = su * cy + cu * sy * cf
    #     glat2 = (np.arctan2(d, c) + np.pi) % (2*np.pi) - np.pi
    #     c = cu * cy - su * sy * cf
    #     x = np.arctan2(sy * sf, c)
    #     c = ((-3. * c2a + 4.) * f + 4.) * c2a * f / 16.
    #     d = ((e * cy * c + cz) * sy * c + y) * sa
    #     glon2 = ((glon1 + x - (1. - c) * d * f + np.pi) % (2*np.pi)) - np.pi    
    #     baz = (np.arctan2(sa, b) + np.pi) % (2 * np.pi)
    #     glon2 *= 180./np.pi
    #     glat2 *= 180./np.pi
    #     baz *= 180./np.pi
     
    #     return (glon2, glat2, baz)
     
    # def equi(self,centerlon, centerlat, radius, *args, **kwargs):
    #     glon1 = centerlon
    #     glat1 = centerlat
    #     X = []
    #     Y = []
    #     for azimuth in range(-180, 180):
    #         glon2, glat2, baz = self.shoot(glon1, glat1, azimuth, radius)
    #         X.append(glon2)
    #         Y.append(glat2)
    #     X.append(X[0])
    #     Y.append(Y[0])
    #     #~ m.plot(X,Y,**kwargs) #Should work, but doesn't...
    #     if self.lon_mid == 180:
    #         X = np.array(X)+180

    #     X,Y = self(X,Y)
    #     self.plot(X,Y,**kwargs)


    # def ellipse(self, x0, y0, a, b, n, phi=0,ax=None, **kwargs):
    #     """
    #     Draws a polygon centered at ``x0, y0``. The polygon approximates an
    #     ellipse on the surface of the Earth with semi-major-axis ``a`` and 
    #     semi-minor axis ``b`` degrees longitude and latitude, made up of 
    #     ``n`` vertices.

    #     For a description of the properties of ellipsis, please refer to [1].

    #     The polygon is based upon code written do plot Tissot's indicatrix
    #     found on the matplotlib mailing list at [2].

    #     Extra keyword ``ax`` can be used to override the default axis instance.

    #     Other \**kwargs passed on to matplotlib.patches.Polygon

    #     RETURNS
    #         poly : a maptplotlib.patches.Polygon object.

    #     REFERENCES
    #         [1] : http://en.wikipedia.org/wiki/Ellipse


    #     """
    #     # if (self.lon_mid==0)&(x0>180):
    #     #     x0 = x0-360
    #     # if (self.lon_mid==180)&(x0>180):
    #     #     x0 = x0-720
    #         # print 'I have subtracted'
    #     print self.lon_mid
    #     a = a/abs(np.cos(np.deg2rad(y0)))
    #     ax = kwargs.pop('ax', None) or self._check_ax()
    #     g = pyproj.Geod(a=self.rmajor, b=self.rminor)
    #     # Gets forward and back azimuths, plus distances between initial
    #     # points (x0, y0)
    #     azf, azb, dist = g.inv([x0, x0], [y0, y0], [x0+a, x0], [y0, y0+b])
    #     tsid = dist[0] * dist[1] # a * b

    #     # Initializes list of segments, calculates \del azimuth, and goes on 
    #     # for every vertex
    #     seg = []
    #     AZ = np.linspace(azf[0], 360. + azf[0], n)
    #     for i, az in enumerate(AZ):
    #         # Skips segments along equator (Geod can't handle equatorial arcs).
    #         # az =+ phi
    #         if np.allclose(0., y0) and (np.allclose(90., az) or
    #             np.allclose(270., az)):
    #             continue
    #         # In polar coordinates, with the origin at the center of the 
    #         # ellipse and with the angular coordinate ``az`` measured from the
    #         # major axis, the ellipse's equation  is [1]:
    #         #
    #         #                           a * b
    #         # r(az) = ------------------------------------------
    #         #         ((b * cos(az))**2 + (a * sin(az))**2)**0.5
    #         #
    #         # Azymuth angle in radial coordinates and corrected for reference
    #         # angle.
    #         azr = 2. * np.pi / 360. * (az + 90.+phi)
    #         A = dist[0] * np.sin(azr)
    #         B = dist[1] * np.cos(azr)
    #         r = tsid / (B**2. + A**2.)**0.5
    #         lon, lat, azb = g.fwd(x0, y0, az, r)
    #         x, y = self(lon, lat)
    #         # Add segment if it is in the map projection region.
    #         if x < 1e20 and y < 1e20:
    #             seg.append((x, y))
    #     segx = np.array(zip(*seg)[0])
    #     if not (all(item >= 0 for item in segx) or all(item < 0 for item in segx)):
    #         return
    #     segy = zip(*seg)[1]
    #     #     segx = [abs(x) for x in segx]
    #     seg = zip(segx,segy)
    #     poly = Polygon(seg, **kwargs)
    #     ax.add_patch(poly)

    #     # Set axes limits to fit map region.
    #     self.set_axes_limits(ax=ax)

    #     return poly


    #         # az +=phi
    #         # # Skips segments along equator (Geod can't handle equatorial arcs).
    #         # if np.allclose(0., y0) and (np.allclose(90., az) or
    #         #     np.allclose(270., az)):
    #         #     continue
    #         # # In polar coordinates, with the origin at the center of the 
    #         # # ellipse and with the angular coordinate ``az`` measured from the
    #         # # major axis, the ellipse's equation  is [1]:
    #         # #                           a * b
    #         # # r(az) = ------------------------------------------
    #         # #         ((b * cos(az))**2 + (a * sin(az))**2)**0.5
    #         # # Azymuth angle in radial coordinates and corrected for reference
    #         # # angle.
    #         # azr = 2. * np.pi / 360. * (az + 90.)
    #         # A = dist[0] * np.sin(azr)
    #         # B = dist[1] * np.cos(azr)
    #         # r = tsid / (B**2. + A**2.)**0.5
    #         # lon, lat, azb = g.fwd(x0, y0, az, r)
    #         # x, y = self(lon, lat)




    def streamline_plot(self,subsample=1,levs=10):
        local_root = '/Users/paulchamberlain/Data/'
        mat = scipy.io.loadmat(os.path.join(local_root,'agva_mean_streamfunction_1000dbar.mat'))
        XC = mat['long'][0][:]
        YC = mat['latg'][0][:]
        streamline=mat['mean_streamfunction_1000dbar'][:,:]
# load in the data
        lon = np.ones(XC.shape)
        lon[:291]=XC[69:]
        lon[291:360]=XC[0:69]
        shifted = np.ones(streamline.shape)
        shifted[:,:291]=streamline[:,69:]
        shifted[:,291:360]=streamline[:,0:69]
        shifted,lon = shiftgrid(1,shifted,lon,start=True)
        shifted,lon = shiftgrid(180+self.lon_mid,shifted,lon,start=False)
# These data do not start at lon = 0, so need to shift them around
        # streamline_mask = np.ma.masked_invalid(shifted)
# Add cyclic point so there are no gaps in the plot
        if self.region!='polar':
            mask = self.mask.mask(lon,YC)
            shifted = np.ma.array(shifted,mask=mask)
        shifted, lon = cartopy.util.add_cyclic_point(shifted,lon)
# need to have this after the mask so that the dimensions line up
        levels = np.linspace(np.nanmin(shifted),np.nanmax(shifted),levs)
# Calculate levels to be plotted
        self.ax.contourf(lon[::subsample],YC[::subsample],shifted[::subsample,::subsample],linewidths=2, label='Gray Streamlines',levels=levels,transform=ccrs.PlateCarree())

    def ssh_plot(self):
        local_root = '/Users/paulchamberlain/Data/SOSE/'
        ssh = np.load(os.path.join(local_root,'sshave.npy'))
        mat = scipy.io.loadmat(os.path.join(local_root,'grid.mat'))
        plot_data = np.ma.masked_equal(ssh,0.)
        XC = mat['XC'][:,0]
        YC = mat['YC'][0,:]
        try:
            plot_data,XC = shiftgrid(180+self.lon_mid,plot_data,XC,start=False)
        except ValueError:
            print 'value error in shiftgrid'
        if self.region!='polar':
            mask = self.mask.mask(XC,YC)
            plot_data = np.ma.array(plot_data,mask=mask)
        lev = np.max([abs(np.nanmin(plot_data)),abs(np.nanmax(plot_data))])
        levels = np.linspace(-lev,lev,20)
        cs = self.ax.contour(XC,YC,plot_data,levels,linewidths=0.5,colors='k',animated=True,transform=ccrs.PlateCarree(),alpha=0.4)
        ca = self.ax.contourf(XC,YC,plot_data,levels,transform=ccrs.PlateCarree(),cmap=plt.cm.RdBu_r,animated=True,alpha=0.4)
        plt.colorbar(ca)
        plt.clabel(cs,inline=1,fontsize=6)


    def bathy(self,color=ocm.cbathy,subsample=10,delta=False,cbar = True,f_over_h=False,fill=True,levs=None):
        local_root = '/Users/paulchamberlain/Data/SOSE/'
        mat = scipy.io.loadmat(os.path.join(local_root,'grid.mat'))
        XC = np.linspace(0,360,6*360)[::subsample]
# we add this because we need a perfectly spaced x grid (SOSE is not) or else it causes problems.
        YC = mat['YC'][0,::subsample]
        plot_data = mat['Depth'][::subsample,::subsample].T/1000.
        plot_data = np.ma.masked_equal(plot_data,0)
# mask out the land
        if delta:
            plot_data = np.sum(np.gradient(plot_data,YC,XC),axis=0)
# calculate the 2d spatial derivative. This will be slightly wrong because the longitude spacing changes with latitude, but we will ignore this for the time being
# this throws an error when 360 is our of the domain, hence the exception handling
        if f_over_h:
            depth_data = np.ma.masked_less_equal(depth_data,100)*1000
# mask depths to below a reasonable value so that we get some resolution in our streamline plots
            dummy,lat_grid = np.meshgrid(XC,YC)
# compute latitude grid 
            omega = 7.2921 * 10**-5
            f = 2*omega*np.sin(np.deg2rad(lat_grid))
#compute grid of f using latitudes
            plot_data = -f/depth_data

            plot_data = np.log(plot_data)
#plot the log so that we can actually see differences in streamlines
        try:
            plot_data,XC = shiftgrid(180+self.lon_mid,plot_data,XC,start=False)
        except ValueError:
            print 'value error in shiftgrid'
# this throws an error when 360 is our of the domain, hence the exception handling
        if self.region!='polar':
            mask = self.mask.mask(XC,YC)
            plot_data = np.ma.array(plot_data,mask=mask)
        if levs:
            levels = np.arange(np.nanmin(plot_data),np.nanmax(plot_data),levs)
        if fill:
            im = self.ax.contourf(XC,YC,plot_data,cmap=color,transform=ccrs.PlateCarree(),alpha = 0.3)
        else:
            im = self.ax.contour(XC,YC,plot_data,cmap=color,levels=levs,transform=ccrs.PlateCarree(),alpha = 0.3)
        if cbar:
            plt.colorbar(im,label='Depth (km)')


    # def buoyancy_plot(self):
    #     file = os.path.join(datadir,'averageBF.npy')
    #     mat = scipy.io.loadmat(os.path.join(datadir,'grid.mat'))
    #     b = np.load(file)
    #     b_mask = np.ma.array(b,mask=((b==0.)|(b==1.)))
    #     b_mask[b_mask<-200]=-200
    #     b_mask[b_mask>200]=200
    #     XC = mat['XC'][:,0]
    #     YC = mat['YC'][0,:]
    #     if self.lon_mid==0:
    #         b_mask, XC = shiftgrid(180, b_mask, XC, start=False)
    #     XC,YC = np.meshgrid(XC,YC)
    #     XX,YY = self(XC,YC)
    #     self.pcolormesh(XX[::self.subsample,::self.subsample],YY[::self.subsample,::self.subsample],b_mask[::self.subsample,::self.subsample],cmap=plt.cm.PuOr,alpha=0.4)
    #     plt.colorbar(label='$W/m^2$',fraction=0.046, pad=0.04)

    # def plot_features(self,fontsz=16):
    #     for feature in map_features: 
    #         XX,YY = self(feature[1],feature[2])
    #         plt.annotate(feature[0], xy=(XX,YY),  xycoords='data',xytext=(XX,YY),color='teal',fontsize=fontsz,fontweight='bold') # ,arrowprops=dict(arrowstyle="fancy", color='b'))

    # def co2flux_plot(self):
    #     df = pd.read_csv(os.path.join(datadir,'sumflux_2006c.txt'), skiprows = 51,sep=r"\s*")
    #     x = np.arange(df.LON.min(),df.LON.max()+5,5)
    #     y = np.arange(df.LAT.min(),df.LAT.max()+4,4)
    #     XC,YC = np.meshgrid(x,y)
    #     CO2 = np.zeros([len(y),len(x)])
    #     di = df.iterrows()
    #     for i in range(len(df)):
    #         row = next(di)[1]
    #         CO2[(row['LON']==XC)&(row['LAT']==YC)] = row['FLUXGMSW06']
    #     if self.lon_mid==0:
    #         CO2, x = shiftgrid(180, CO2, x, start=False)
    #     CO2_mask = np.ma.masked_equal(CO2,0.)
    #     XC,YC = np.meshgrid(x,y)
    #     XX,YY = self(XC,YC)
    #     self.pcolormesh(XX,YY,CO2_mask,cmap=plt.cm.PRGn)
    #     plt.colorbar(label='CO2 Flux $gm C/m^2/yr$',fraction=0.046, pad=0.04)

    # def temp_plot(self,level='bottom',colormap=plt.cm.Greys_r):
    #     if level is 'top':
    #         level=0
    #     elif level is 'bottom':
    #         level = 18
    #     temp = np.load(os.path.join(datadir,'averageT_'+str(level)+'.npy'))
    #     mat = scipy.io.loadmat(os.path.join(datadir,'grid.mat'))
    #     temp_mask = np.ma.masked_equal(temp,0.)
    #     temp_mask[temp_mask>2]=2
    #     XC = mat['XC'][:,0]
    #     YC = mat['YC'][0,:]
    #     if self.lon_mid==0:
    #         temp_mask, XC = shiftgrid(180, temp_mask, XC, start=False)
    #     XC,YC = np.meshgrid(XC,YC)
    #     XX,YY = self(XC,YC)
    #     clevs = np.linspace(np.nanmin(temp_mask),np.nanmax(temp_mask),40)
    #     self.contourf(XX[::self.subsample,::self.subsample],YY[::self.subsample,::self.subsample],temp_mask[::self.subsample,::self.subsample],cmap=colormap,levels=clevs,animated=True,alpha=.75)
    #     plt.colorbar(label='Temperature $(^\circ C)$',fraction=0.046, pad=0.04,ticks=[-1.5, -.5, .5,1.5])




    # def eaverage(self,level='bottom',day=0):
    #     subsample = self.subsample
    #     if level is 'top':
    #         level=0
    #     elif level is 'bottom':
    #         level = 23
    #     with open(os.path.join(datadir,'low_res_sose_mat.pickle'),'rb') as handle: 
    #         mat = pickle.load(handle)
    #     XC = mat['XC'][:,0]
    #     YC = mat['YC'][0,:]
    #     with open(os.path.join(datadir,'Eaverage_'+str(level)+'.pickle'),'rb') as handle: 
    #         plot_data = pickle.load(handle)    
    #     if self.lon_mid==0:
    #         plot_data, XC = shiftgrid(180, plot_data, XC, start=False)
    #     XC,YC = np.meshgrid(XC,YC)
    #     XX,YY = self(XC,YC)
    #     self.contourf(XX[::subsample,::subsample],YY[::subsample,::subsample],plot_data[::subsample,::subsample],cmap=plt.cm.Greys,norm = LogNorm(),animated=True)
    #     plt.colorbar(label='J/kg',fraction=0.046, pad=0.04)



    # def eke(self,level='bottom',day=0):
    #     subsample = self.subsample
    #     if level is 'top':
    #         level=0
    #     elif level is 'bottom':
    #         level = 23
    #     with open(os.path.join(datadir,'low_res_sose_mat.pickle'),'rb') as handle: 
    #         mat = pickle.load(handle)
    #     XC = mat['XC'][:,0]
    #     YC = mat['YC'][0,:]
    #     with open(os.path.join(datadir,'EKE_'+str(level)+'.pickle'),'rb') as handle: 
    #         data = pickle.load(handle)    
    #     plot_data = data[day,:,:]
    #     if self.lon_mid==0:
    #         plot_data, XC = shiftgrid(180, plot_data, XC, start=False)
    #     XC,YC = np.meshgrid(XC,YC)
    #     XX,YY = self(XC,YC)
    #     self.contourf(XX[::subsample,::subsample],YY[::subsample,::subsample],plot_data[::subsample,::subsample],cmap=plt.cm.Greys,norm = LogNorm(),animated=True,alpha=0.6)
    #     plt.colorbar(label='Mean Eddy Kinetic Energy (J/kg)',fraction=0.046, pad=0.04,shrink=0.5)



    # def etotal(self,level='bottom',day=0):
    #     subsample = self.subsample
    #     if level is 'top':
    #         level=0
    #     elif level is 'bottom':
    #         level = 23
    #     with open(os.path.join(datadir,'low_res_sose_mat.pickle'),'rb') as handle: 
    #         mat = pickle.load(handle)
    #     XC = mat['XC'][:,0]
    #     YC = mat['YC'][0,:]
    #     with open(os.path.join(datadir,'Etotal_'+str(level)+'.pickle'),'rb') as handle: 
    #         data = pickle.load(handle)    
    #     plot_data = data[day,:,:]
    #     if self.lon_mid==0:
    #         plot_data, XC = shiftgrid(180, plot_data, XC, start=False)
    #     XC,YC = np.meshgrid(XC,YC)
    #     XX,YY = self(XC,YC)
    #     self.contourf(XX[::subsample,::subsample],YY[::subsample,::subsample],plot_data[::subsample,::subsample],cmap=plt.cm.Greys,norm = LogNorm(),animated=True)
    #     plt.colorbar(label='J/kg',fraction=0.046, pad=0.04)



    # def quiver_plot(self,level='bottom'):
    #     subsample = self.subsample*4+5
    #     if level is 'top':
    #         level=0
    #     elif level is 'bottom':
    #         level = 23
    #     U = np.load(os.path.join(datadir,'averageU_'+str(level)+'.npy'))
    #     V = np.load(os.path.join(datadir,'averageV_'+str(level)+'.npy'))
    #     mat = scipy.io.loadmat(os.path.join(datadir,'grid.mat'))
    #     XC = mat['XC'][:,0]
    #     YC = mat['YC'][0,:]
    #     if self.lon_mid==0:
    #         U, dummy = shiftgrid(180, U, XC, start=False)
    #         V, XC = shiftgrid(180, V, XC, start=False)
    #     XC,YC = np.meshgrid(XC,YC)
    #     XX,YY = self(XC,YC)
    #     V_mask = np.ma.masked_equal(V,0.)
    #     U_mask = np.ma.masked_equal(U,0.)
    #     self.quiver(XX[::subsample,::subsample],YY[::subsample,::subsample],U_mask[::subsample,::subsample],V_mask[::subsample,::subsample],linewidths=0.5)

    # def SOSO(self,markersz = 15):
    #     x = [a[1] for a in SOSO_weddel_coord.values()]
    #     x1 = [a[1] for a in SOSO_drake_coord.values()]
    #     y = [a[0] for a in SOSO_weddel_coord.values()]
    #     y1 = [a[0] for a in SOSO_drake_coord.values()]
    #     xpos, ypos = self(x,y)
    #     self.plot(xpos,ypos,'*',markersize = markersz, color = 'gold',label = 'Weddell Sound Sources')
    #     xpos, ypos = self(x1,y1)
    #     self.plot(xpos,ypos,'s',markersize = markersz/2, color = 'gold',label = 'Dimes Sound Sources') #,label = 'Sound Sources')

    # def LOP_plot(self):
    #     soso_df = self.soso_df[(self.soso_df.Date==self.date)&(self.soso_df.Cruise==self.cruise)]
    #     for it in soso_df.iterrows():
    #         if it[1]['PosQC'] in SOSO_weddel_coord: 
    #             # print SOSO_weddel_coord[it[1]['PosQC']][0]
    #             # print SOSO_weddel_coord[it[1]['PosQC']][1]
    #             self.equi(SOSO_weddel_coord[it[1]['PosQC']][1],SOSO_weddel_coord[it[1]['PosQC']][0], self.toa_to_dist(it[1]['Observation']),lw=2.,label='source '+it[1]['PosQC'])
    #         if it[1]['PosQC'] in SOSO_drake_coord: 
    #             self.equi(self, SOSO_drake_coord[it[1]['PosQC']][1],SOSO_drake_coord[it[1]['PosQC']][0], self.toa_to_dist(it[1]['Observation']),lw=2.,label='source '+it[1]['PosQC'])

    # def covariance_plot(self):
    #     w = df[df.Date==self.date].eig1.values[0]
    #     v = df[df.Date==self.date].eig2.values[0]
    #     print w
    #     print v

    #     xcoord = df[df.Date==self.date].Lon.values[0]
    #     ycoord = df[df.Date==self.date].Lat.values[0]
    #     angle = np.degrees(np.arctan(v[1,np.argmax(w)]/v[0,np.argmax(w)]))
    #     print angle
    #     poly = self.ellipse(xcoord, ycoord, 2*max(w)*np.sqrt(5.991),2*min(w)*np.sqrt(5.991),400, phi=angle, facecolor='c', zorder=10,
    #         alpha=0.3)
        
    # def region_plot(self,region=[]):
    #     regions=[(-60,20,-76,-60,'weddell'),(140,280,-76,-60,'ross')]
    #     if region:
    #         regions = [x for x in regions if x[4] in region]
    #     for region in regions:
    #         lons = np.array([region[1],region[1]]+np.linspace(region[1],region[0],200).tolist()+[region[0],region[0],region[1]])
    #         lons[lons>180] = lons[lons>180]-360
    #         print lons
    #         lats = [region[2],region[3]]+([region[3]]*200)+[region[3],region[2],region[2]]
    #         self.plot(lons.tolist(),lats,'k--',latlon=True,zorder=11,linewidth=7)


    # def show(self):
    #     plt.show()