from GeneralUtilities.Plot.Cartopy.eulerian_plot import BaseCartopy

class RegionalBase(BaseCartopy):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.ax.set_extent([self.llcrnrlon,self.urcrnrlon,self.llcrnrlat,self.urcrnrlat], crs=ccrs.PlateCarree())
        self.finish_map()

class SOSECartopy(RegionalBase):
    llcrnrlon=-180.
    llcrnrlat=-80.
    urcrnrlon=180.
    urcrnrlat=-25
    def __init__(self,*args,**kwargs):
        print('I am plotting antarctic region')
        super().__init__(*args,**kwargs)

class CreteCartopy(RegionalBase):
    llcrnrlon=20.
    llcrnrlat=30
    urcrnrlon=30
    urcrnrlat=40
    def __init__(self,*args,**kwargs):
        print('I am plotting Crete')
        super().__init__(*args,**kwargs)

class KonaCartopy(RegionalBase):
    llcrnrlon=-157
    llcrnrlat=18.8
    urcrnrlon=-155.6
    urcrnrlat=20.35
    def __init__(self,*args,**kwargs):
        print('I am plotting Kona')
        super().__init__(*args,**kwargs)

        self.ax.set_extent([llcrnrlon,urcrnrlon,llcrnrlat,urcrnrlat], crs=ccrs.PlateCarree())
        self.finish_map()

class PuertoRicoCartopy(RegionalBase):
    llcrnrlon=-68.5 
    llcrnrlat=16
    urcrnrlon=-65
    urcrnrlat=22.5
    def __init__(self,*args,**kwargs):
        print('I am plotting Puerto Rico')
        super().__init__(*args,**kwargs)

class MobyCartopy(RegionalBase):
    center_lat = 20.8
    center_lon = -157.2
    llcrnrlon=(self.center_lon-3)
    llcrnrlat=(self.center_lat-3)
    urcrnrlon=(self.center_lon+3)
    urcrnrlat=(self.center_lat+3)

    def __init__(self,*args,**kwargs):
        print('I am plotting Moby')
        super().__init__(*args,**kwargs)
        self.ax.scatter(center_lon,center_lat,500,marker='*',color='Red',zorder=10)

class GOMCartopy(RegionalBase):
    llcrnrlon=-100.
    llcrnrlat=20.5
    urcrnrlon=-81.5
    urcrnrlat=30.5
    def __init__(self,*args,**kwargs):
        print('I am plotting GOM')
        super().__init__(*args,**kwargs)

class CCSCartopy(RegionalBase):
    llcrnrlon=-135.
    llcrnrlat=20
    urcrnrlon=-105
    urcrnrlat=55
    def __init__(self,*args,**kwargs):
        print('I am plotting GOM')
        super().__init__(*args,**kwargs)

class DrakePassageCartopy(RegionalBase):
    llcrnrlon=-130.
    llcrnrlat=-70.
    urcrnrlon=0.
    urcrnrlat=-35.
    def __init__(self,*args,**kwargs):
        print('I am plotting GOM')
        super().__init__(*args,**kwargs)