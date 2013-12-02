# -*- coding:utf-8 -*-
"""
@author Geir Arne Waagbø
@see http://code.google.com/p/pywrfplot/
"""
import matplotlib as M
M.use('agg')
import matplotlib.pyplot as plt
from numpy import arange
from mpl_toolkits.basemap import Basemap 
from Params import *
from Utils import *
import numpy as N
import pdb
from mesowestscripts import WRFterrain

import sys
sys.path.append('/uufs/chpc.utah.edu/common/home/u0737349/lawsonpython/')
    
# plotting
height, width = (10,10)
plt.rc('text',usetex=False)
#fonts = {'family':'Computer Modern','size':16}
#plt.rc('font',**fonts)

def mapDomains(gl,terrain,includePMSL=False,drawxs=False):
    """Creates a map of the outer domain, showing
    the location of the inner domains.
    If includePMSL is True contours of sea level pressure is
    plotted as well (if a WPS-file is found)
    
    """
    plt.figure(figsize=(width,height))
    
    nc1 = openWRF(1)
    nc2 = openWRF(2)
    nc3 = openWRF(3)
    nc4 = openWRF(4)
    m = _getMapForNC(nc1)

    if includePMSL:
        met1 = openWPS(1)
        if met1 is not None:
            PMSL = met1.variables['PMSL'][0]
            _Nx1,_Ny1,_Nz1,longitude1,latitude1,_dx,_dy,_x,_y = getDimensions(nc1)
            x, y = m(longitude1,latitude1)
            cs = plt.contour(x,y, PMSL/100.,pmsl_int,colors='black')
            plt.clabel(cs, inline=1,  fmt='%1.0f', fontsize=11)
        else:
            print 'Could not find wps-file, PMSL not plotted'
    
    if (nc1 is not None):
        _plotBorder(nc1,m,gl[0])
    if (nc2 is not None):
        _plotBorder(nc2,m,gl[1])
    if (nc3 is not None):
        _plotBorder(nc3,m,gl[2])
    if (nc4 is not None):
        _plotBorder(nc4,m,gl[3])
   
    if terrain: 
        TC = (500,4500,500) # Terrain contours
        terrain, tlon, tlat = WRFterrain()
        topodat = N.where(terrain < 1,None,terrain)
        topodatm = N.ma.masked_values(topodat, None)
        tlonx, tlatx = m(tlon,tlat)
        m.contourf(tlonx, tlatx, topodatm, cmap=M.cm.Greys,extend='both')
        cb = plt.colorbar(orientation='horizontal')
        cb.set_label('Terrain(m)')
 
    # Draw cross-sections
    if drawxs:
        for xs in drawxs.keys():
            Alat, Alon = drawxs[xs]['startll']
            Blat, Blon = drawxs[xs]['endll']
            #Axpt, Aypt = m(Alon,Alat)
            #Bxpt, Bypt = m(Blon,Blat)
            m.drawgreatcircle(Alon,Alat,Blon,Blat,del_s=1,linewidth=2)
    
    fname = outdir+'domains.png'
    print fname
    plt.savefig(fname)#,bbox_inches='tight',pad_inches=0.5)
    plt.clf()

def mapWind(nest,time):
    """Creates a map of the domain, showing
    wind barbs for the given time
    """
    nc = openWRF(nest)
    nc1 = openWRF(nest+1)
    Nx,Ny,_Nz,longitude,latitude,_dx,_dy,_x,_y = getDimensions(nc)
    m = _getMapForNC(nc,False,_getDL(nest),100)
    _makeDots(m)
    u10 = nc.variables['U10'][time,:,:]
    v10 = nc.variables['V10'][time,:,:]
    # Use data from every 10th grid point
    windx = 1+Nx/10
    windy = 1+Ny/10
    lat10 = np.zeros((windy,windx))
    lon10 = np.zeros((windy,windx))
    uwind = np.zeros((windy,windx))
    vwind = np.zeros((windy,windx))
    for j in range(windy):
        for i in range(windx):
            uwind[j,i] = 0.5*(u10[j*10,i*10]+u10[j*10,i*10+1])
            #print 'u: ' + str(uwind[j,i]) 
            vwind[j,i] = 0.5*(v10[j*10,i*10]+v10[j*10+1,i*10])
            #print 'v: ' + str(vwind[j,i]) 
            lat10[j,i] = latitude[j*10,i*10]
            lon10[j,i] = longitude[j*10,i*10]

    x10,y10  = m(lon10,lat10)
    plt.barbs(x10,y10,uwind,vwind,barb_increments=barb_increments,linewidth=1.0,color='green')

    if (nc1 is not None):
        _plotBorder(nc1,m,'black')
    plt.show()  
    plt.close()

def mapTerrain(nest=1,includePMSL=False,includeTemp=False):
    """Creates a map of the domain, showing
    filled contours of the terrain
    """
    nc = openWRF(nest)
    met = openWPS(nest)
    nc1 = openWRF(nest+1)
    _Nx,_Ny,_Nz,longitude,latitude,_dx,_dy,_x,_y = getDimensions(nc)

    m = _getMapForNC(nc,False,_getDL(nest),100)
    x, y = m(longitude,latitude)
    _makeDots(m)
    
    if (includeTemp and met is not None):
        ST = met.variables['TT'][0][0,:,:]
        cs = plt.contour(x,y, ST-T_zero,arange(-48,50,2.),colors='blue',linestyles='solid')
        plt.clabel(cs, inline=1,  fmt='%1.0f', fontsize=12)
    if (includePMSL and met is not None):
        levels = arange(960.,1040.,1.)
        PMSL = met.variables['PMSL'][0]
        cs = plt.contour(x,y, PMSL/100.,levels,colors='black',label='Pressure')
        plt.clabel(cs, inline=1,  fmt='%1.0f', fontsize=11)

    heightground = (nc.variables['PH'][0][0,:,:] + nc.variables['PHB'][0][0,:,:])/g
    plt.contourf(x, y, heightground, levels=arange(0,max_h,100),cmap=cmap_grey,extend='max')
    plt.colorbar()
    if (nc1 is not None):
        _plotBorder(nc1,m,'black')
    #plt.show()  
    plt.savefig(outdir+naming+'domains_terrain.png')
    plt.close()

def mapCloud(nest,time):
    """Creates a map of the domain, showing
    filled contours of the cloud water
    """
    nc = openWRF(nest)
    _Nx,_Ny,_Nz,longitude,latitude,_dx,_dy,_x,_y = getDimensions(nc)

    m = _getMapForNC(nc,False,_getDL(nest),100)
    x, y = m(longitude,latitude)

    qcloud = 1000.0*np.sum(nc.variables['QCLOUD'][time,:,:,:],axis=0)
    plt.contourf(x, y, qcloud, cmap=cmap_red)
    plt.colorbar()
    #plt.show()  
    plt.close()

def _getMapForNC(nc,fill=False,dl=5,at=1000):
    Nx,Ny,_Nz,longitude,latitude,dx,dy,_x,_y = getDimensions(nc)
    lon = longitude[Ny/2,:]
    lat = latitude[:,Nx/2]
    lon_min = np.round(lon[0],0)
    lon_max = np.round(lon[Nx-1],0)+1
    lat_min = np.round(lat[0],0)
    lat_max = np.round(lat[Ny-1],0)+1
    return _getMap(fill,dx,dy,Nx,Ny,mapResolution,lat_focuspoint,lon_focuspoint,lat_min,lat_max,lon_min,lon_max,dl,at)

def _getMap(fill,dx,dy,Nx,Ny,res,lat,lon,lat_min,lat_max,lon_min,lon_max,dl,at):
    m = Basemap(width=dx*Nx+8e5,height=dy*Ny+5e5,
                resolution=res,area_thresh=at,projection='stere',
                lat_1=lat,lat_2=lat,
                lat_0=lat,lon_0=lon)
    m.drawcoastlines()
    if fill:
        m.fillcontinents(color=LandColor,lake_color=WaterColor)
        m.drawmapboundary(fill_color=WaterColor)
    m.drawstates()
    m.drawcountries() 
    m.drawparallels(np.arange(lat_min,lat_max,dl),labels=[1,0,0,0])
    m.drawmeridians(np.arange(lon_min,lon_max,dl),labels=[0,0,0,1])
    return m

def _plotBorder(nc,map,gridlabel,color='black'):
    Nx,Ny,_Nz,longitude,latitude,_dx,_dy,_x,_y = getDimensions(nc)
    x,y = map(longitude,latitude)
    plt.plot(x[0,:],y[0,:],color,lw=2)
    
    # Label domain
    xl = len(x[0,:])
    midpt = len(y[0,:])/2
    plt.annotate(gridlabel,color=color,size=21,xy=(x[0,-(0.12*xl)],y[0,midpt]),bbox=dict(fc='white'),alpha=1,va='center',ha='left')
    plt.plot(x[:,0],y[:,0],color,lw=2)
    plt.plot(x[Ny-1,:],y[Ny-1,:],color,lw=2)
    plt.plot(x[:,Nx-1],y[:,Nx-1],color,lw=2)

def _makeDots(m):
    x_f,y_f = m(lon_focuspoint, lat_focuspoint)
    plt.plot(x_f,y_f,'ro')
    if (lon_rg != -1 and lat_rg != -1):
            x_rg, y_rg = m(lon_rg, lat_rg)
            plt.plot(x_rg,y_rg,'ro')

def _getDL(nest):
    # Tick increment for latitude/longitude
    if nest==4:
        dl = 0.2
    elif nest == 3:
        dl = 1.
    elif nest == 2:
        dl = 2.
    else:
        dl = 4.       
    return dl    
