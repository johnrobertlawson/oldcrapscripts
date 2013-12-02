# -*- coding:utf-8 -*-
"""
@author Geir Arne Waagbø
@see http://code.google.com/p/pywrfplot/
"""
import matplotlib as M
import matplotlib.pyplot as plt
from numpy import arange
import pdb
import numpy as N
import math

from pywrfplotParams import *
from pywrfplotUtils import *
import customcmaps

"""
Is it possible to change the cross-section to e.g. NW to SE...?
"""

### USER OPTIONS for figures
height, width = (9,17)

def xzCloudPlot(nest,datetuple,latlon, plotWind, plotTemp=False,plotRH=False,plotTheta=True, plotW=False):
    nc = openWRF(nest)
    
    plt.rc('text',usetex=True)
    fonts = {'family':'Computer Modern','size':16}
    plt.rc('font',**fonts)
    
    time = convert_time(nest,datetuple) # JRL function to convert time stamp to time
    Nx,Ny,Nz,longitude,_lats,_dx,_dy,x_nr,y_nr = getDimensions(nc)
    
    # Find nearest grid points along a line between two lat/lon   
    latxs = (latlon[0][0],latlon[1][0])
    lonxs = (latlon[0][1],latlon[1][1])
    lats = _lats.astype(N.float32)
    lons = longitude.astype(N.float32)
    # x-y position of start and end points
    stxy = getXYxs(lons[Ny/2,:],lats[:,Nx/2],latlon[0])
    endxy = getXYxs(lons[Ny/2,:],lats[:,Nx/2],latlon[1])
    # Cross-section is length of dx or dy, whichever is longer
    #dx = abs(stxy[0]-endxy[0])
    #dy = abs(stxy[1]-endxy[1])
    #dxs = [dx if dx > dy else dy][0]
    
    # Joe Kington method
    hypotenuse = int(N.hypot(endxy[0]-stxy[0], endxy[1]-stxy[1]))
    xx = N.linspace(stxy[0],endxy[0],hypotenuse)
    yy = N.linspace(stxy[1],endxy[1],hypotenuse)
    xint, yint = xx.astype(int), yy.astype(int)
    angle = N.arctan((yy[-1]-yy[0])/(xx[-1]-xx[0])) # In radians

    #x_nr and y_nr are the x and y locs for x-section centre
    heightground_x,heighthalf_xz = _getHeightJRL(nc,nest, time, yy, xx, Nz, hypotenuse)    
    #print 'Model height: ' + str(heightground_x[x_nr])
    
    #pdb.set_trace()
    
    #xslats = N.linspace(latxs[0],latxs[1],dxs)
    #xslons = N.linspace(lonxs[0],lonxs[1],dxs)
    #latf = lats.flatten()
    #lonf = lons.flatten()
    #lat2d, lon2d = N.meshgrid(latf,lonf,copy=False)
    #for ns in range(0,dxs):
        # Convert negative lons to positive
    #    lon_2 = xslons[ns] 
    #    ijll = []
    #    if lon_2 < 0:
    #        lon_2 =+ 360
    #    dlatlon = (lon2d-lon_2)**2 + (lat2d-xslats(ns))**2
    #    ijll.append(N.where(dlatlon==dlatlon.min()).min())

    #pdb.set_trace() 
    
    if plotTheta==True:
        theta_xz = nc.variables['T'][time,:,yy.astype(N.int),xx.astype(N.int)] + T_base 
    if plotW:
        w_xz = nc.variables['W'][time,:-1,yy.astype(N.int),xx.astype(N.int)]
    #    theta = 
    #if plot
    #    P = nc.variables['P'][time,:,y_nr,:] + nc.variables['PB'][time,:,y_nr,:] 
    #T = theta*(P/P_bot)**kappa # Temperatur i halvflatene (Kelvin)
    #rho = P/(R*T) #[kg/m3]   

    #qcloud_xz = 1000.0*nc.variables['QCLOUD'][time,:,y_nr,:]*rho # regner om til g/m3
    #qrain_xz = 1000.0*nc.variables['QRAIN'][time,:,y_nr,:]*rho 
    #qsnow_xz = 1000.0*nc.variables['QSNOW'][time,:,y_nr,:]*rho 
    
    # This is a very lazy/bad way to deal with the staggering
    if plotWind:
        #u_xz = nc.variables['U'][time,:,y_nr,:-1] # Staggering causes one more in the 4th dimension
        #v_xz = nc.variables['V'][time,:,y_nr,:] 
        u_xz = nc.variables['U'][time,:,yy.astype(N.int),xx.astype(N.int)] # NOT SURE HERE [:-1] or -1
        v_xz = nc.variables['V'][time,:,yy.astype(N.int),xx.astype(N.int)]
    if plotWind == 'parallel':
        urot = (N.cos(angle)*u_xz) + (N.sin(angle)*v_xz)
    if plotWind == 'normal':
        vrot = (N.cos(angle)*v_xz) - (N.sin(angle)*v_xz)
    if plotWind == 'magnitude':
        wind_xz = N.sqrt(u_xz**2 + v_xz**2)
        #planewind = N.sqrt(urot**2 + vrot**2)
        # Decompose into plane-parallel
        #planewindf = N.asarray([math.atan2(u,v) for v,u in zip(v_xz.flatten(),u_xz.flatten())])
        #planewind = N.reshape(planewindf,wind_xz.shape)

    hgl = nc.DX/N.cos(angle) # Hypotenuse grid length along cross-section
    hgl_ticks = [(x*hgl)/1000.0 for x in range(0,hypotenuse)] # In km
    hgl_labels = ["%2.1f" %x for x in hgl_ticks]
    grid = N.swapaxes(N.repeat(N.array(hgl_ticks).reshape(hypotenuse,1),Nz,axis=1),0,1)
       
    # Colorbar custom load
    my_cmap = M.colors.LinearSegmentedColormap('my_colormap',customcmaps.jet,256)
    norml = M.colors.BoundaryNorm(w_wind_levels,256)
     
    #pdb.set_trace()
    plt.figure(figsize=(width, height))
    plt.set_cmap(cmap_red)
    #plt.axis([0,Nx-1,z_min,z_max]) #JRL: edited in z_min
    plt.axis([0,hypotenuse,z_min,z_max]) #JRL: edited in z_min
    #print u'Cloud water red, snow blue, rain green ($g/m^3$)'
    #grid = N.reshape(N.tile(arange(hypotenuse),Nz),(Nz,-1))
    #plt.contourf(grid, heighthalf_xz, qcloud_xz, alpha=0.9,levels=xz_cloudwater_levels, cmap=cmap_red)#
    #plt.colorbar()
    #plt.contourf(grid, heighthalf_xz, qrain_xz, alpha=0.6,levels=xz_rain_levels, cmap=cmap_green)#
    #plt.colorbar()
    #plt.contourf(grid, heighthalf_xz, qsnow_xz, alpha=0.6,levels=xz_snow_levels,cmap=cmap_blue)# 
    #plt.colorbar()
    #plt.contourf(grid, heighthalf_xz, wind_xz, alpha=0.6, levels=xz_wind_levels,cmap=plt.cm.jet)
    if plotWind == 'parallel':
        plt.contourf(grid, heighthalf_xz, urot, alpha=0.6, levels=u_wind_levels,cmap=plt.cm.jet)
    if plotWind == 'normala':
        plt.contourf(grid, heighthalf_xz, vrot, alpha=0.6, levels=u_wind_levels,cmap=plt.cm.jet)
    if plotWind == 'magnitude':
        plt.contourf(grid, heighthalf_xz, wind_xz, alpha=0.6, levels=xz_wind_levels,cmap=plt.cm.jet)
    #plt.contourf(grid, heighthalf_xz, theta, alpha=0.6, levels=xz_theta_levels,cmap=plt.cm.jet)
    #pdb.set_trace()
    if plotW:
        cfW = plt.contourf(grid, heighthalf_xz, w_xz, alpha=0.6, levels=w_wind_levels,cmap=my_cmap,norm=norml,extend='both')
        #    cfW.cmap.set_under('blue')   
        #    cfW.cmap.set_over('red')
    #plt.contourf(grid, heighthalf_xz, w_xz, alpha=0.6,cmap=plt.cm.jet)
    if plotW == False:
        cb = plt.colorbar(orientation="vertical")
    else:
        cb = plt.colorbar(cfW)
    if plotWind:
        cb.set_label(plotWind + ' wind speed (m/s)')    
    if plotW:
        cb.set_label('Vertical wind speed (m/s)')    
    if plotTemp:
        temp_int = arange(-80.0,50.0,2.0)
        cs = plt.contour(grid, heighthalf_xz, T-T_zero, temp_int,colors='black',linestyles='solid')#linewidths=4
        plt.clabel(cs, inline=1,  fmt='%1.0f', fontsize=12,colors='black')
    if plotRH:
        rh = _getRH(nc,time,-1,y_nr,T,P)
        rh_int = arange(90.,111.,5.)
        cs = plt.contour(grid, heighthalf_xz,rh , rh_int, colors='grey')
        plt.clabel(cs, inline=1,  fmt='%1.0f', fontsize=12, colors='grey')
    if plotTheta:
        theta_int = arange(260,400,2)
        cs = plt.contour(grid, heighthalf_xz, theta_xz, theta_int, colors='black')
        plt.clabel(cs, inline=1,  fmt='%1.0f', fontsize=12, colors='grey')

    plt.plot(hgl_ticks,heightground_x,color='black')
    plt.fill_between(hgl_ticks,heightground_x,0,facecolor='lightgrey')
    
    #pdb.set_trace()
    labeldelta = 15
    plt.xticks(np.arange(0,hypotenuse,labeldelta),hgl_labels[::labeldelta])
    plt.yticks(np.arange(z_min,z_max,dz)) # JRL: edited in z_min
    plt.xlabel('Distance along cross-section (km)')
    plt.ylabel(u'Height above sea level (m)')
    plt.title('Potential temperature at ' + str(datetuple[3]) + 'Z on ' + str(datetuple[:3]))
    #plt.show()
    fname = outdir+naming+str(nest)+'_'+str(datetuple[2])+str(datetuple[3])+'Z_'+'W'+'_xsection.png'
    plt.savefig(fname)
    with open(outdir+"xz_timelist.txt","a") as myfile:
        myfile.write(fname + '\n')
    plt.close()        

def yzCloudPlot(nest,time,plotTemp=True,plotRH=False):
    nc = openWRF(nest)
    #time = convert_time(nest,time) # JRL function to convert time stamp to time
    Nx,Ny,Nz,_longs,latitude,_dx,_dy,x_nr,y_nr = getDimensions(nc)
    
    heightground_y,heighthalf_yz = _getHeight(nc,nest, time, -1, Ny, Nz, x_nr,-1)    
    print 'Model height: ' + str(heightground_y[y_nr])

    theta = nc.variables['T'][time,:,:,x_nr] + T_base 
    P = nc.variables['P'][time,:,:,x_nr] + nc.variables['PB'][time,:,:,x_nr] 
    T = theta*(P/P_bot)**kappa # Temperatur i halvflatene (Kelvin)
    rho = P/(R*T) #[kg/m3]

    qcloud_yz = 1000.0*nc.variables['QCLOUD'][time,:,:,x_nr]*rho # regner om til g/m3
    qrain_yz = 1000.0*nc.variables['QRAIN'][time,:,:,x_nr]*rho 
    #qsnow_yz = 1000.0*nc.variables['QSNOW'][time,:,:,x_nr]*rho 
   
    plt.figure()
    plt.set_cmap(cmap_red)
    plt.axis([0,Ny-1,z_min,z_max]) # JRL: Edited in z_min
    print u'Cloud water red, snow blue, rain green ($g/m^3$)'
    grid = np.reshape(np.tile(arange(Ny),Nz),(Nz,-1))
    plt.contourf(grid, heighthalf_yz, qcloud_yz, alpha=0.9,levels=xz_cloudwater_levels, cmap=cmap_red)#
    plt.colorbar()
    plt.contourf(grid, heighthalf_yz, qrain_yz, alpha=0.6,levels=xz_rain_levels, cmap=cmap_green)#
    plt.colorbar()
    plt.contourf(grid, heighthalf_yz, qsnow_yz, alpha=0.6,levels=xz_snow_levels,cmap=cmap_blue)# 
    plt.colorbar()

    if plotTemp:
        cs = plt.contour(grid, heighthalf_yz, T-T_zero, temp_int,colors='black',linestyles='solid')#linewidths=4
        plt.clabel(cs, inline=1,  fmt='%1.0f', fontsize=12,colors='black')
    if plotRH:
        rh = _getRH(nc,nest,time,x_nr,-1,T,P)
        rh_int = arange(90.,111.,5.)
        cs = plt.contour(grid, heighthalf_yz,rh , rh_int,colors='grey')
        plt.clabel(cs, inline=1,  fmt='%1.0f', fontsize=12,colors='grey')
    plt.plot(arange(Ny),heightground_y,color='black')
    plt.fill_between(arange(Ny),heightground_y,0,facecolor='lightgrey')
    plt.xticks(np.arange(0,Ny,8),np.round(latitude[::8,Nx/2], 1), fontsize='small')
    plt.yticks(np.arange(z_min,z_max,dz), fontsize='small') #JRL: edited in z_min
    plt.xlabel('Breddegrad')
    plt.ylabel(u'Høyde [m]')
    plt.show()
    plt.close()        

def _getHeight(nc,nest,time,Nx,Ny,Nz,x_nr,y_nr):
    #time = convert_time(nest,time) # JRL function to convert time stamp to time
    # Calculate height above sea level for mass levels
    # NB! Either x_nr must be -1 or y_nr must be -1
    # Note: geopotential defined at full levels (w-levels)
    #       Must interpolate to find geopotential at half levels (u-levels)
    geopot = (nc.variables['PH'][time,0:Nz,y_nr,:] + nc.variables['PHB'][time,0:Nz,y_nr,:]) if y_nr!=-1 else \
             (nc.variables['PH'][time,0:Nz,:,x_nr] + nc.variables['PHB'][time,0:Nz,:,x_nr])
    mu = (nc.variables['MU'][time,y_nr,:]+nc.variables['MUB'][time,y_nr,:]) if y_nr!=-1 else \
         (nc.variables['MU'][time,:,x_nr]+nc.variables['MUB'][time,:,x_nr])
    znw = nc.variables['ZNW'][time,0:Nz] # full (w) levels
    znu = nc.variables['ZNU'][time,0:Nz] # half (u,mass) levels

    heighthalf = np.zeros((Nz,Nx if y_nr!=-1 else Ny))# height in meters
    for i in arange(Nx if y_nr!=-1 else Ny):
        pfull = mu[i]*znw+P_top
        phalf = mu[i]*znu+P_top
        for k in arange(Nz):
            heighthalf[k,i]=interp(geopot[:,i],pfull[:],phalf[k])/g
    heightground = geopot[0,:]/g
    return heightground,heighthalf

def _getRH(nc,nest,time,x_nr,y_nr,T,P):
    #time = convert_time(nest,time) # JRL function to convert time stamp to time
    es_w = es(T-T_zero) # metningstrykk (Pascal)
    qsat = eps*es_w/(P-0.378*es_w)
    #es_ice = es_w*(T/T_zero)**2.66
    #qsat_ice = eps*es_ice/(P-0.378*es_w)
    qvapor = nc.variables['QVAPOR'][time,:,y_nr,:] if y_nr!=-1 else nc.variables['QVAPOR'][time,:,:,x_nr]
    rh = 100.*qvapor*(1-qsat)/(qsat*(1-qvapor))
    #rh_ice = 100.*qvapor*(1-qsat_ice)/(qsat_ice*(1-qvapor))
    return rh

def _getHeightJRL(nc,nest,time,xx,yy,Nz,hypotenuse):    
    xx = xx.astype(N.int)                                                             
    yy = yy.astype(N.int)                                                             
    #time = convert_time(nest,time) # JRL function to convert time stamp to time
    # Calculate height above sea level for mass levels
    # NB! Either x_nr must be -1 or y_nr must be -1
    # Note: geopotential defined at full levels (w-levels)
    #       Must interpolate to find geopotential at half levels (u-levels)
    geopot = nc.variables['PH'][time,0:Nz,xx,yy] + nc.variables['PHB'][time,0:Nz,xx,yy]
    mu = nc.variables['MU'][time,xx,yy] + nc.variables['MUB'][time,xx,yy]
    znw = nc.variables['ZNW'][time,0:Nz] # full (w) levels
    znu = nc.variables['ZNU'][time,0:Nz] # half (u,mass) levels

    heighthalf = np.zeros((Nz,hypotenuse))# height in meters
    for i in arange(hypotenuse):
        pfull = mu[i]*znw+P_top
        phalf = mu[i]*znu+P_top
        for k in arange(Nz):
            heighthalf[k,i]=interp(geopot[:,i],pfull[:],phalf[k])/g
    heightground = geopot[0,:]/g
    return heightground,heighthalf
