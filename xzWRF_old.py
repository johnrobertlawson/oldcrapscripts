# -*- coding:utf-8 -*-
"""
@author Geir Arne Waagbø
@see http://code.google.com/p/pywrfplot/
"""
import matplotlib.pyplot as plt
from numpy import arange
import pdb
import numpy as N

from pywrfplotParams import *
from pywrfplotUtils import *

"""
Is it possible to change the cross-section to e.g. NW to SE...?
"""

### USER OPTIONS for figures
height, width = 9, 17

def xzCloudPlot(nest,datetuple,latlon,plotTemp=False,plotRH=False,plotTheta=True, plotWind=True):
    nc = openWRF(nest)

    time = convert_time(nest,datetuple) # JRL function to convert time stamp to time
    Nx,Ny,Nz,longitude,_lats,_dx,_dy,x_nr,y_nr = getDimensions(nc)
    
    #x_nr and y_nr are the x and y locs for x-section centre
    heightground_x,heighthalf_xz = _getHeight(nc,nest, time, Nx, -1, Nz, -1, y_nr)    
    print 'Model height: ' + str(heightground_x[x_nr])
    
    if plotTheta==True:
        theta = nc.variables['T'][time,:,y_nr,:] + T_base 
    #if plot
    #    P = nc.variables['P'][time,:,y_nr,:] + nc.variables['PB'][time,:,y_nr,:] 
    #T = theta*(P/P_bot)**kappa # Temperatur i halvflatene (Kelvin)
    #rho = P/(R*T) #[kg/m3]   

    #qcloud_xz = 1000.0*nc.variables['QCLOUD'][time,:,y_nr,:]*rho # regner om til g/m3
    #qrain_xz = 1000.0*nc.variables['QRAIN'][time,:,y_nr,:]*rho 
    #qsnow_xz = 1000.0*nc.variables['QSNOW'][time,:,y_nr,:]*rho 
    
    # This is a very lazy/bad way to deal with the staggering
    if plotWind == True:
        u_xz = nc.variables['U'][time,:,y_nr,:-1] # Staggering causes one more in the 4th dimension
        v_xz = nc.variables['V'][time,:,y_nr,:] 
        wind_xz = N.sqrt(u_xz**2 + v_xz**2)

    plt.figure(figsize=(width, height))
    plt.set_cmap(cmap_red)
    plt.axis([0,Nx-1,z_min,z_max]) #JRL: edited in z_min
    #print u'Cloud water red, snow blue, rain green ($g/m^3$)'
    grid = np.reshape(np.tile(arange(Nx),Nz),(Nz,-1))
    #plt.contourf(grid, heighthalf_xz, qcloud_xz, alpha=0.9,levels=xz_cloudwater_levels, cmap=cmap_red)#
    #plt.colorbar()
    #plt.contourf(grid, heighthalf_xz, qrain_xz, alpha=0.6,levels=xz_rain_levels, cmap=cmap_green)#
    #plt.colorbar()
    #plt.contourf(grid, heighthalf_xz, qsnow_xz, alpha=0.6,levels=xz_snow_levels,cmap=cmap_blue)# 
    #plt.colorbar()
    plt.contourf(grid, heighthalf_xz, wind_xz, alpha=0.6, levels=xz_wind_levels,cmap=plt.cm.jet)
    #plt.contourf(grid, heighthalf_xz, theta, alpha=0.6, levels=xz_theta_levels,cmap=plt.cm.jet)
    cb = plt.colorbar(orientation="vertical")
    cb.set_label('Wind magnitude speed (m/s)')    

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
        cs = plt.contour(grid, heighthalf_xz, theta, theta_int, colors='black')
        plt.clabel(cs, inline=1,  fmt='%1.0f', fontsize=12, colors='grey')

    plt.plot(arange(Nx),heightground_x,color='black')
    plt.fill_between(arange(Nx),heightground_x,0,facecolor='lightgrey')
    plt.xticks(np.arange(0,Nx,20),np.round(longitude[Ny/2,::8], 1), fontsize='small')
    plt.yticks(np.arange(z_min,z_max,dz), fontsize='small') # JRL: edited in z_min
    plt.xlabel('Longitude')
    plt.ylabel(u'Height ASL [m]')
    plt.title('Potential temperature at ' + str(datetuple[3]) + 'Z on ' + str(datetuple[:3]))
    #plt.show()
    fname = outdir+naming+str(nest)+'_'+str(datetuple[2])+str(datetuple[3])+'Z_xsection.png'
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
    
