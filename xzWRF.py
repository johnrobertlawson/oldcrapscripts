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
import sys
import cPickle as pickle

sys.path.append('/uufs/chpc.utah.edu/common/home/u0737349/lawsonpython/')  

from pywrfplotParams import *
from pywrfplotUtils import *
import customcmaps

### USER OPTIONS for figures
plt.rc('text',usetex=True)
fonts = {'family':'Computer Modern','size':16}
plt.rc('font',**fonts)
height, width = (9,17)

# Plot x-z cross-sections
# XS is dictionary with variables and cross-section paths to loop over.

def xzCloudPlot(nest,datetuple,XS,saveoutput=0):
    nc = openWRF(nest)
   
    time = convert_time(nest,datetuple) # JRL function to convert time stamp to time
    Nx,Ny,Nz,longitude,_lats,_dx,_dy,x_nr,y_nr = getDimensions(nc)
    lats = _lats.astype(N.float32)
    lons = longitude.astype(N.float32)
    
    for xs in XS.keys():
        
        # x-y position of start and end points
        stxy = getXYxs(lons[Ny/2,:],lats[:,Nx/2],XS[xs]['startll'])
        endxy = getXYxs(lons[Ny/2,:],lats[:,Nx/2],XS[xs]['endll'])
        
        # Joe Kington method for slice through data
        hyppts = int(N.hypot(endxy[0]-stxy[0], endxy[1]-stxy[1])) # Number of pts along hypotenuse
        xx = N.linspace(stxy[0],endxy[0],hyppts)
        yy = N.linspace(stxy[1],endxy[1],hyppts)
        xint, yint = xx.astype(int), yy.astype(int)
        angle = N.arctan((yy[-1]-yy[0])/(xx[-1]-xx[0])) # In radians

        # Get terrain heights?
        heightground_x,heighthalf_xz = _getHeightJRL(nc,nest, time, yy, xx, Nz, hyppts)    
        
        # Get potential temperature for contouring
        theta = nc.variables['T'][time,:,yint,xint] + T_base 

        # Set up plot
        # Length of hypotenus in km
        hlen = (1/1000.0) * N.sqrt((-1.0*hyppts*nc.DX*N.cos(angle))**2 + (hyppts*nc.DY*N.sin(angle))**2)
        # Hypotenuse grid length ticks: locs/labels along cross-section
        hgl_ticks = N.arange(0,hlen,hlen/hyppts)
        hgl_labels = ["%3.0f" %x for x in hgl_ticks]
        grid = N.swapaxes(N.repeat(N.array(hgl_ticks).reshape(hyppts,1),Nz,axis=1),0,1)
           
        # Load custom colorbar
        jet_sw = M.colors.LinearSegmentedColormap('jet_sw',customcmaps.jet_sw,256)
        w_jet_sw = M.colors.LinearSegmentedColormap('w_jet_sw',customcmaps.w_jet_sw,256)
        Wnorm = M.colors.BoundaryNorm(w_wind_levels,256)
        Unorm = M.colors.BoundaryNorm(u_wind_levels,256)

        # PLOTTING 
        for v in XS[xs]['vars']:
            plt.figure(figsize=(width, height))
            if nc.DX == nc.DY: # Make sure boxes are square...
                plt.axis([0,(hyppts*nc.DX)-1,z_min,z_max]) 
      
            # This is a very lazy/bad way to deal with the staggering
            u_xz = nc.variables['U'][time,:,yint,xint]
            v_xz = nc.variables['V'][time,:,yint,xint]
            
            if (v == 'parawind'):
                cfdata = (N.cos(angle)*u_xz) - (N.sin(angle)*v_xz)
                cf = plt.contourf(grid, heighthalf_xz, cfdata, alpha=0.6, levels=u_wind_levels,cmap=jet_sw,norm=Unorm,extend='both')
                #varCBlabel = 'Rotated zonal wind speed (ms${-1}$)'
                varCBlabel = 'Wind speed (ms$^{-1}$)'
            elif (v == 'perpwind'): # Currently reversed to give negative wind for easterly.
                cfdata = (N.cos(angle)*v_xz) + (N.sin(angle)*u_xz)
                cf = plt.contourf(grid, heighthalf_xz, cfdata*-1, alpha=0.6, levels=u_wind_levels,cmap=jet_sw,norm=Unorm,extend='both')
                #varCBlabel = 'Rotated perpendicular wind speed (ms${-1}$)'
                varCBlabel = 'Wind speed (ms$^{-1}$)'
            elif (v == 'windmag'):
                cfdata = N.sqrt(u_xz**2 + v_xz**2)
                cf = plt.contourf(grid, heighthalf_xz, cfdata, alpha=0.6, levels=xz_wind_levels,cmap=jet_sw,norm=Unorm)
            elif (v == 'uwind'):
                pass
            elif (v == 'vwind'):
                pass
            elif (v == 'W'):
                cfdata = nc.variables['W'][time,:-1,yint,xint]
                cf = plt.contourf(grid, heighthalf_xz, cfdata, alpha=0.6, levels=w_wind_levels,cmap=w_jet_sw,norm=Wnorm,extend='both')
                varCBlabel = 'Vertical wind speed (ms${-1}$)'
            elif (v == 'T'):
                pass
            elif (v == 'RH'):
                pass

            # Plot theta
            theta_int = N.arange(260.0,400.0,2.0)
            cs = plt.contour(grid, heighthalf_xz, theta, theta_int, colors='black')
            #manlocs = [(120.0,y) for y in N.arange(2200.0,6000.0,300.0)]
            #plt.clabel(cs, inline=1,  fmt='%3u', fontsize=12, colors='black',manual=manlocs)
            plt.clabel(cs, inline=1,  fmt='%3.0f', fontsize=12, colors='black')

            plt.plot(hgl_ticks,heightground_x,color='black')
            plt.fill_between(hgl_ticks,heightground_x,0,facecolor='lightgrey')
            
            labeldelta = 15
            #plt.xticks(N.arange(0,hyppts,labeldelta),hgl_labels[::labeldelta])
            plt.yticks(N.arange(z_min,z_max,dz)) # JRL: edited in z_min
            plt.xlabel('Distance along cross-section (km)')
            plt.ylabel(u'Height above sea level (m)')
            
            fname = '_'.join(('xs',xs,v,str(datetuple[1]),str(datetuple[2]),str(datetuple[3]))) + '.png'
            plt.savefig(outdir+fname,bbox_inches='tight',pad_inches=0.5)
            plt.close()
            
            # Draw colorbar the first time through
            try:
                with open(outdir+v+'CB.png'): pass
            except IOError:
                # Plot colorbar
                fig = plt.figure(figsize=(8,6))
                cbar_ax = fig.add_axes([0.15,0.05,0.7,0.02])
                cb = plt.colorbar(cf,cax=cbar_ax,orientation='horizontal')
                cb.set_label(varCBlabel)
                plt.savefig(outdir+v+'CB.png',bbox_inches='tight')
            
            # Save all data to pickle file
            if saveoutput:
                D = {'day':datetuple[2], 'hour':datetuple[3], 'theta':theta,
                    'grid':grid, 'heighthalf_xz':heighthalf_xz, 'hgl_ticks':hgl_ticks,
                    'heightground_x':heightground_x, 'hyppts':hyppts, 'hgl_labels':hgl_labels,
                    'xs':xs, 'v':v, 'data':cfdata, 'DX':nc.DX}
                pickledir = '/uufs/chpc.utah.edu/common/home/u0737349/dsws/thesis/xspickle/' + ens + '/'
                fname = '_'.join(('data',xs,v,str(datetuple[1]),str(datetuple[2]),str(datetuple[3]))) + '.p'
                with open(pickledir+fname,'wb') as p:
                    pickle.dump(D,p)
        
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
