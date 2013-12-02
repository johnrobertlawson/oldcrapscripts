# -*- coding:utf-8 -*- 
#Adapted from code by Geir Arne Waagbø
#By John Lawson (ISU)
#http://code.google.com/p/pywrfplot/

import matplotlib as M
import matplotlib.pyplot as plt
from numpy import arange
import pdb
import numpy as N
import math
import sys
import os
import cPickle as pickle

sys.path.append('/home/jrlawson/pythoncode/general/')  

from Params import *
from Utils import *
import customcmaps
import generalmet

### USER OPTIONS for figures
plt.rc('text',usetex=False)
#fonts = {'family':'Computer Modern','size':16}
fonts = {'family':'sans-serif','sans-serif':'Helvetica'}
plt.rc('font',**fonts)
height, width = (9,17)

# Plot x-z cross-sections
# XS is dictionary with variables and cross-section paths to loop over.

def xzPlot(nest,datetuple,XS,saveoutput=0):
    outdirxs = outdir + 'xs/'
    
    nc = openWRF(nest)
    time = convert_time(nest,datetuple) # JRL function to convert time stamp to time
    Nx,Ny,Nz,longitude,_lats,_dx,_dy,x_nr,y_nr = getDimensions(nc)
    lats = _lats.astype(N.float32)
    lons = longitude.astype(N.float32)
    
    # Create folder for plot if necessary
    try:
        os.stat(outdirxs)
    except:
        os.makedirs(outdirxs)
    
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
        heightground_x,heighthalf_xz = _getHeight(nc,nest, time, yy, xx, Nz, hyppts)    
        
        # Set up plot
        # Length of hypotenus in km
        hlen = (1/1000.0) * N.sqrt((-1.0*hyppts*nc.DX*N.cos(angle))**2 + (hyppts*nc.DY*N.sin(angle))**2)
        # Hypotenuse grid length ticks: locs/labels along cross-section
        hgl_ticks = N.arange(0,hlen,hlen/hyppts)
        hgl_labels = [r"%3.0f" %x for x in hgl_ticks]
        grid = N.swapaxes(N.repeat(N.array(hgl_ticks).reshape(hyppts,1),Nz,axis=1),0,1)

        # PLOTTING 
        for v in XS[xs]['vars']:
            plt.figure(figsize=(width, height))
            if nc.DX == nc.DY: # Make sure boxes are square...
                plt.axis([0,(hyppts*nc.DX/1000.0)-1,z_min,z_max+dz]) 
            else:
                print r"Domain not square - fix me" 
                break

            if (v == 'parawind'):
                u_xz = nc.variables['U'][time,:,yint,xint]
                v_xz = nc.variables['V'][time,:,yint,xint]
                cfdata = (N.cos(angle)*u_xz) - (N.sin(angle)*v_xz)
                cf = plt.contourf(grid, heighthalf_xz, cfdata, alpha=0.6, levels=u_wind_levels,cmap=jet_sw,norm=Unorm,extend='both')
                #varCBlabel = 'Rotated zonal wind speed (ms${-1}$)'
                varCBlabel = r'Wind speed (ms$^{-1}$)'
            if (v == 'perpwind'): # Currently reversed to give negative wind for easterly.
                u_xz = nc.variables['U'][time,:,yint,xint]
                v_xz = nc.variables['V'][time,:,yint,xint]
                cfdata = (N.cos(angle)*v_xz) + (N.sin(angle)*u_xz)
                cf = plt.contourf(grid, heighthalf_xz, cfdata*-1, alpha=0.6, levels=u_wind_levels,cmap=jet_sw,norm=Unorm,extend='both')
                #varCBlabel = 'Rotated perpendicular wind speed (ms${-1}$)'
                varCBlabel = r'Wind speed (ms$^{-1}$)'
            if (v == 'windmag'):
                u_xz = nc.variables['U'][time,:,yint,xint]
                v_xz = nc.variables['V'][time,:,yint,xint]
                cfdata = N.sqrt(u_xz**2 + v_xz**2)
                cf = plt.contourf(grid, heighthalf_xz, cfdata, alpha=0.6, levels=xz_wind_levels,cmap=jet_sw,norm=Unorm)
            if (v == 'uwind'):
                u_xz = nc.variables['U'][time,:,yint,xint]
                jet_sw = M.colors.LinearSegmentedColormap('jet_sw',customcmaps.jet_sw,256)
                Unorm = M.colors.BoundaryNorm(u_wind_levels,256)
                cf = plt.contourf(grid, heighthalf_xz, u_xz, alpha=0.6, levels=u_wind_levels,cmap=jet_sw,norm=Unorm,extend='both')
                varCBlabel = r'Wind speed (ms${-1}$)'
            if (v == 'uwind_sr'):
                u_xz = nc.variables['U'][time,:,yint,xint]
                storm_spd,storm_dir = XS[xs]['storm_motion'] # Get u-component from this
                zonal_storm_motion, dummy = generalmet.decompose_wind(storm_spd,storm_dir)
                u_sr = u_xz - zonal_storm_motion
                jet_sw = M.colors.LinearSegmentedColormap('jet_sw',customcmaps.jet_sw,256)
                Unorm = M.colors.BoundaryNorm(u_wind_levels,256)
                cf = plt.contourf(grid, heighthalf_xz, u_sr, alpha=0.6, levels=u_wind_levels,cmap=jet_sw,norm=Unorm,extend='both')
                varCBlabel = r'Storm-relative wind speed (ms${-1}$)'
            if (v == 'vwind'):
                v_xz = nc.variables['V'][time,:,yint,xint]
            if (v == 'W'):
                cfdata = nc.variables['W'][time,:-1,yint,xint]
                Wnorm = M.colors.BoundaryNorm(w_wind_levels,256)
                w_jet_sw = M.colors.LinearSegmentedColormap('w_jet_sw',customcmaps.w_jet_sw,256)
                cf = plt.contourf(grid, heighthalf_xz, cfdata, alpha=0.6, levels=w_wind_levels,cmap=w_jet_sw,norm=Wnorm,extend='both')
                varCBlabel = r'Vertical wind speed (ms${-1}$)'
            if (v == 'T'):
                pass
                #cfdata = nc.variables['T'][time,:-1,yint,xint]
            if (v == 'RH'):
                pass
                #cfdata = nc.variables['W'][time,:-1,yint,xint]
            if (v == 'theta'):
                theta = nc.variables['T'][time,:,yint,xint] + T_base 
                theta_int = N.arange(260.0,400.0,2.0)
                cs = plt.contour(grid, heighthalf_xz, theta, theta_int, colors='black')
                plt.clabel(cs, inline=1,  fmt='%3.0f', fontsize=12, colors='black')
            if (v == 'pertpres'):
                pertpres = nc.variables['P'][time,:,yint,xint]
                #cf = plt.contourf(grid, heighthalf_xz, cfdata, alpha=0.6, levels=w_wind_levels,cmap=w_jet_sw,norm=Wnorm,extend='both') 
                cf = plt.contourf(grid, heighthalf_xz, pertpres) 
                varCBlabel = r'Perturbation pressure (Pa)'

            if (v == 'qvapor'):
                qvapor = nc.variables['QVAPOR'][time,:,yint,xint]
                #cf = plt.contourf(grid, heighthalf_xz, cfdata, alpha=0.6, levels=w_wind_levels,cmap=w_jet_sw,norm=Wnorm,extend='both') 
                cf = plt.contourf(grid, heighthalf_xz, qvapor) 
                varCBlabel = r'Water vapor mixing ratio (kg kg$^{-1}$)'

            plt.plot(hgl_ticks,heightground_x,color='black')
            plt.fill_between(hgl_ticks,heightground_x,0,facecolor='lightgrey')
            
            labeldelta = 15
            #plt.xticks(N.arange(0,hyppts,labeldelta),hgl_labels[::labeldelta])
            plt.yticks(N.arange(z_min,z_max+dz,dz)) # JRL: edited in z_min
            plt.xlabel(r'Distance along cross-section (km)')
            plt.ylabel(r'Height above sea level (m)')
            
            fname = '_'.join(('xs',xs,v,"%02d" %datetuple[1],"%02d" %datetuple[2],"%02d" %datetuple[3],"%02d" %datetuple[4])) + '.png'
            plt.savefig(outdirxs+fname,bbox_inches='tight',pad_inches=0.5)
            plt.close()
            
            # Draw colorbar the first time through
            try:
                with open(outdirxs+v+'CB.png'): pass
            except IOError:
                # Plot colorbar
                fig = plt.figure(figsize=(8,6))
                cbar_ax = fig.add_axes([0.15,0.05,0.7,0.02])
                cb = plt.colorbar(cf,cax=cbar_ax,orientation='horizontal')
                cb.set_label(varCBlabel)
                plt.savefig(outdirxs+v+'CB.png',bbox_inches='tight')
            
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
        
def _getHeight(nc,nest,time,xx,yy,Nz,hypotenuse):    
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
