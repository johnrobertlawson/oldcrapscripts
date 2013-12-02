# Plot ensemble member diagnostics at a location
# e.g. CAPE, shear

# Imports
import numpy as N
import matplotlib as M
M.use('Agg')
import matplotlib.pyplot as plt
import pdb
from netCDF4 import Dataset
import sys
sys.path.append('/home/jrlawson/pythoncode/general/')

import wrf_tools
import generalmet
import gridded_data
from Utils import convert_time

# plotting
height, width = (10,10)
outdir = '/home/jrlawson/pythoncode/plotting/output/2006052300/'
plt.rc('text',usetex=True)
fonts = {'family':'serif','size':16}
plt.rc('font',**fonts)

# Options
CAPE = 0
shear_0_3 = 1
shear_0_6 = 1
T2 = 1
Q2 = 1

# Location lat/lon
loc = (41.31,-96.36)

# Domain
dom = 1

# Time
init_t_string = '2006052300'
timetuple = (2006,05,23,12,0,0)
nicetime = ''.join(["%02u" %n for n in timetuple[0:4]]) 

# Ensembles
ensembles = 'all'
if ensembles == 'all':
    ensnames = ['C0'] + ['E' + str(n) for n in range(1,11)]

# Observational data plotted too
obs = 1
obspath = '/home/jrlawson/data/sounding/oax_'+nicetime+'.txt' 

# Directory of netCDF files
ncdir = '/ptmp/jrlawson/predictability/'+init_t_string+'/wrfout/'
ncfname = 'wrfout_d0'+str(dom)+'_2006-05-23_00:00:00'

# Functions
def return_CAPE(ncpath,ens,loc,timetuple,lats,lons):
    nc = Dataset(ncpath,'r')


def return_shear(lev1,lev2,ncpath,ens,loc,timetuple,lats,lons,dom):
    nc = Dataset(ncpath,'r')
    time = convert_time(dom,timetuple)
    x,y,exactlat,exactlon = gridded_data.getXY(lats,lons,loc[0],loc[1])
    u = nc.variables['U'][time,:,y,x][:]
    v = nc.variables['V'][time,:,y,x][:]
    #wind = N.sqrt(u**2 + v**2)
    geopot = nc.variables['PH'][time,:,y,x][:] + nc.variables['PHB'][time,:,y,x][:]
    Z = geopot/9.81
    zAGL = Z - nc.variables['HGT'][time,y,x] # Height AGL
    
    # Compute the height, AGL, of the u/v surfaces
    zAGL_wind = N.array([(z0+z1)/2.0 for z0,z1 in zip(zAGL[:-1],zAGL[1:])])
    # Find wind at lev2 via interpolation
    lev2 *= 1000 # get into metres
    u2 = N.interp(lev2,zAGL_wind,u)
    v2 = N.interp(lev2,zAGL_wind,v)

    if lev1 == 0:
        u1 = nc.variables['U10'][time,y,x]
        v1 = nc.variables['V10'][time,y,x]
    else:
        # Find wind at lev1 via interpolation
        pass 
    shear = N.sqrt((u2-u1)**2 + (v2-v1)**2)
    return shear

def load_sounding_data(obspath):
    var_names = ['pres','hght','temp','dwpt','relh','mixr','drct','sknt','thta','thte','thtv']
    D = {}
    D = N.genfromtxt(obspath,unpack=True,dtype={'names':var_names,'formats':['f4']*len(var_names)},filling_values=-9999,)
    return D

def return_obs_shear(lev1,lev2,obspath):
    
    # Load data from text file
    D = load_sounding_data(obspath)
    # Interpolate wspd and wdir for lev2 from sounding
    # Wdir interpolation might be rubbish if it crosses north...
    wspd2 = N.interp(lev2,D['hght'],D['sknt']*0.514444)
    wdir2 = N.interp(lev2,D['hght'],D['drct'])
    u2,v2 = generalmet.decompose_wind(wspd2,wdir2)

    if lev1 == 0:
        u1,v1 = generalmet.decompose_wind(D['sknt'][0]*0.514444,D['drct'][0])
    else:
        pass
    # Compute shear
    shear = N.sqrt((u1-v1)**2 + (u2-v2)**2)
    return shear

def plot_shear(lev1,lev2):
    fig = plt.figure()
    for i,ens in enumerate(ensnames):
        ncpath = ncdir+ens+'/'+ncfname
        ens_shear = return_shear(lev1,lev2,ncpath,ens,loc,timetuple,lats,lons,dom)
        plt.bar(i,ens_shear)
    if obs == 1:
        # Include observed shear
        obs_shear = return_obs_shear(lev1,lev2,obspath)
        plt.bar(len(ensnames),obs_shear,color='black')
        xlabelnames = ensnames + ['OBS']
        xticklocs = [0.5+n for n in range(len(ensnames)+1)]
    else:
        xlabelnames = ensnames
        xticklocs = [0.5+n for n in range(len(ensnames))]
    fname = '_'.join(('shear',str(lev1),str(lev2),nicetime))
    plt.xticks(xticklocs,xlabelnames)
    plt.ylabel('Shear '+str(lev1)+'--'+str(lev2)+r'\,km (s$^{-1}$)')
    fig.savefig(outdir+fname, bbox_inches='tight', pad_inches=0.5)
    fig.clf()

def plot_Q2(loc,timetuple,lats,lons,dom):
    fig = plt.figure()
    for i,ens in enumerate(ensnames):
        # Load data from nc file
        time = convert_time(dom,timetuple)
        x,y,exactlat,exactlon = gridded_data.getXY(lats,lons,loc[0],loc[1])
        ncpath = ncdir+ens+'/'+ncfname
        nc = Dataset(ncpath,'r')
        ens_Q2 = nc.variables['Q2'][time,y,x]*1000
        plt.bar(i,ens_Q2)
    if obs == 1:
        # Include observed shear
        D = load_sounding_data(obspath)
        obs_Q2 = D['mixr'][0]
        plt.bar(len(ensnames),obs_Q2,color='black')
        xlabelnames = ensnames + ['OBS']
        xticklocs = [0.5+n for n in range(len(ensnames)+1)]
    else:
        xlabelnames = ensnames
        xticklocs = [0.5+n for n in range(len(ensnames))]
    fname = '_'.join(('Q2',nicetime))
    plt.xticks(xticklocs,xlabelnames)
    plt.ylabel('Mixing ratio at surface (g kg$^{-1}$)')
    fig.savefig(outdir+fname, bbox_inches='tight', pad_inches=0.5)
    fig.clf()

def plot_T2(loc,timetuple,lats,lons,dom):
    fig = plt.figure()
    for i,ens in enumerate(ensnames):
        # Load data from nc file
        time = convert_time(dom,timetuple)
        x,y,exactlat,exactlon = gridded_data.getXY(lats,lons,loc[0],loc[1])
        ncpath = ncdir+ens+'/'+ncfname
        nc = Dataset(ncpath,'r')
        ens_T2 = nc.variables['T2'][time,y,x] - 273.15
        plt.bar(i,ens_T2)
#    if obs == 1:
#        # Include observed shear
#        D = load_sounding_data(obspath)
#        obs_T2 = D['thta'][0]
#        plt.bar(len(ensnames),obs_T2,color='black')
#        xlabelnames = ensnames + ['OBS']
#        xticklocs = [0.5+n for n in range(len(ensnames)+1)]
#    else:
    xlabelnames = ensnames
    xticklocs = [0.5+n for n in range(len(ensnames))]
    fname = '_'.join(('T2',nicetime))
    plt.gca().autoscale(enable=True,axis='y',tight=True)
    plt.xticks(xticklocs,xlabelnames)
    plt.ylabel('2m potential temperature (C)')
    fig.savefig(outdir+fname, bbox_inches='tight', pad_inches=0.5)
    fig.clf()

  

# First, get lat and lon from WRF (using control)
nc = Dataset(ncdir+'C0/'+ncfname,'r')
lats,lons = wrf_tools.latlon_1D(nc)

if shear_0_3:
    plot_shear(0,3)

if shear_0_6:
    plot_shear(0,6)

if T2:
    plot_T2(loc,timetuple,lats,lons,dom)

if Q2:
    plot_Q2(loc,timetuple,lats,lons,dom)
