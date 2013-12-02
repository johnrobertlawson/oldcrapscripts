# Run this script with and without uintah by changing pywrfplotParams variable 'uintah'.

### IMPORTS
import numpy as N
import math
import matplotlib as M
M.use('Agg')
import matplotlib.pyplot as plt
import pdb
import cPickle as pickle

# User imports
import sys

import skewT as ST
from Params import T_zero,T_base,kappa,barb_increments,P_bot,outdir,ens
from Utils import gamma_s,td,e,openWRF,getDimensions,convert_time

sys.path.append('/home/jrlawson/pythoncode/general/')
import unix_tools
import generalmet
import gridded_data

# OPTIONS
# plotting
height, width = (10,10)
plt.rc('text',usetex=True)  
fonts = {'family':'serif','size':16}
plt.rc('font',**fonts) 

# Time tuple
timetuple = (2006,05,23,12,0,0)
#timetuple = (2006,05,24,00,0,0)
# Time for filenames etc
nicetime = ''.join(["%02u" %n for n in timetuple[0:4]]) 

# File with verification sounding data
obsfile = '/home/jrlawson/data/sounding/oax_'+nicetime+'.txt'
# Lat/long to compare model soundings are:
obs_latlon = (41.31,-96.36)
# Pick domain
dom = 1

# Generate pickle files of WRF output soundings
save_wrf_output = 1 
# Pickle files are saved here:
pickledir = '/home/jrlawson/data/sounding/WRFoutput/'
picklef = pickledir + 'WRFsounding_' + nicetime + '_' + ens + '.p'
unix_tools.createfolder(pickledir)

outdirfull = outdir + 'soundings/'
unix_tools.createfolder(outdirfull)

### FUNCTIONS
def gettime():
    t = convert_time(dom,timetuple)
    return t

def dewpoint(T,RH): # Lawrence 2005 BAMS?
    #T in C
    #RH in 0-100 format
    es = 6.11 * (10**((7.5*T)/(237.7+T)))
    e = es * RH/100.0
    alog = 0.43429*N.log(e) - 0.43429*N.log(6.11)
    Td = (237.7 * alog)/(7.5-alog) 
    #pdb.set_trace()
    return Td

var_names = ['pres','hght','temp','dwpt','relh','mixr','drct','sknt','thta','thte','thtv']
D = {}
D['obs'] = N.genfromtxt(obsfile,unpack=True,dtype={'names':var_names,'formats':['f4']*len(var_names)},filling_values=-9999,)

# Top and bottom limit of Skew T plots
P_b = 105000.0
P_t = 20000.0

# Time iterable
#hrIter = [x for x in xrange(0,24)] # Dec 1 hours
#dayIter = [1 for x in xrange(0,24)] # Dec the first
#hrIter.extend([x for x in xrange(0,7)]) # Dec 2
#dayIter.extend([2 for x in xrange(0,7)]) # Dec the second
#hrStr = ["%02d" %t for t in hrIter]
#dayStr = ["%02d" %t for t in dayIter]

# Create figure
fig = plt.figure(figsize=(width,height))
ST._isotherms()
ST._isobars()
ST._dry_adiabats()
ST._moist_adiabats()

# Observational data
pres = D['obs']['pres']*100 #hPa
temp = D['obs']['temp']
dewp = D['obs']['dwpt']
wspd = generalmet.convert_kt2ms(D['obs']['sknt'])
uwind,vwind = generalmet.decompose_wind(wspd,D['obs']['drct'])

"""
# Thin out 5000m and less
thinme = N.where(D['obs']['obs0'] < 5000)[0][::5]
leaveme = N.where(D['obs']['obs0'] > 5000)
uwindlow = D['obs']['obs13'][thinme]
vwindlow = D['obs']['obs14'][thinme]
plow = pres[thinme]
uwindhigh = D['obs']['obs13'][leaveme]
vwindhigh = D['obs']['obs14'][leaveme]
phigh = pres[leaveme]
uwind = N.hstack((uwindlow, uwindhigh))
vwind = N.hstack((vwindlow, vwindhigh))
pthin = N.hstack((plow,phigh))
"""
thin_locs = gridded_data.thinned_barbs(pres)
# replace pres with pthin below to have thinned values
ST._windbarbs_real(uwind[thin_locs],vwind[thin_locs],pres[thin_locs],delta=1,n=35,color='black')
ST._temperature_real(temp,pres,color='black',linestyle='solid')
ST._dewpoint_real(dewp,pres,color='black',linestyle='dashed')  
del pres, temp, dewp, thin_locs

# Plot
plt.axis([-20,50,P_b,P_t])
plt.xlabel(r'Temperature ($^{\circ}$C) at 1000\,hPa')
xticks = N.arange(-20,51,5)
plt.xticks(xticks,['' if tick%10!=0 else str(tick) for tick in xticks])
plt.ylabel('Pressure (hPa)')
yticks = N.arange(P_bot,P_t-1,-10**4)
ytix = ["%4u" %(p/100.0) for p in yticks]
plt.yticks(yticks,ytix)

#yticks = N.arange(P_bot,P_t-1,-10**4)
#plt.yticks(yticks,yticks/100)

# Legend
ax = plt.gca()
handles,labels = ax.get_legend_handles_labels()
obsArt = plt.Line2D((0,1),(0,0), color='black')
wrfArt = plt.Line2D((0,1),(0,0), color='blue')
plt.legend([obsArt,wrfArt],['Observed','WRF Output'])

#plt.savefig(outdir+'kslcSkewT.png',bbox_inches='tight')

### ADD MODEL DATA
nc = openWRF(dom) # Pick domain here
Nx,Ny,Nz,lons,lats,dx,dy,x,y = getDimensions(nc)

del x,y # We want our own locations

time = gettime()
x,y, exactlat, exactlon = gridded_data.getXY(lats[:,Nx/2],lons[Ny/2,:],obs_latlon[0],obs_latlon[1])
P = nc.variables['P'][time,:,y,x] + nc.variables['PB'][time,:,y,x]
elev = nc.variables['HGT'][0,y,x]

thin_locs = gridded_data.thinned_barbs(P)

ST._windbarbs(nc,time,y,x,P,thin_locs,n=45,color='blue')
ST._temperature(nc,time,y,x,P,linestyle='solid',color='blue')
ST._dewpoint(nc,time,y,x,P,linestyle='dashed',color='blue')

plt.xticks(xticks,['' if tick%10!=0 else str(tick) for tick in xticks])
plt.yticks(yticks,ytix)

plt.savefig(outdirfull+'verifSkewT_'+ens+'_'+nicetime+'.png',bbox_inches='tight',pad_inches=0.1)
plt.clf()

if save_wrf_output == 1:
    u,v = ST.return_data('wind',nc,time,y,x,thin_locs)
    T = ST.return_data('temp',nc,time,y,x,thin_locs,P=P)
    Td = ST.return_data('dwpt',nc,time,y,x,thin_locs,P=P)
    
    WRFdict = {'u':u,'v':v,'T':T,'Td':Td,'P':P}
    with open(picklef,'wb') as p:
        pickle.dump(WRFdict,p)
