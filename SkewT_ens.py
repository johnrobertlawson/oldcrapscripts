# Load in each ensemble member WRF output for a location
# Find mean and std of temp/dewp/wind speed/wind direction
# Plot composite Skew Ts of this data

import cPickle as pickle
import numpy as N
import sys
import pdb
import matplotlib as M
M.use('Agg')
import matplotlib.pyplot as plt

sys.path.append('/home/jrlawson/pythoncode/general/')   
import gridded_data
import skewT as ST
from Params import T_zero,T_base,kappa,barb_increments,P_bot
#from Utils import gamma_s,td,e,getDimensions,convert_time
import generalmet

# OPTIONS
# Time tuple
timetuple = (2006,05,23,12,0,0)
#timetuple = (2006,05,24,00,0,0)
# Time for filenames etc
nicetime = ''.join(["%02u" %n for n in timetuple[0:4]])

# Plotting options
dir = '/home/jrlawson/pythoncode/plotting/output/2006052300/'
height, width = (10,10)
plt.rc('text',usetex=True)  
fonts = {'family':'serif','size':16}
plt.rc('font',**fonts) 


# Here we go
ensnames = ['C0'] + ['E' + str(n) for n in range(1,11)]
vars = ('u','v','T','Td')

D = {}
for ens in ensnames:
    if ens == 'C0':
        print "Skipping control"
        continue
    else:
        print "Beginning ", ens
        
    pickledir = '/home/jrlawson/data/sounding/WRFoutput/'
    picklef = pickledir + 'WRFsounding_' + nicetime + '_' + ens + '.p'
    
    try:
        with open(picklef,'rb') as f:
            E_D = pickle.load(f)
    except IOError:
        print "No data there!"
    finally:
        for z in vars:
            D[z] = gridded_data.dstack_loop(E_D[z],D,z)
        if ens == 'E1':
            P = E_D['P']            

data = {'std':{},'mean':{}}
# Calculate values
for z in vars:
    data['mean'][z] = N.mean(D[z],axis=2)[0]
    data['std'][z] = N.std(D[z],axis=2)[0]

# Plot
# Top and bottom limit of Skew T plots
P_b = 105000.0
P_t = 20000.0


# Create figure
fig = plt.figure(figsize=(width,height))
ST._isotherms()
ST._isobars()
ST._dry_adiabats()
ST._moist_adiabats()

thin_locs = gridded_data.thinned_barbs(P) 

ST._windbarbs_real(data['mean']['u'],data['mean']['v'],P[thin_locs],delta=1,n=35,color='black')
ST._temperature_real(data['mean']['T'],P,color='black',linestyle='solid')
ST._dewpoint_real(data['mean']['Td'],P,color='black',linestyle='dashed')  

plt.axis([-20,50,P_b,P_t])
plt.xlabel(r'Temperature ($^{\circ}$C) at 1000\,hPa')
xticks = N.arange(-20,51,5)
plt.xticks(xticks,['' if tick%10!=0 else str(tick) for tick in xticks])
plt.ylabel('Pressure (hPa)')
yticks = N.arange(P_bot,P_t-1,-10**4)
ytix = ["%4u" %(p/100.0) for p in yticks]
plt.yticks(yticks,ytix)

fname = 'meanSkewT_'+nicetime+'.png'
plt.savefig(dir+fname,bbox_inches='tight')
plt.clf()

for z in vars:
    plt.figure()
    plt.gca().invert_yaxis()
    if (z=='T') or (z=='Td'):
        plt.plot(data['std'][z],P)
        plt.savefig(dir+'std_'+z+'_'+nicetime+'.png')
        plt.clf()
    elif z=='u':
        continue
    else:
        wspd,wdir = generalmet.combine_wind_components(data['std']['u'],data['std']['v'])
        plt.plot(wspd,P[thin_locs])
        plt.savefig(dir+'std_wspd_'+nicetime+'.png')
        plt.clf()
        plt.figure()
        plt.gca().invert_yaxis()
        plt.plot(wdir,P[thin_locs])
        plt.savefig(dir+'std_wdir_'+nicetime+'.png')
        
