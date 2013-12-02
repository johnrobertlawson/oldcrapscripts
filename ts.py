import matplotlib.pyplot as plt
import numpy as N
import pdb
import os
import cPickle as pickle
import sys

sys.path.append('/uufs/chpc.utah.edu/common/home/u0737349/lawsonpython/') 
from pywrfplotParams import *
from pywrfplotUtils import *
import generalmet

### USER OPTIONS
height, width = 11, 16
plt.rc('text',usetex=True)
fonts = {'family':'Computer Modern','size':16}
plt.rc('font',**fonts)  
N.set_printoptions(precision=3,suppress=True)

### LAZY SETTINGS
#00z 1st occurs at forecast time 42
fchr = 42

def timeseries(nest, loc): 
    # columns are:
    # col 1: forecast time in hr
    # col 3,4: nearest grid to stn
    # col 5: 2m Temp (K)
    # col 7: U wind
    # col 8: V
    # col 9: sfc P (Pa)
   
    ts_fname = directory+loc+'.d0'+str(nest)+'.TS'
    """
    txt_fname = directory+loc+'.d0'+str(nest)+'.txt'
    #pdb.set_trace()
    try:
        # print "Looking for .txt time series"
        testfile = open(txt_fname)
    except IOError:
        print "Converting .TS to .txt"
        os.system('cp ' + ts_fname +' ' + txt_fname)
    finally:
        print "Opening .txt time series"
        ts_hour,ix,iy,t,u,v,psfc = N.loadtxt(txt_fname,'r',usecols=(1,3,4,5,7,8,9),skiprows=1, delimiter='\t')
    """

    ts_hour,ix,iy,t,u,v,psfc = N.loadtxt(ts_fname,usecols=(1,3,4,5,7,8,9),skiprows=1,unpack=True)
    
    # Create wind magnitude
    #wind = N.sqrt(u**2 + v**2)
    wspd, wdir = generalmet.combine_wind_components(u,v)
    
    # Save to pickle for other plots
    picklepath = '/uufs/chpc.utah.edu/common/home/u0737349/dsws/thesis/timeseries/'

    savedict = {'wspd':wspd, 'wdir':wdir, 'ts_hour':ts_hour, 'locnum':loc, 'uintah':uintah}
    picklename = picklepath + loc + uintah + '_WRF.p'
    pickle.dump(savedict,open(picklename, 'wb'))

    # Need to create list of times in human format for plot
    plt.figure(figsize=(width,height))
    plt.plot(ts_hour,wspd)
    plt.plot([42,42],[0,10])
    plt.plot([54,54],[0,10])
    plt.plot([66,66],[0,10])
    #plt.plot(humantime,wind

    #plt.xticks(
    #plt.yticks(
    plt.xlabel('Time')
    plt.ylabel('Wind speed (m/s)')
    #plt.title('Meteogram of wind speed with time, WRF output')
    #plt.show()
    plt.savefig(outdir+naming+str(nest)+'_'+loc+'_timeseries.png')
    plt.close()


