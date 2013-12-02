## Adaptation of pywrf by John Lawson
import numpy as N
import os

import mapWRF
import xzWRF_MCS
#import tzWRF
import skewT
#import ts
from Params import outdir, naming 

### Options for plotting
domains_and_nests = 0  # Add ability to choose time for MSLP
xz_cross_sections = 0 # Add ability to plot top-down location of xz x-sec
xz_animation = 0 # Create animation of xz
tz_cross_sections = 0
skew_T = 1 # Add convert time thing
time_series = 0 # Times are currently 'hours since start'

# Select time for plots
xs_time = (2006,5,24,5,30,0)
skewT_time = (2006,5,23,12,0,0)
# Do multiple plots
#tlist = N.arange(0,24,1)
#timelist = [(2011,12,1,t,0,0) for t in tlist] # need a more elegant way to choose multiple day ranges?

# Select locations for plots
loclist = ['loc'+str(n) for n in range(1,10)]
# Grid labels for domains
gl = (r'12\,km',r'4\,km',r'1.3\,km')
gl3 = (r'3\,km')

# Dictionaries for cross-section:
XS = {}
XS['zonalwind'] = {'startll': (40.81, -96.7), 'endll': (40.81, -93.65),'vars':('uwind','uwind_sr','W','pertpres','qvapor','theta'),'storm_motion':(14,270)} # Lincoln to S of DMX
# ^^ Storm motion as if it's wind ^^ 
if domains_and_nests:
    # Send cross-section data to plot on domains?
    mapWRF.mapDomains(gl3,terrain=False,drawxs=XS)
    #mapWRF.mapTerrain(includePMSL=True)
if xz_cross_sections:
    # xzWRF needs time as a six-element tuple: 
    # (year,month,day,hour,min,sec)
    #for t in xs_time: # wind has: normal, parallel, magnitude options.
    xzWRF_MCS.xzPlot(nest=1, datetuple=xs_time, XS=XS)
if skew_T:
    skewT.skewTPlot(nest=1, timetuple=skewT_time)
if tz_cross_sections:
    tzWRF.tzCloudPlot(nest=3) # JRL: need to fix WRF time step issue
if time_series:
    for l in loclist:
        ts.timeseries(nest=3,loc=l)
if xz_animation:
    print "Creating loop."
    os.system('convert -delay 50 -loop 0 @'+outdir+'xz_timelist.txt '+outdir+naming+'loop.gif')
    
