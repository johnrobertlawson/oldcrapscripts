## Adaptation of pywrf by John Lawson
import numpy as N
import os

import mapWRF
import xzWRF
import tzWRF
import skewT
import ts
from pywrfplotParams import outdir, naming 

### Options for plotting
domains_and_nests = 0  # Add ability to choose time for MSLP
xz_cross_sections = 1 # Add ability to plot top-down location of xz x-sec
xz_animation = 0 # Create animation of xz
tz_cross_sections = 0
skew_T = 0 # Add convert time thing
time_series = 0 # Times are currently 'hours since start'

# Select time for plots
timetuple = (2011,12,1,0,0,0)
# Do multiple plots
tlist = N.arange(0,24,1)
timelist = [(2011,12,1,t,0,0) for t in tlist] # need a more elegant way to choose multiple day ranges?
partdaylist = range(0,7)
timelist.extend([(2011,12,2,t,0,0) for t in partdaylist])

# Select locations for plots
loclist = ['loc'+str(n) for n in range(1,10)]
# Grid labels for domains
gl = (r'12\,km',r'4\,km',r'1.3\,km')
# Start lat/lon and end lat/lon for xsections
#startll = (40.9336, -111.8921) # Centerville
startll = (40.85, -112.3) # GSL
endll = (41.4298, -110.1205) # Church Butte
#startll = (41.759, -110.132) # North of Church Butte
#endll = (40.820, -110.108) # Near Kings Peak
#startll = (41.759, -110.963) 
#endll = (40.820, -110.963) # further west

# Dictionaries for cross-section:
XS = {}
XS['westeast'] = {'startll': (40.85, -112.3), 'endll': (41.4298, -110.1205),'vars':('W','parawind')} # Centerville to Church Butte
XS['northsouth'] = {'startll':(42, -110.55), 'endll':(40.5, -110.373),'vars':('W','perpwind')} # North of Church Butte to Uintahs

if domains_and_nests:
    # Send cross-section data to plot on domains?
    mapWRF.mapDomains(gl)
    #mapWRF.mapTerrain(includePMSL=True)
if xz_cross_sections:
    #os.system('rm '+outdir+'xz_timelist.txt') # This removes the file ready to write another
    #os.system('touch '+outdir+'xz_timelist.txt')
    # xzWRF needs time as a six-element tuple: 
    # (year,month,day,hour,min,sec)
    for t in timelist: # wind has: normal, parallel, magnitude options.
        xzWRF.xzCloudPlot(nest=3, datetuple=t, XS=XS)
if skew_T:
    skewT.skewTPlot(nest=3, timetuple=timetuple)
if tz_cross_sections:
    tzWRF.tzCloudPlot(nest=3) # JRL: need to fix WRF time step issue
if time_series:
    for l in loclist:
        ts.timeseries(nest=3,loc=l)
if xz_animation:
    print "Creating loop."
    os.system('convert -delay 50 -loop 0 @'+outdir+'xz_timelist.txt '+outdir+naming+'loop.gif')
    
