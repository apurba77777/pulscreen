
import os, sys
import numpy as np
import argparse as ap
from funlib.utils import *
from funlib.basicfns import *
from funlib.plotfns import *



#	-----------------------------------------------------------------------
#	Script for FRB polarization analysis
#
#								AB, September 2026
#	-----------------------------------------------------------------------

def get_args():
    
	#   Read command-line arguments for pulse generation

	parser = ap.ArgumentParser(
		description = "Arguments for pulse processing"
	)

	#   Input files, identifiers & parameters
	parser.add_argument("--inpkl", help = "pickle file containing the simulated pulse", type = str, default = None)
	parser.add_argument("--derm", help = "Derotating RM", type = float, default = 0.0)
	parser.add_argument('--rmrange', nargs=2, type=float, help="Time range to measure RM", default=None)	
	parser.add_argument("--rmbins", help = "Bins for RM measurement", type = int, default = 1)
	parser.add_argument("--rmax", help = "Maximum RM to search", type = float, default = 1.0e3)
	parser.add_argument('--nrange', nargs=2, type=float, help="Time range to measure noise", default=None)
	parser.add_argument('--prange', nargs=2, type=float, help="Time range to plot", default=None)

	parser.add_argument("--showfdf", help = "Show Faraday depth function", action='store_true')

	#	Actions
	parser.add_argument("--stksds", help = "Plot Stoke dynamic spectra", action='store_true')
	parser.add_argument("--polr", help = "Plot polarization properties", action='store_true')
	parser.add_argument("--dpa", help = "Plot PA profile", action='store_true')

	args = parser.parse_args()

	return args
#   --------------------------------------------------------------------------------------------------------------------------


#	--------------------------	Read inputs	-------------------------------

argus   = get_args()

#   Read and complain about input parameters
if (argus.inpkl == None):
    print(" Missing pickle! You might want to make one...")
    sys.exit()
else:
	dsfile	=	open(f"{argus.inpkl}.pkl",'rb')
	dsdata	=	pkl.load(dsfile)
	dsfile.close()
    


#	-------------------------	Analyze pulse	-------------------------------


nchan	=	len(dsdata.fmhzarr)
nmsrange=	(min(dsdata.tmsarr), max(dsdata.tmsarr))
if (argus.nrange!=None):
	nmsrange	= argus.nrange

rmrange=	(min(dsdata.tmsarr), max(dsdata.tmsarr))
if (argus.rmrange!=None):
	rmrange	= argus.rmrange

#	Estimate Noise spectra
noisespec	=	estimate_noise(dsdata.dspec4, dsdata.tmsarr, nmsrange[0], nmsrange[1]) 
noistks		=	np.sqrt(np.nansum(noisespec**2,axis=1))/len(dsdata.fmhzarr)

#	Estimate RM
rmarr		=	[]
rmwinms		=	(rmrange[1] - rmrange[0]) / argus.rmbins
for i in range(0, argus.rmbins):	
	msl			=	rmrange[0] + i*rmwinms
	msr			=	msl + rmwinms
	resrmt		=	estimate_rm(dsdata.dspec4, dsdata.fmhzarr, dsdata.tmsarr, noisespec, msl, msr, argus.rmax, 1.0, 0, 0, showplt=argus.showfdf)
	#print(msl, msr, resrmt[0], resrmt[1])
	rmarr.append([msl, msr, (msl+msr)/2, resrmt[0], resrmt[1]])

rmarr		=	np.array(rmarr)


#	Correct fot the given RM

corrdspec	=	unfarot(dsdata.dspec4, dsdata.fmhzarr, argus.derm)
tsdata		=	calc_profiles(corrdspec, dsdata.fmhzarr, dsdata.tmsarr, noisespec, 0, 0)

pltrange	=	[np.amin(dsdata.tmsarr),np.amax(dsdata.tmsarr)]
if (argus.prange!=None):
	pltrange	= argus.prange


if (argus.stksds):	
	plot_iquvt(None,corrdspec,tsdata.iquvt,dsdata.fmhzarr,dsdata.tmsarr,pltrange,[5.0,8.0])
	
if(argus.polr):
	plot_ilvpadst(None,noistks,corrdspec,tsdata,dsdata.fmhzarr,dsdata.tmsarr,pltrange,rmarr,[4.0,7.2])

if(argus.dpa):
	plot_dpa(None,noistks,tsdata,dsdata.tmsarr,[4.0,4.0],5)
        





































