
import os, sys
import matplotlib as mpl
import numpy as np
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
from funlib.utils import *

mpl.rcParams['pdf.fonttype']	= 42
mpl.rcParams['ps.fonttype'] 	= 42
mpl.rcParams['savefig.dpi'] 	= 600
mpl.rcParams['font.family'] 	= 'sans-serif'
mpl.rcParams['font.size']		= 8
mpl.rcParams["xtick.major.size"]= 3
mpl.rcParams["ytick.major.size"]= 3

#   ------------------------------------------------------------------------------
#	Functions for simulating a pulse
#
#								AB, September 2026
#	-------------------------------------------------------------------------------


def gauspuls(fmhzarr, tmsarr, gcomps):
	
	#	Generate dynamic spectrum for Gaussian pulses
	
    gpdspec	=	np.zeros((4, fmhzarr.shape[0], tmsarr.shape[0]), dtype=float)
    ngp		=	len(gcomps)
    plsarr	=	np.zeros(len(tmsarr), dtype=float)

    print(f"Generating a {ngp} component pulse...\n")

    for g in range(0,ngp):

                nrmarr	=	gcomps[g]['peak']*( (fmhzarr/np.nanmedian(fmhzarr))**gcomps[g]['spec'] )
                plsarr	=	np.exp( -(tmsarr - gcomps[g]['cent'])**2 / (2*(gcomps[g]['wms']**2)) )
                pa_arr 	= 	np.deg2rad(gcomps[g]['padeg'] + (tmsarr - gcomps[g]['cent'])*gcomps[g]['dpadt'])
        
                for c in range(0,len(fmhzarr)):
                    gpdspec[0, c]	=	gpdspec[0, c] + nrmarr[c]*plsarr	                        		    # I
                    gpdspec[1, c]   = 	gpdspec[0, c] * gcomps[g]['lfrac'] * np.cos(2 * pa_arr)                 # Q
                    gpdspec[2, c]   =	gpdspec[0, c] * gcomps[g]['lfrac'] * np.sin(2 * pa_arr)                 # U
                    gpdspec[3, c]   =   gpdspec[0, c] * gcomps[g]['vfrac']    		                            # V

    return (gpdspec)
#	--------------------------------------------------------------------------------


def frotateds(fmhzarr, dspec, rm):
	
	#	Apply Faraday oration to Stokes dynamic spectrum 
	
    lm2arr	=	(ccC*1.0e-8 / fmhzarr)**2
    lm20 	= 	np.nanmedian(lm2arr)
        
    for c in range(0,len(fmhzarr)):
        pa_farr 	= 	rm*(lm2arr[c] - lm20)						# Faraday rotation
        dspec[1, c] = 	dspec[1, c] * np.cos(2 * pa_farr)           # Q
        dspec[2, c] =	dspec[1, c] * np.sin(2 * pa_farr)           # U

    return (dspec)
#	--------------------------------------------------------------------------------


def scatterds(dspec, fmhzarr, tmsarr, taums, scindex):
        
    # Scatter a given dynamic spectrum
	
    scdspec = np.zeros(dspec.shape, dtype=float)
    taucms 	= taums * ((fmhzarr / np.nanmedian(fmhzarr)) ** scindex)

    for c in range(len(fmhzarr)):
        irfarr = np.heaviside(tmsarr, 1.0) * np.exp(-tmsarr / taucms[c]) / taucms[c]
        for stk in range(4): 
            scdspec[stk, c] = np.convolve(dspec[stk, c], irfarr, mode='same')

    print(f"--- Scattering time scale = {taums:.2f} ms, {np.nanmin(taucms):.2f} ms to {np.nanmax(taucms):.2f} ms")

    return (scdspec)
#	--------------------------------------------------------------------------------


def addnoise(dspec, fmhzarr, tmsarr):
        
    # Add noise to a dynamic spectrum
	
    for stk in range(4): 
        dspec[stk] = dspec[stk] + np.random.normal(loc=0.0, scale=1.0, size=(fmhzarr.shape[0], tmsarr.shape[0]))

    return (dspec)
#	--------------------------------------------------------------------------------













































































