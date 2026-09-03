
import os, sys
import matplotlib as mpl
import numpy as np
import argparse as ap
import yaml as ym
from funlib.genfns import *
from funlib.utils import *


#	-----------------------------------------------------------------------
#	Simulating a pulse profile
#
#								AB, Sep 2026
#                                                               
#	-----------------------------------------------------------------------

def gen_args():
    
    #   Read command-line arguments for pulse generation

    parser = ap.ArgumentParser(
        description = "Arguments for pulse generation"
    )

    #   Input files, identifiers & parameters
    parser.add_argument("--infile", help = "YAML file with input params", type = str, default = None)
    parser.add_argument("--pname", help = "Name of the generated pulse", type = str, default = "lightning")
        
    args = parser.parse_args()

    return args
#   --------------------------------------------------------------------------------------------------------------------------





#	--------------------------	Read inputs	-------------------------------

argus   = gen_args()

#   Read and complain about input parameters
if (argus.infile == None):
    print(" Missing input YAML file! Please provide one...")
    sys.exit()
else:
    with open(argus.infile+'.yml', 'r') as infl:
        pars = ym.load(infl, Loader=ym.SafeLoader)
        print(" Inputs provided -- \n")
        print(ym.dump(pars, sort_keys=False))


if (os.path.exists(pars['WorkDir'])):
    print("Found ",pars['WorkDir'])
else:
    print("Creating ",pars['WorkDir'])
    os.system("mkdir "+pars['WorkDir'])



#	-------------------------	Generate pulse	-------------------------------


fmhzarr		=	np.arange( pars['CenFreq'] - (pars['NChan']*pars['DfMhz'])/2.0 , pars['CenFreq'] + (pars['NChan']*pars['DfMhz'])/2.0, \
                                        pars['DfMhz'], dtype=float)

tmsarr		=	np.arange( -pars['TwinMs'], pars['TwinMs'], pars['TresMs'], dtype=float )

#	Generate inital dynamic spectrum with Gassian components
dynspec 	=	gauspuls(fmhzarr, tmsarr, pars['ComPars'])

for scrn in pars['ScrPars']:
    if (scrn['scrtype']=='fr'):
        #   Apply Faraday Rotation
        dynspec     =   frotateds(fmhzarr, dynspec, scrn['rm'])
    elif (scrn['scrtype']=='sc'):
        #	Scatter the dynamic spectrum 
        dynspec		=	scatterds(dynspec, fmhzarr, tmsarr, scrn['tms'], scrn['scind'])
    else:
        print("Unknown screen type. Doing nothing...")

dynspec     =   addnoise(dynspec, fmhzarr, tmsarr)



#   ------------------------    Pickle the pulse    ---------------------------


#	'Pickle' the simulated FRB and save it to the disk
fakefrb		=	simfrb(argus.pname, fmhzarr, tmsarr, pars['ScrPars'], pars['FrefMhz'], pars['ComPars'], dynspec)      

frbfile		=	open(f"{pars['WorkDir']}/{argus.pname}.pkl",'wb')             # Create the data directory, keep all simulated frbs 
pkl.dump(fakefrb, frbfile)		
frbfile.close()





#   Show the Pulse 
fig, axs = plt.subplots(5, figsize=(5, 10))
fig.suptitle('Scattered Dynamic Spectrum')

# Plot the mean across all frequency channels (axis 0)
axs[0].plot(np.nanmean(dynspec[0,:], axis=0), markersize=2 ,label='I')
axs[0].plot(np.nanmean(dynspec[1,:], axis=0), markersize=2, label='Q')
axs[0].plot(np.nanmean(dynspec[2,:], axis=0), markersize=2, label='U')
axs[0].plot(np.nanmean(dynspec[3,:], axis=0), markersize=2, label='V')

# Plot the 2D scattered dynamic spectrum
axs[1].imshow(dynspec[0], aspect='auto', interpolation='none', origin='lower', cmap='plasma')
axs[2].imshow(dynspec[1], aspect='auto', interpolation='none', origin='lower', cmap='plasma')
axs[3].imshow(dynspec[2], aspect='auto', interpolation='none', origin='lower', cmap='plasma')
axs[4].imshow(dynspec[3], aspect='auto', interpolation='none', origin='lower', cmap='plasma')
axs[4].set_xlabel("Time (samples)")
axs[4].set_ylabel("Frequency (MHz)")

plt.tight_layout()
plt.show()




























































