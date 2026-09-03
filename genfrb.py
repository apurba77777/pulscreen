
import os, sys
import matplotlib as mpl
import numpy as np
import argparse as ap
import yaml as ym
from funlib.genfns import *
from funlib.utils import *
from funlib.plotfns import *
from funlib.basicfns import *


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
noisespec	=	estimate_noise(dynspec, tmsarr, min(tmsarr), max(tmsarr)) 
tsdata		=	calc_profiles(dynspec, fmhzarr, tmsarr, noisespec, 0, 0)
plot_iquvt(None, dynspec, tsdata.iquvt, fmhzarr, tmsarr, (min(tmsarr), max(tmsarr)), [5.0,8.0])






























































