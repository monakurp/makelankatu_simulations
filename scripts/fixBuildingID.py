#!/usr/bin/env python
import sys
import argparse
import numpy as np
from utilities import writeLog
from mapTools import *
from netcdfTools import *
import scipy.ndimage.measurements as snm


'''
Description:
checks if any grid points are missing a type if a building is not present. 
Also forces any incorrect surface types to follow PIDS.

Author: Jani Stromberg
        jani.stromberg@helsinki.fi
        University of Helsinki
'''

#==========================================================#
parser = argparse.ArgumentParser( prog='fixBuildingID.py', 
                                  description='''Defines a unique building_id for each building.''')
parser.add_argument("-f","--building_id", type=str, \
  help="Name of building id input raster data file.")
args = parser.parse_args()

#==========================================================#
# Renaming ... that's all
building_id = args.building_id

# Read in the data. The output raster will inherit Rdict properties not shown here.
Rdict = readNumpyZTile(building_id)
R = Rdict['R']
Rdims = np.array(np.shape(R))
ROrig = Rdict['GlobOrig']
dPx = Rdict['dPx']

# Give an individual id for each building
Rmod = np.copy( R )
Rmod[Rmod>0] = 1
Rmod[Rmod<=0] = 0 # snm.label functions with 0s and 1s
labeled_array, num_features = snm.label( Rmod )
if np.sum( labeled_array < 0 ) > 0:
  print('Negative building_id! Exiting.')
  sys.exit(1)
labeled_array[labeled_array==0] = -9999
R = labeled_array.astype(int)
del Rmod, labeled_array, num_features

# Save file
Rdict["R"]=R
Rdict["Labels"] = ['']
saveTileAsNumpyZ(building_id, Rdict)


print(" ...{} saved successfully.".format(building_id))
