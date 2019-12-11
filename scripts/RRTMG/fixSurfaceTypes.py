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
parser = argparse.ArgumentParser( prog='fixSurfaceTypes.py', 
                                  description='''Fills missing types for PALM from input raster data.''')
parser.add_argument("-f0","--vegetation", type=str, \
  help="Name of vegetation type input raster data file.")
parser.add_argument("-f1","--pavement", type=str, \
  help="Name of pavement type input raster data file.")
parser.add_argument("-f2","--water", type=str, \
  help="Name of water type input raster data file.")
parser.add_argument("-f3","--buildingmask", type=str, \
  help="Name of building mask input raster data file.")
parser.add_argument("-f4","--building", type=str, \
  help="Name of building type input raster data file.")
parser.add_argument("-f5","--soil", type=str, \
  help="Name of soil type input raster data file.")
parser.add_argument("-f6","--building_id", type=str, \
  help="Name of building id input raster data file.")
args = parser.parse_args()
writeLog(parser, args)

#==========================================================#
# Renaming ... that's all
vegetation = args.vegetation
pavement = args.pavement
water = args.water
buildingmask = args.buildingmask
building = args.building
soil = args.soil
building_id = args.building_id

# Read in the data. The output raster will inherit Rdict1 properties not shown here.
Rdict1 = readNumpyZTile(vegetation)
R1 = Rdict1['R']
R1dims = np.array(np.shape(R1))
R1Orig = Rdict1['GlobOrig']
dPx1 = Rdict1['dPx']
nPx = np.shape(R1)
dPx = dPx1

Rdict2 = readNumpyZTile(pavement)
R2 = Rdict2['R']
R2dims = np.array(np.shape(R2))
R2Orig = Rdict2['GlobOrig']
dPx2 = Rdict2['dPx']
Rdict2 = None


Rdict3 = readNumpyZTile(water)
R3 = Rdict3['R']
R3dims = np.array(np.shape(R3))
R3Orig = Rdict3['GlobOrig']
dPx3 = Rdict3['dPx']
Rdict3 = None

Rdict4 = readNumpyZTile(buildingmask)
R4 = Rdict4['R']
R4dims = np.array(np.shape(R4))
R4Orig = Rdict4['GlobOrig']
dPx4 = Rdict4['dPx']
Rdict4 = None

Rdict5 = readNumpyZTile(building)
R5 = Rdict5['R']
R5dims = np.array(np.shape(R5))
R5Orig = Rdict5['GlobOrig']
dPx5 = Rdict5['dPx']
Rdict5 = None 

Rdict6 = readNumpyZTile(soil)
R6 = Rdict6['R']
R6dims = np.array(np.shape(R6))
R6Orig = Rdict6['GlobOrig']
dPx6 = Rdict6['dPx']
Rdict6 = None

Rdict7 = readNumpyZTile(building_id)
R7 = Rdict7['R']
R7dims = np.array(np.shape(R7))
R7Orig = Rdict7['GlobOrig']
dPx7 = Rdict7['dPx']
Rdict7 = None



if( (R1Orig == R2Orig).all() and (R2Orig == R3Orig).all() and (R3Orig == R4Orig).all() and 
    (R4Orig == R5Orig).all() and (R5Orig == R6Orig).all() and (R6Orig == R7Orig).all()):
  print(' Excellent! The origos match.')
else:
  print(' The tiles do not have identical origos. Exiting.')
  sys.exit(1)

if( (R1dims == R2dims).all() and (R2dims == R3dims).all() and (R3dims == R4dims).all() and 
    (R4dims == R5dims).all() and (R5dims == R6dims).all() and (R6dims == R7dims).all()):
  print(' Excellent! The dimensions match.')
else:
  print(' The tiles do not have identical dimensions. Exiting.')
  sys.exit(1)



print("Filling missing types ...")



# Prioritize building height, all other types (except soil) changed to fill values if there's overlap.
idx = R4 > 0
R1[idx] = -127; R2[idx] =-127; R3[idx] = -127

# Edit building type to fit the more accurate building height layer.
idx = R4 < 1
R5[idx] = -127

# Assume any holes are due to differences between building height and other layers, so handling these as paved surfaces.
idx = (R1 < 1) & (R2 < 1) & (R3 < 1) & (R4 < 1)
R2[idx] = 1

# Give an individual id for each building
Rmod = np.copy( R7 )
Rmod[Rmod>0] = 1
Rmod[Rmod<=0] = 0 # snm.label functions with 0s and 1s
labeled_array, num_features = snm.label( Rmod )
if np.sum( labeled_array < 0 ) > 0:
  print('Negative building_id! Exiting.')
  sys.exit(1)
labeled_array[labeled_array==0] = -127
R7 = labeled_array.astype(int)
del Rmod, labeled_array, num_features

# PIDS classification-----------------------------
# User can make their own changes to these conditions as they see fit.
# You can see original classification in the layer 'Labels'. i.e. Rdict['Labels']
# You can check the PIDS classes at: https://palm.muk.uni-hannover.de/trac/wiki/doc/app/iofiles/pids#no1

for i in range(nPx[0]):
    for j in range(nPx[1]):
        vegetype = R1[i,j]
        pavetype = R2[i,j]
        watertype = R3[i,j]
        buildtype = R5[i,j]
        soiltype = R6[i,j]
        buildid = R7[i,j]

        if (pavetype == 1):
            R2[i,j] = 2
        elif (pavetype == 2):
            R2[i,j] = 3
        elif (pavetype == 3):
            R2[i,j] = 12
        elif (pavetype == 4 or pavetype == 6):
            R2[i,j] = 10
        elif (pavetype == 5):
            R2[i,j] = 14
        elif (pavetype == 7 or pavetype == 12):
            R2[i,j] = 4
        elif (pavetype == 8 or pavetype == 15):
            R2[i,j] = 13
        elif (pavetype == 9):
            R2[i,j] = 7
        elif (pavetype == 10):
            R2[i,j] = 9
        elif (pavetype == 11):
            R2[i,j] = 8
        elif (pavetype == 13):
            R2[i,j] = 15
        elif (pavetype == 14):
            R2[i,j] = 5
        elif (pavetype == 16):
            R2[i,j] = 1


        if (vegetype == 1 or vegetype == 2):
            R1[i,j] = -127
            R2[i,j] = 1
        elif (vegetype == 3):
            R1[i,j] = 1
        elif (vegetype == 4):
            R1[i,j] = 16
        elif (vegetype == 5):
            R1[i,j] = 2
        elif (vegetype == 6 or vegetype == 7 or vegetype == 8 or vegetype == 9):
            R1[i,j] = 7
            
        if (buildtype == 1):
            R5[i,j] = 2
        elif (buildtype == 2):
            R5[i,j] = 5
        
        #Pavement/vegetation type must have soil type defined under it. Approximating as medium porosity soil.
        if (soiltype < 1 and (vegetype > 0 or pavetype > 0)):
            R6[i,j] = 2
        if (pavetype > 0 and watertype > 0 and vegetype > 0):
            R3[i,j] = -127

        #Currently, point can only belong to vegetation, pavement or water. Favour pavement.
        if (pavetype > 0 and watertype > 0):
            R3[i,j] = -127
        if (pavetype > 0 and vegetype > 0):
            R1[i,j] = -127
        if (vegetype > 0 and watertype > 0):
            R3[i,j] = -127
        
            
        
#------------------------------------------------

# Info about building type overlayed onto building mask from building height layer.
idx = R5 > 0
R4[idx] = R5[idx]

# Update labels according to PIDS

vegetation_labels = ['Bare soil (no vegetation)', 'Crops, mixed farming', 'Short grass', 'Evergreen needleleaf trees', 
                     'Deciduous needleleaf trees', 'Evergreen broadleaf trees', 'Deciduous broadleaf trees', 'Tall grass', 'Desert', 
                     'Tundra', 'Irrigated crops', 'Semidesert', 'Ice caps and glaciers', 'Bogs and marshes', 'Evergreen shrubs', 
                     'Deciduous shrubs', 'Mixed forest / woodland', 'Interrupted forest']

pavement_labels = ['Unknown pavement (asphalt/concrete mix)', 'Asphalt (asphalt concrete)', 'Concrete (portland concrete)', 
                   'Sett', 'Paving stoned', 'cobblestone', 'Metal', 'Wood', 'Gravel', 'Fine gravel', 'Pebblestone', 'Woodchips', 
                   'Tartan (sports)', 'Artificial turf (sports)', 'Clay (sports)', 'Building (dummy)']

water_labels = ['Lake', 'River', 'Ocean', 'Pond', 'Fountain']


building_labels = ['Residential (-1950)', 'Residential (1951-2000)', 'Residential (2001-)', 'Office (-1950)', 
                   'Office (1951-2000)', 'Office (2001-)']

soil_labels = ['Coarse', 'Medium', 'Medium-fine', 'Fine', 'Very fine', 'Organic']

for i in range (len(vegetation_labels)):
    vegetation_labels[i] = ('{}').format(vegetation_labels[i])
    if (i < len(pavement_labels)):
        pavement_labels[i] = ('{}').format(pavement_labels[i])
        if (i < len(building_labels)):
            building_labels[i] = ('{}').format(building_labels[i])
            if (i < len(soil_labels)):
                soil_labels[i] = ('{}').format(soil_labels[i])
                if (i < len(water_labels)):
                  water_labels[i] = ('{}').format(water_labels[i])
    

print(" ... done.\n")

# Write output data file
print(" Writing output file...")
# Save as Numpy Z file.


Rdict1["R"]=R1
Rdict1["dPx"]=dPx
Rdict1["Labels"] = vegetation_labels
saveTileAsNumpyZ( vegetation, Rdict1 )

Rdict1["R"]=R2
Rdict1["Labels"] = pavement_labels
saveTileAsNumpyZ( pavement, Rdict1 )

Rdict1["R"]=R3
Rdict1["Labels"] = water_labels
saveTileAsNumpyZ( water, Rdict1 )

Rdict1["R"]=R4
Rdict1["Labels"] = building_labels
saveTileAsNumpyZ( buildingmask, Rdict1 )

Rdict1["R"]=R6
Rdict1["Labels"] = soil_labels
saveTileAsNumpyZ( soil, Rdict1 )

Rdict1["R"]=R7
Rdict1["Labels"] = ['']
saveTileAsNumpyZ( building_id, Rdict1 )


print(" ...{}, {}, {}, {}, {}, {} saved successfully.".format(vegetation,pavement,water,buildingmask,soil, building_id))
