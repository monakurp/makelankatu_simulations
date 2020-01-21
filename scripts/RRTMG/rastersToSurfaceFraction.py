#!/usr/bin/env python
import sys
import argparse
import numpy as np
from utilities import writeLog
from mapTools import *
from netcdfTools import *


'''
Description:
Generates surface fractions for each surface layer (vegetation, pavement, water) for when there's overlap between layers.

Author: Jani Stromberg
        jani.stromberg@helsinki.fi
        University of Helsinki
'''

#==========================================================#
parser = argparse.ArgumentParser(prog='rastersToSurfacefraction.py', description='''Writes surface_fractions for PALM from input raster data.''')
parser.add_argument("-f0","--vegetation", type=str, \
  help="Name of vegetation type input raster data file.")
parser.add_argument("-f1","--pavement", type=str, \
  help="Name of pavement type input raster data file.")
parser.add_argument("-f2","--water", type=str, \
  help="Name of water type input raster data file.")
parser.add_argument("-fo", "--fileout", type=str, \
  help="Name of the output 3D data file. ")
args = parser.parse_args()
writeLog(parser, args)

#==========================================================#
# Renaming ... that's all
vegetation = args.vegetation
pavement = args.pavement
water = args.water
dz       = 1
fileout  = args.fileout

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

if( (R1Orig == R2Orig).all() and (R2Orig == R3Orig).all() ):
  print(' Excellent! The origos match.')
else:
  print(' The tiles do not have identical origos. Exiting.')
  sys.exit(1)

if( (R1dims == R2dims).all() and (R2dims == R3dims).all() ):
  print(' Excellent! The dimensions match.')
else:
  print(' The tiles do not have identical dimensions. Exiting.')
  sys.exit(1)


# Calculate the shape of the new 3D array
nPx3D = nPx; dPx3D = dPx
dPx3D = np.append(dPx, dz)
nPx3D = np.append(nPx3D, 3)


# Fill the 3D surface fraction array
nsurface_frac = np.zeros([nPx3D[1], nPx3D[0], nPx3D[2]])
nPc=[nPx3D[1], nPx3D[0], nPx3D[2]]
dPc=[dPx3D[1], dPx3D[0], dPx3D[2]]
print(" 3D grid array dimensions [x,y,z]: {}, {}, {}".format(*nPc))
print(" Generating fractions ...")

#ToDo Fix: for some reason, where there's overlap between all three layers, PALM calculates the summed fractions
#as >1. 
for i in range(nPc[1]):
  for j in range(nPc[0]):
    vege  = R1[i,j] #nsurface_frac[:,:,0]
    pavem = R2[i,j] #nsurface_frac[:,:,1]
    water = R3[i,j] #nsurface_frac[:,:,2]
    


    if (vege > 0 and pavem > 0):
        nsurface_frac[i,j,0] = 0.5
        nsurface_frac[i,j,1] = 0.5
    elif (vege > 0 and water > 0):
        nsurface_frac[i,j,0] = 0.5
        nsurface_frac[i,j,2] = 0.5
    elif (pavem > 0 and water > 0):
        nsurface_frac[i,j,1] = 0.5
        nsurface_frac[i,j,2] = 0.5
    elif (vege > 0):
        nsurface_frac[i,j,0] = 1
    elif (pavem > 0):
        nsurface_frac[i,j,1] = 1
    elif (water > 0):
        nsurface_frac[i,j,2] = 1
    else:
        nsurface_frac[i,j,:] = -9999.0



print(" ... done.\n")

# Write output data file
print(" Writing output file...")
  
# Save as Numpy Z file.

Rdict1["R"]=nsurface_frac
Rdict1["dPx"]=dPx3D
saveTileAsNumpyZ( fileout, Rdict1 )
print(" ...{} saved successfully.".format(fileout))
