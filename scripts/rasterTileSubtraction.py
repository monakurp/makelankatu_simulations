#!/usr/bin/env python
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mapTools import *

'''
Description: Subtract a raster tile from another raster tile.


Author: Mona Kurppa
        mona.kurppa@helsinki.fi
        University of Helsinki
'''

#==========================================================#
parser = argparse.ArgumentParser( prog = 'rasterTileSubtraction.py')
parser.add_argument( "-f", "--file", type = str,\
  help = "Name of the raster data file that will be modified." )
parser.add_argument( "-fs", "--filesub", type = str, \
  help = "Name of the raster data file that will be subtracted from the original file." )
parser.add_argument( "-fm", "--filemask", type = str, default = None, \
  help = "Name of the raster data file applied as a mask (optional)." )
parser.add_argument( "-mvl", "--minimum_value_limit", type = float,\
  help = "Values smaller than this limit are omitted in the output file.", default = 0.0 )
parser.add_argument( "-ntp", "--negative_to_positive", action="store_true", \
  help = "Set negative values to positive in the end.", default = False ) 
parser.add_argument( "-fo", "--fileout", type = str,\
  help = "Name of the output .npz-file." )
parser.add_argument( "-p", "--printOn",\
  help="Print the resulting raster data.", action = "store_true", default = False)
parser.add_argument("-pp", "--printOnly",\
  help="Only print the resulting data. Don't save.", action="store_true", default = False)
args = parser.parse_args()
#==========================================================#

# Rename some args for brevity.
filein    = args.file
filesub   = args.filesub
filemask  = args.filemask
fileout   = args.fileout
printOn   = args.printOn
printOnly = args.printOnly
min_lim   = args.minimum_value_limit
ntp       = args.negative_to_positive

filenames = [filein, filesub, filemask, fileout]

# Read in the data. The output raster will inherit Rdict1 properties not shown here.
Rdict = readNumpyZTile( filein )
R     = Rdict['R']
Rdims = np.array( np.shape( R ) )
ROrig = Rdict['GlobOrig']
dPx   = Rdict['dPx']

Rdict_sub = readNumpyZTile( filesub )
R_sub     = Rdict_sub['R']
Rdims_sub = np.array( np.shape( R_sub ) )
ROrig_sub = Rdict_sub['GlobOrig']
dPx_sub   = Rdict_sub['dPx']

if ( filemask is None ):
  R_mask     = 0 * np.copy(R) + 1
  Rdims_mask = np.copy( Rdims )
  ROrig_mask = np.copy( ROrig )
  dPx_mask   = np.copy( dPx )
else:
  Rdict_mask = readNumpyZTile( filemask )
  R_mask     = Rdict_mask['R']
  Rdims_mask = np.array( np.shape( R_mask ) )
  ROrig_mask = Rdict_mask['GlobOrig']
  dPx_mask   = Rdict_mask['dPx']

if( ( ROrig == ROrig_sub ).all() and ( ROrig == ROrig_mask ).all() ):
  print(' Excellent! The origos match.')
else:
  print(' The tiles do not have identical origos. Exiting.')
  i = 0
  for origo in [ ROrig, ROrig_sub, ROrig_mask]:
    print('{}: {}'.format( filenames[i], origo ) )
    i += 1
  sys.exit(1)
  
if( ( Rdims == Rdims_sub ).all() and ( Rdims == Rdims_mask ).all() ):
  print(' Excellent! The dimensions match.')
else:
  print(' The tiles do not have identical dimensions. Exiting.')
  i = 0
  for dims in [ Rdims, Rdims_sub, Rdims_mask]:
    print('{}: {}'.format( filenames[i], dims ) )
    i += 1
  sys.exit(1)
  
# Substract the raster to be subtracted from the original raster
R_mod = R - R_sub

# Create a minimum value mask and filter
min_lim_mask = R < min_lim 
R_mod[min_lim_mask] = 0

# Use the mask file to further filter the original raster
R_mod[R_mask==0] = 0

# Further set mask value to zero where the original raster map is set to zero
R_mask[R_mod==0] = 0

# Negative values to positive
if ntp:
  nvp_mask = ( R_mod<0 ) & ( R_mask>0 )
  R_mod[nvp_mask] *= -1.0
  R_sub[nvp_mask] += R_mod[nvp_mask]

Rdict['R'] = R_mod
Rdict_sub['R'] = R_sub
Rdict_mask['R'] = R_mask


if ( not printOnly ):
  if ( fileout is None ):
    print('Name of the output .npz-file missing! Exiting.')
    sys.exit(1)
  saveTileAsNumpyZ( fileout, Rdict )
  if ( filemask is not None ):
    saveTileAsNumpyZ( '{}_mod.npz'.format( filemask[0:-4] ), Rdict_mask  )
  if ntp:
    saveTileAsNumpyZ( '{}_mod.npz'.format( filesub[0:-4] ), Rdict_sub )

if ( printOn or printOnly ):
  i = 0
  for Ri in [R, R_sub, R_mask, R_mod]:
    fig = plt.figure()
    ax = fig.add_subplot(111)
    if i<3:
      im = ax.imshow(Ri)
    else:
      im = ax.imshow(Ri, cmap='RdBu_r', vmin=-5, vmax=5)
      cbar = plt.colorbar( im )	
      ax.grid(True)
    i += 1
  plt.show()
