#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from mapTools import applyMargins, saveTileAsNumpyZ

file_root = ['topo_root.npz','oro_root.npz','vege_root.npz', 'landuse_root.npz']

# Zero or non-zero margin widths as ratios (0-1): [L,R,B,T]
mw = [0.0, 0.1, 0.0, 0.0]

# Margin ramp widths as ratios (0-1): [L,R,B,T]
mr = [0.0, 0.02, 0.0, 0.0]

# Margins heights: [L,R,B,T].
mh = [0.0, 0.0, 0.0, 0.0]

for i in range( len( file_root ) ):

  root = np.load( file_root[i] )
  Rdict = dict( root )
  root.close()

  R = Rdict['R']

  # Apply margins
  Rf = applyMargins( R , mw, mr, mh )

  # Save file
  Rdict['R'] = Rf
  saveTileAsNumpyZ( file_root[i], Rdict )

#  fig = plt.figure()
#  plt.imshow( Rf, interpolation='None' )
#  plt.show()




