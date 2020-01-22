#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 17 11:00:10 2019

@author: stromjan
"""
import numpy as np
from mapTools import *

filename = '/home/stromjan/Maps/npz_files/Topography_H3D.npz'
filename2 = '/home/stromjan/Maps/npz_files/Land_cover_types_H3D.npz'

Rdict = readNumpyZTile(filename)
Rdict2 = readNumpyZTile(filename2)

R = Rdict['R']
GlobOrig = Rdict['GlobOrig']
Rdims = np.array(np.shape(R))

R2 = Rdict2['R']
GlobOrig2 = Rdict2['GlobOrig']
Rdims2 = np.array(np.shape(R2))


# Get indices of container location
y_container = GlobOrig[0] - 6675954
x_container = 25497343 - GlobOrig[1]

#container height 3 meters
R[y_container-2,x_container] += 3
R[y_container-1,x_container] += 3
R[y_container-1,x_container-1] += 3
R[y_container-3,x_container-1] += 3
R[y_container-2,x_container-1] += 3
R[y_container-2,x_container-2] += 3
R[y_container-2,x_container-3] += 3
R[y_container-3,x_container-2] += 3
R[y_container-3,x_container-3] += 3
R[y_container-3,x_container-4] += 3
R[y_container-4,x_container-3] += 3
R[y_container-4,x_container-4] += 3
R[y_container-4,x_container-5] += 3
R[y_container-5,x_container-4] += 3
R[y_container-5,x_container-5] += 3


#Land cover type set as 6 (building)
R2[y_container-2,x_container] = 6
R2[y_container-1,x_container] = 6
R2[y_container-1,x_container-1] = 6
R2[y_container-3,x_container-1] = 6
R2[y_container-2,x_container-1] = 6
R2[y_container-2,x_container-2] = 6
R2[y_container-2,x_container-3] = 6
R2[y_container-3,x_container-2] = 6
R2[y_container-3,x_container-3] = 6
R2[y_container-3,x_container-4] = 6
R2[y_container-4,x_container-3] = 6
R2[y_container-4,x_container-4] = 6
R2[y_container-4,x_container-5] = 6
R2[y_container-5,x_container-4] = 6
R2[y_container-5,x_container-5] = 6


Rdict['R'] = R
Rdict2['R'] = R2

fout = '/home/stromjan/Maps/npz_files/Topography_H3D_Container.npz'
fout2 = '/home/stromjan/Maps/npz_files/Land_cover_types_H3D_Container.npz'

saveTileAsNumpyZ( fout, Rdict )
saveTileAsNumpyZ( fout2, Rdict2 )

