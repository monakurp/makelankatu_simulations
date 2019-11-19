import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.measurements as snm
from mapTools import saveTileAsNumpyZ
import netCDF4 as nc
import os

mydir = os.path.expanduser('~/makelankatu_simulations/input_data_to_palm')
os.chdir((mydir))

domains = ['root','parent','child']

di = 0
for d in domains:
  
  datain = np.load('Container/building_id_{}.npz'.format(domains[di] ) ) 
  Rdict = dict( datain )
  datain.close()
  
  Rmod = Rdict['R']
  Rmod[Rmod==-127] = 0

  labeled_array, num_features = snm.label( Rmod )
  labeled_array[labeled_array==0.0] = -127
  
  Rdict['R'] = labeled_array
  saveTileAsNumpyZ( 'building_id_{}.npz'.format( domains[di] ), Rdict )
  
  ncname = 'PIDS_STATIC'
  if di>0:
    ncname = 'PIDS_STATIC_N0{}'.format( di+1 )
  ncdata = nc.Dataset( ncname, 'r+'  )
  ncdata.variables['building_id'][:] = labeled_array[::-1,:]
  ncdata.close()
  
  di += 1
  

  
