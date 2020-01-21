import numpy as np
import matplotlib.pyplot as plt
import scipy.ndimage.measurements as snm
import netCDF4 as nc
import os

'''
Description:
Create a building_id for each individual building
'''

mydir = os.path.expanduser('~/makelankatu_simulations/input_data_to_palm')
os.chdir((mydir))

domains = ['root','parent','child']

di = 0
for d in domains:
  
  # Read in the data
  datain = np.load('building_id_{}.npz'.format(domains[di] ) )
  Rdict = dict( datain )
  datain.close()
  
  # Choose the raster array and temporally set all fill value (i.e. -127) to zero
  Rmod = Rdict['R']
  Rmod[Rmod<0] = 0

  # Label each building using scipy.ndimage.measurements.snm function
  labeled_array, num_features = snm.label( Rmod )

  # Change zeros back to fill values
  labeled_array[labeled_array==0.0] = -9999
  
  # Save the modified array in the .npz file
  Rdict['R'] = labeled_array
  np.savez_compressed( 'building_id_{}.npz'.format( domains[di] ), **Rdict )
  
  # Save the modified array in the PIDS_STATIC file
  ncname = 'PIDS_STATIC'
  if di>0:
    ncname = 'PIDS_STATIC_N0{}'.format( di+1 )
  ncdata = nc.Dataset( ncname, 'r+'  )
  ncdata.variables['building_id'][:] = labeled_array[::-1,:]
  ncdata.close()
  
  di += 1
  

  
