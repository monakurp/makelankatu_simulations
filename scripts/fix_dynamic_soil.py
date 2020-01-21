import netCDF4 as nc
import numpy as np

filename = '/home/monakurp/makelankatu_simulations/input_data_to_palm/cases/20170609_morning/PIDS_DYNAMIC_real_no_vp'

datain = nc.Dataset( filename, 'r+' )

soil_m = datain['init_soil_m'][:]
soil_t = datain['init_soil_t'][:]

for k in range( np.shape( soil_m )[0] ):
  soil_m[k,:,:] = np.max( soil_m[k,:,:,] )
  soil_t[k,:,:] = np.max( soil_t[k,:,:,] )
  
datain['init_soil_m'][:] = soil_m
datain['init_soil_t'][:] = soil_t

datain.close()
