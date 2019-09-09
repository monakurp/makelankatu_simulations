import numpy as np
import os

folder = '/home/monakurp/makelankatu_simulations/source_data'

os.chdir(( folder ))

file_T_mean = 'dailymean_T_Kaisaniemi.csv'
outfile = 'HDD_Kaisaniemi.csv'

# mean temperature
daily_T = np.genfromtxt(file_T_mean, skip_header=1, delimiter=',', usecols=(0,1,2,5) )
HDDi = np.zeros( len( daily_T ) )
for i in range( len( HDDi ) ):
  if daily_T[i,1] <= 6:
    if daily_T[i,-1] <= 10.0:
      HDDi[i] = 17.0 - daily_T[i,-1]
  else:
    if daily_T[i,-1] <= 12.0:
      HDDi[i] = 17.0 - daily_T[i,-1]
      
HDD = np.cumsum(HDDi)
tosave = np.hstack((daily_T,HDD.reshape(len(HDD),1)))
   
      
np.savetxt( outfile, tosave, delimiter=',' )
