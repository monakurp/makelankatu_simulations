import numpy as np
import netCDF4 as nc
import os
import sys
import datetime as dtime
import matplotlib.pyplot as plt

os.chdir(os.path.expanduser("~"))
os.chdir(('makelankatu_simulations'))
np.set_printoptions( precision=10 )

parameter = True
variable  = False

# -------------------------------------------------------------------------------------------------#
# PROVIDE THESE:
# -------------------------------------------------------------------------------------------------#

# Simulation time: start of drone measurements - ~15min --> end of drone measurements
sim_year    = 2017
sim_month   = 6 # 12
sim_day     = 9 # 7
sim_time    = 'evening'
datestr     = '{}{:02d}{:02d}'.format( sim_year, sim_month, sim_day )
adchem_type = 'FRES' # or 'orig'

flow_spinup_min = 15

if ( datestr=='20170609' and sim_time=='morning' ):
  start_hour = 7
  start_min  = 15 - flow_spinup_min
  end_hour   = 9
  end_min    = 15
  plusUTC    = 3
  
if ( datestr=='20170609' and sim_time=='evening' ):
  start_hour = 20
  start_min  = 10
  end_hour   = 21
  end_min    = 15
  plusUTC    = 3

elif ( datestr=='20170614' and sim_time=='morning' ):
  start_hour = 6
  start_min  = 55
  end_hour   = 9
  end_min    = 0
  plusUTC    = 3
  
elif ( datestr=='20171207' and sim_time=='morning' ):
  start_hour = 7
  start_min  = 5
  end_hour   = 9
  end_min    = 15
  plusUTC    = 2


orig     = [60.1663312, 24.873065935] # 900 m shift to the left, was 24.8895075

# PALM grid

grid_type = 'real' # real or test

if grid_type == 'real':
  nx = 768-1
  ny = 768-1
  nz = 80
  
  dx =  9.0 
  dy =  9.0
  dz =  6.0
  
  dz_stretch_factor = 1.03 # was 1.04
  dz_stretch_level = 300.0 # was 400.0
  dz_max = 18.0
  
  nx_N03 = 576-1
  ny_N03 = 576-1
  nz_N03 = 144 # was 96

  dx_N03 = 1.0
  dy_N03 = 1.0
  dz_N03 = 1.0   
  
  dt = 3600.0

elif grid_type == 'test':
  nx = 15
  ny = 15
  nz = 24   
  
  dx =  10.0 
  dy =  10.0
  dz =  10.0
  
  dz_stretch_factor = 1.00
  dz_stretch_level = 60.0
  dz_max = 6.0
  
  nx_N03 = 19
  ny_N03 = 19
  nz_N03 = 60

  dx_N03 = 2.0
  dy_N03 = 2.0
  dz_N03 = 2.0 
  
  dt = 360.0

ntimesteps = end_hour - start_hour + 1
if end_min > 0:
  ntimesteps += 1

# -------------------------------------------------------------------------------------------------#

#%% Filenames

fname_full_backup = 'source_data/cases/{}_{}/meps_mbr0_full_backup_2_5km_{}T00Z.nc'.format( \
                                       datestr, sim_time, datestr )
fname_subset      = 'source_data/cases/{}_{}/meps_subset_2_5km_{}T00Z.nc'.format( \
                                       datestr, sim_time, datestr )
fname_kivenlahti = '../Katumetro/kivenlahti_mast_{}{}{}.txt'.format( sim_day, sim_month, sim_year )
fname_kumpula = '../Katumetro/smeardata_{}{:02d}{:02d}120000.txt'.format( sim_year, sim_month, sim_day )


#%% Corners of the simulation domain

radius_earth = 6371229.0 # m
npoints_interpolation = 4

lats = np.zeros( [npoints_interpolation], dtype=float ); lats[0] = orig[0]
lons = np.zeros( [npoints_interpolation], dtype=float ); lons[0] = orig[1]

r_lat      = radius_earth * np.cos( orig[0]*np.pi/180.0 )
circle     = 2.0 * np.pi * radius_earth
circle_lat = 2.0 * np.pi * r_lat

for i in range( len( lats )-1 ):
  lats[i+1] = lats[0] + (i+1)/ (npoints_interpolation-1) * (ny+1)*dy / circle * 360.0
  lons[i+1] = lons[0] + (i+1)/ (npoints_interpolation-1) * (nx+1)*dx / circle_lat * 360.0
  
print('latitudes:')
print("".join(['{}, '.format(i)*npoints_interpolation for i in lats])) 
print('')
print('longitudes:')
print("".join(['{}, '.format(i) for i in lons]) *npoints_interpolation)
  
#%% Read in variables:

full_backup = nc.Dataset( fname_full_backup )
subset = nc.Dataset( fname_subset )

# Check also that the files are for the same time and coordinates!

time = np.array( full_backup['time'] ) # "seconds since 1970-01-01 00:00:00 +00:00"
timestring = dict()
if ( np.array( subset['time'][:] ) != time[:] ).any():
  sys.exit('Files do not have the same time stamps!')
for t in range( len( time ) ):
  timestring[t] = dtime.datetime.utcfromtimestamp( time[t] )
  if ( timestring[t].hour + plusUTC ) == start_hour:  
    ss = t
  if ( timestring[t].hour + plusUTC ) == end_hour: 
    ee = t
if end_min > 0:
  ee += 1
    
longitude = np.array( full_backup['longitude'] )
if ( np.array( subset['longitude'][:] ) != longitude[:] ).any():
  sys.exit('Files do not have the same longitudes!')
longitude = np.reshape( longitude, ( len(lats),len(lons) ) )
  
latitude = np.array( full_backup['latitude'] )
if ( np.array( subset['latitude'][:] ) != latitude[:] ).any():
  sys.exit('Files do not have the same latitudes!')
latitude = np.reshape( latitude, ( len(lats),len(lons) ) )

ensemble = np.array( full_backup['ensemble_member'] )
if ( np.array( subset['ensemble_member'][:] ) != ensemble[:] ).any():
  sys.exit('Files do not have the same ensemble member!')

hybrid = np.array( full_backup['hybrid'] ); hybrid = np.flipud( hybrid )
Z_0 = np.array( subset['surface_geopotential'] )
Z_0 = np.squeeze( np.squeeze( Z_0, axis=2 ), axis=1 )
Z_0 = Z_0[ss:ee+1,:,:]
Z_0 = np.reshape( Z_0, ( ntimesteps, len(lats), len(lons) ) )

T = np.array( full_backup['air_temperature_ml'] )  # in K
T = np.squeeze( T, axis=2 )
T = T[ss:ee+1,::-1,:,:]
T = np.reshape( T, ( ntimesteps, len(hybrid), len(lats),len(lons) ) )

q = np.array( full_backup['specific_humidity_ml'] )  # in kg/kg
q = np.squeeze( q, axis=2 )
q = q[ss:ee+1,::-1,:,:]
q = np.reshape( q, ( ntimesteps, len(hybrid), len(lats),len(lons) ) )

u = np.array( full_backup['x_wind_ml'] )  # in m/s
u = np.squeeze( u, axis=2 ) 
u = u[ss:ee+1,::-1,:,:]
u = np.reshape( u,( ntimesteps, len(hybrid), len(lats),len(lons) ) )

v = np.array( full_backup['y_wind_ml'] )  # in m/s
v = np.squeeze( v, axis=2 )
v = v[ss:ee+1,::-1,:,:]
v = np.reshape( v,( ntimesteps, len(hybrid), len(lats),len(lons) ) )

w = np.array( full_backup['upward_air_velocity_ml'] )  # in m/s
w = np.squeeze( w, axis=2 )
w = w[ss:ee+1,::-1,:,:]
w = np.reshape( w,( ntimesteps, len(hybrid), len(lats),len(lons) ) )

ap = np.array( full_backup['ap'] )
b  = np.array( full_backup['b'] )
ps = np.array( subset['surface_air_pressure'] )
ps = ps[ss:ee+1,:,:]
ps = np.reshape( ps,( ntimesteps, len(lats),len(lons) ) )

full_backup.close()
subset.close()

#%% Define the pressure and height levels

# p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)"

p = np.zeros( [ ntimesteps, len(hybrid), len(lats), len(lons) ] )

for t in range( ntimesteps ):
  for h in range( len( hybrid ) ):
    p[t,h,:,:] = ap[h] + b[h] * ps[t,:,:]
p = p[:,::-1,:,:]
    
# Z_2 = Z_1 - T_avg * R/g0 * ln(p_2/p_1)
R = 287.0 # J/kg/K
g0 = 9.80665 # m/s2
T_avg = 0.5 * ( T[:,0:-1,:,:] + T[:,1::,:,:] )

# geopotential height
Z = np.zeros( [ ntimesteps, len(hybrid), len(lats), len(lons) ] ) + np.nan

for t in range( ntimesteps ):
  for h in range( len( hybrid ) - 1 ):
    if h==0:
      Z[t,h,:,:] = Z_0[t,:,:]/g0
    else:
      Z[t,h,:,:] = Z[t,h-1,:,:] + T_avg[t,h-1,:,:] * R/g0 * np.log( p[t,h-1,:,:]/p[t,h,:,:] )
Z[:,-1,:,:]= Z[:,-2,:,:]
      
#%% Potential temperature

pt = np.zeros( [ ntimesteps, len(hybrid), len(lats), len(lons) ] ) + np.nan
for t in range( ntimesteps ):
  for h in range( len( hybrid ) ):
    pt[t,h,:,:] = T[t,h,:,:] * ( 100000.0 / p[t,h,:,:] ) ** 0.286
    
#%% Kivenlahti data

# 0 = year	
# 1 = month
# 2 = day
# 3 = hour(utc+3)
# 4 = minute
# 5 = T_2
# 6 = U_2
# 7 = WD_2
# 8 = T_26
# 9 = U_26
# 10 = WD_26
# 11 = T_49
# 12 = U_49
# 13 = WD_49
# 14 = T_92
# 15 = U_92
# 16 = WD_92
# 17 = T_141
# 18 = U_141
# 19 = WD_141
# 20 = T_217
# 21 = U_217
# 22 = WD_217
# 23 = T_265
# 24 = U_265
# 25 = WD_265
# 26 = T_327
# 27 = U_327
# 28 = WD_327

kivenlahti = np.loadtxt( fname_kivenlahti, skiprows=1 )
kivenlahti[kivenlahti==-999] = np.nan

WDk = np.zeros([ntimesteps,8])
Uk = np.zeros([ntimesteps,8])
zk = np.array( [2, 26, 49, 92, 141, 217, 265, 327] ) 
cols = [ 6, 9, 12, 15, 18, 21, 24, 27 ] 

h = start_hour
t = 0
while t < ntimesteps:
  ind = ( kivenlahti[:,3]==h )
  for zi in range( 8 ):
    Uk[t,zi] = np.nanmean( kivenlahti[ind,cols[zi]] )
    WDk[t,zi] = np.nanmean( kivenlahti[ind,cols[zi]+1] )

  t += 1
  h += 1
  
#%% Kumpula data

kumpula = np.loadtxt( fname_kumpula )
Uku = np.zeros( ntimesteps ) 
WDku = np.zeros( ntimesteps ) 

h = start_hour
t = 0
while t < ntimesteps:
  ind = ( kumpula[:,3] == h )
  Uku[t] = np.nanmean( kumpula[h,6] )
  WDku[t] = np.nanmean( kumpula[h,7] )
  
  h += 1
  t += 1
  
#%% Plot
plt.close('all')
  
#WDh = np.arctan2( -u, -v ) * 180.0 / np.pi
WDh = 270.0 - ( np.arctan2( v, u ) * 180.0 / np.pi )
WDh[WDh>360] = WDh[WDh>360]-360
Uh = np.sqrt( u**2 + v**2 )

fig1 = plt.figure()
fig1.text( 0.5, 0.02, 'U (m/s)', ha='center', fontsize=12 )
fig2 = plt.figure()
fig2.text( 0.5, 0.02, 'WD (deg)', ha='center', fontsize=12 )
for figi in [fig1, fig2]:
  figi.subplots_adjust( hspace=0.4, wspace=0.1, left=0.1, right=0.95, bottom=0.11)
  
for t in range( ntimesteps ):
  ax1 = fig1.add_subplot(1, ntimesteps, t+1)
  ax2 = fig2.add_subplot(1, ntimesteps, t+1)
  for j in range( len( lats ) ):
    for i in range( len( lons ) ):
      if j==0 and i==0 and t==0:
        label = 'Harmonie'
      else:
        label='_nolegend_'
      ax1.plot( Uh[t,:,j,i], Z[t,:,j,i], label=label )
      ax2.plot( WDh[t,:,j,i], Z[t,:,j,i], label=label )
      
  if t==0:
    label='Kivenlahti'
  else:
    label='_nolegend_'    
  ax1.plot( Uk[t,:], zk, 'r*', ms=10, label=label )   
  ax2.plot( WDk[t,:], zk, 'r*', ms=10, label=label  )  
  
  if t==0:
    label='Kumpula'
  else:
    label='_nolegend_'    
  ax1.plot( Uku[t], [32], 'b*', ms=10, label=label)
  ax2.plot( WDku[t], [32], 'b*', ms=10, label=label)
  
  ax1.set_xlim( 0, 10 )
  
  for axi in [ax1, ax2]:
    axi.set_ylim( 0, 400 )
    axi.grid(True)
    if t==0:
      axi.set_ylabel( 'z (m)' )
    else:
      axi.set_yticklabels([])
  ax1.set_xticks([0, 2, 4, 6, 8])
  ax2.set_xticks([90,180,270,360])
  
for figi in [fig1, fig2]:
  figi.legend( loc=9, ncol=3)
  
fig1.savefig( '/home/monakurp/Katumetro/Figures/compare_harmonie_U_{}{:02d}{:02d}_{}.png'.format( 
  sim_year, sim_month, sim_day, sim_time ), format='png', dpi=300 )
fig2.savefig( '/home/monakurp/Katumetro/Figures/compare_harmonie_WD_{}{:02d}{:02d}_{}.png'.format( 
  sim_year, sim_month, sim_day, sim_time ), format='png', dpi=300 )
   