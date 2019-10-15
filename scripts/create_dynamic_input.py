import numpy as np
import netCDF4 as nc
import os
import sys
import datetime as dtime
from scipy.interpolate import interp1d, interp2d
import matplotlib.pyplot as plt
import hdf5storage
from netcdfTools import createNetcdfVariable

os.chdir(('/home/stromjan/Maps/Scripts/Dynamic'))

parameter = True
variable  = False

# -------------------------------------------------------------------------------------------------#
# PROVIDE THESE:
# -------------------------------------------------------------------------------------------------#

# Simulation time:
sim_year   = 2017
sim_month  = 6 # 12
sim_day    = 9 # 7
start_hour = 7 #7
end_hour   = 9 #9
plusUTC    = 3 # 2 h
sim_time   = 'morning' #'morning'

if sim_month==6 and sim_day==9:
  orig     = [60.1663312, 24.873065935] # 900 m shift to the left
else:
  orig     = [60.1663312,24.8895075] 

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
  nz_N03 = 96

  dx_N03 = 1.0
  dy_N03 = 1.0
  dz_N03 = 1.0    

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

ntimesteps = end_hour - start_hour + 1
dt = 3600.0
# -------------------------------------------------------------------------------------------------#

#%% Filenames

datestr = '{}{:02d}{:02d}'.format( sim_year, sim_month, sim_day )

fname_full_backup = '/home/stromjan/Maps/Scripts/Dynamic/meps_mbr0_full_backup_2_5km_{}T00Z.nc'.format( datestr, sim_time, datestr )
fname_subset      = '/home/stromjan/Maps/Scripts/Dynamic/meps_subset_2_5km_{}T00Z.nc'.format( datestr, sim_time, datestr )
fname_soil        = '/home/stromjan/Maps/Scripts/Dynamic/reanalysis_soil.nc'
fname_adchem      = '/home/stromjan/Maps/Scripts/Dynamic/09062017.mat'

fname_out     = '/home/stromjan/Maps/Scripts/Dynamic/PIDS_DYNAMIC_{}'.format( datestr, sim_time, grid_type )
fname_out_N03 = '/home/stromjan/Maps/Scripts/Dynamic/PIDS_DYNAMIC_{}_N03'.format( datestr, sim_time, grid_type )

dsout     = nc.Dataset( fname_out, 'w' )
dsout_N03 = nc.Dataset( fname_out_N03, 'w' )

# Set global attributes:
for dsouti in [dsout, dsout_N03]:
  dsouti.title          = "PALM dynamic input file for scenario Makelankatu, morning 9 Jun 2017"
  dsouti.institution    = 'INAR, University of Helsinki'
  dsouti.author         = 'Mona Kurppa (mona.kurppa@helsinki.fi)'
  dsouti.references     = "--"
  dsouti.comment        = 'Data retrieved from MEPS data (model Harmonie, data archive: http://thredds.met.no/thredds/catalog/meps25epsarchive/catalog.html)'
  dsouti.origin_lat     = orig[0]
  dsouti.origin_lon     = orig[1]
  dsouti.origin_z       = 0.0
  dsouti.origin_time = "{}-{:02d}-{:02d} {:02d}:00:00 +{:02d}".format(sim_year, sim_month, sim_day, start_hour, plusUTC) #Was previously missing
  dsouti.palm_version   = "6.0"
  dsouti.rotation_angle = 0.0

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
      

#%% Soil temperature and moisture
    
soil = nc.Dataset( fname_soil )

soil_lon = np.array( soil['longitude'] )
soil_lat = np.array( soil['latitude'] )
soil_time = np.array( soil['time'] )

soil_timestr = [dtime.datetime( 1900, 1, 1 ) + dtime.timedelta( hours=float(t1) ) for t1 in soil_time]

zsoil = np.array([0.035, 0.175, 0.64, 1.945]) # https://confluence.ecmwf.int/pages/viewpage.action?pageId=56660259
zsoil_PALM = np.array([0.01, 0.02, 0.04, 0.06, 0.14, 0.26, 0.54, 1.86]) #PALM default soil configuration


if grid_type=='real':
  loni = np.where( ( soil_lon > lons[0] ) & ( soil_lon < lons[-1]+0.25 ) )[0][0]
  lati = np.where( ( soil_lat > lats[0] ) & ( soil_lat < lats[-1]+0.25 ) )[0][0]
else:
  loni = 100
  lati = 119

soil_t = np.zeros( len(zsoil), dtype=float )
soil_m = np.zeros( len(zsoil), dtype=float )

if sim_month < 7:
  ti = 0
else:
  ti = 1

for i in range( len( zsoil ) ):
  soil_t[i]  = soil['stl{}'.format(i+1)][ti,lati,loni]
  soil_m[i] = soil['swvl{}'.format(i+1)][ti,lati,loni]


#soil_t = [282.16, 282.09, 281.95, 281.82, 281.28, 280.51, 278.94, 275.37]
#soil_m = [0.322, 0.322, 0.324, 0.325, 0.331, 0.340, 0.357, 0.403] #Last one was 0.404, maximum for coarse soil is 0.403

#fitting a curve to get soil variables at PALM soil levels
soil_t_fit = np.poly1d(np.polyfit(zsoil,soil_t,2))
soil_m_fit = np.poly1d(np.polyfit(zsoil,soil_m,2))
zsoil = zsoil_PALM
soil_t = soil_t_fit(zsoil)
soil_m = soil_m_fit(zsoil)
soil_m[soil_m > 0.403] = 0.403 #check that soil_m doesn't exceed minimum saturation moisture of soil_types

#%% PALM grid:

# PALM model extents:
lx = (nx+1) * dx
ly = (ny+1) * dy

# z (nz) and zw (nz-1):
z_column = np.zeros( nz+1, dtype=float )
z_column[0] = dz * 0.5
dz_stretched = dz
dz_stretch_level_start = dz_stretch_level
for k in range( 1, len(z_column), 1):
  z_column[k] = z_column[k-1] + dz_stretched
   # do stretching if needed
  if ( ( dz_stretch_level_start <= z_column[k] ) and ( dz_stretch_factor>1.0 ) ):
    dz_stretched = np.minimum( dz_stretched * dz_stretch_factor, dz_max)
    z_column[k] = z_column[k-1] + dz_stretched
z = z_column[0:-1]
zw_column = 0.5 * ( z_column[0:-1] + z_column[1::] )
zw = zw_column[0:-1]
  
# x (nx+1) and xu (nx):
x = np.linspace( 0.5*dx, lx - 0.5*dx, nx+1)
xu = np.linspace( dx, lx - dx, nx)

# y (ny+1) and yv (ny):
y = np.linspace( 0.5*dy, ly - 0.5*dy, ny+1)
yv = np.linspace( dy, ly - dy, ny)

# time:
time_palm = np.arange(0.0, (ntimesteps-1)*dt+1, dt )
if grid_type=='test':
  time_palm = 0.1*time_palm


#%% Save dimensions to the dynamic input file:

# x:
dsout.createDimension( 'x', len(x) )
xv = dsout.createVariable( 'x', 'f4', ('x',) )
xv.units = "m"
xv[:] = x
xv.standard_name = "x coordinate of cell centers"

# xu:
dsout.createDimension( 'xu', len(xu) )
xuv = dsout.createVariable( 'xu', 'f4', ('xu',) )
xuv.units = "m"
xuv[:] = xu
xuv.standard_name = "x coordinate of cell faces"

# y:
dsout.createDimension( 'y', len(y) )
yvar = dsout.createVariable( 'y', 'f4', ('y',) )
yvar.units = "m"
yvar[:] = y
yvar.standard_name = "y coordinate of cell centers"

# yv:
dsout.createDimension( 'yv', len(yv) )
yvv = dsout.createVariable( 'yv', 'f4', ('yv',) )
yvv.units = "m"
yvv[:] = yv
yvv.standard_name = "y coordinate of cell faces"

# z:
dsout.createDimension( 'z', len(z) )
zv = dsout.createVariable( 'z', 'f4', ('z',) )
zv.units = "m"
zv[:] = z
zv.standard_name = "z coordinate of cell centers"

# zw:
dsout.createDimension( 'zw', len(zw) )
zwv = dsout.createVariable( 'zw', 'f4', ('zw',) )
zwv.units = "m"
zwv[:] = zw
zwv.standard_name = "z coordinate of cell faces"

# zsoil:
dsout.createDimension( 'zsoil', len(zsoil) )
zsoilv = dsout.createVariable( 'zsoil', 'f4', ('zsoil',) )
zsoilv.units = "m"
zsoilv.positive = "down"
zsoilv[:] = zsoil
zsoilv.standard_name = "depth_below_land"


# time:
dsout.createDimension( 'time', len(time_palm) )
timev = dsout.createVariable( 'time', 'f4', ('time',) )
timev.units = 'seconds since {:04d}{:02d}{:02d} {}:00 UTC'.format( sim_year, sim_month, sim_day, start_hour)
timev[:] = time_palm
timev.standard_name = "time"
timev.long_name = "time"

#%% Create horizontal homogeneous fields of soil temperature and moisture

init_soil_t = np.zeros( [len(zsoil), ny+1, nx+1], dtype=float )
init_soil_m = np.zeros( [len(zsoil), ny+1, nx+1], dtype=float )

for i in range( len( zsoil ) ):
  init_soil_t[i,:,:] = soil_t[i]
  init_soil_m[i,:,:] = soil_m[i]

#%% Interpolate variable fields to PALM grid

u_palm  = np.zeros( [ ntimesteps, nz,   ny+1, nx  ], dtype=float )
v_palm  = np.zeros( [ ntimesteps, nz,   ny,   nx+1], dtype=float )
w_palm  = np.zeros( [ ntimesteps, nz-1, ny+1, nx+1], dtype=float )
pt_palm = np.zeros( [ ntimesteps, nz,   ny+1, nx+1], dtype=float )
qv_palm = np.zeros( [ ntimesteps, nz,   ny+1, nx+1], dtype=float )

arrays_out = [u_palm, v_palm, w_palm, pt_palm, qv_palm]
arrays_in  = [u,      v,      w,      pt,      q      ]
zs         = [z,      z,      zw,     z,       z      ]
nys        = [ny+1,   ny,     ny+1,   ny+1,    ny+1   ]
nxs        = [nx,     nx+1,   nx+1,   nx+1,    nx+1   ]
names  = ['u',   'v',    'w',    'pt',    'qv'    ]

print('Interpolate:')
for a in range(  len( arrays_out ) ):

  print('   {}'.format(names[a]) )  
  
  # x- and y-axis based on longitude and latitude
  xlon = np.linspace( np.min(longitude), np.max(longitude), nxs[a] )
  ylat = np.linspace( np.min(latitude),  np.max(latitude),  nys[a] )
  
  for t in range( ntimesteps ):
    temporary = np.zeros( [ len(hybrid), nys[a], nxs[a] ], dtype=float )
    Z_temporary = np.zeros( [ len(hybrid), nys[a], nxs[a] ], dtype=float )
    
    # first interpolate in horizontal:
    for k in range( len(hybrid) ):
      fhor = interp2d( longitude, latitude, arrays_in[a][t,k,:,:], kind='cubic', fill_value="extrapolate" )
      temporary[k,:,:] = fhor( xlon, ylat )
      
      fhorZ = interp2d( longitude, latitude, Z[t,k,:,:], kind='cubic', fill_value="extrapolate" )
      Z_temporary[k,:,:] = fhorZ( xlon, ylat )
    
    for j in range( nys[a] ):
      for i in range( nxs[a] ):
        fver = interp1d( Z_temporary[:,j,i], temporary[:,j,i], kind='linear', fill_value="extrapolate" )
        arrays_out[a][t,:,j,i] = fver( zs[a] )
del arrays_in, temporary, Z_temporary
    
u_palm  = arrays_out[0]
v_palm  = arrays_out[1]
w_palm  = arrays_out[2]
pt_palm = arrays_out[3]
qv_palm = arrays_out[4]   

del arrays_out

#%% Dynamic input variables

arrays = [u_palm, v_palm, w_palm, pt_palm, qv_palm]
names  = ['u',   'v',    'w',    'pt',    'qv'    ]
zs     = ['z',   'z',    'zw',   'z',     'z'     ]
ys     = ['y',   'yv',   'y',    'y',     'y'     ]
xs     = ['xu',  'x',    'x',    'x',     'x'     ]
units  = ["m/s", "m/s",  "m/s",  "K",     "kg/kg" ]

# -------------------------------------------------------------------------------------------------#
# initialisation:
# -------------------------------------------------------------------------------------------------#

# init_atmosphere_XX (z/zw)
lnames = ["initial wind component in x direction", "initial wind component in y direction",\
          "initial wind component in z direction", "initial potential temperature",\
          "initial specific humidity"]
          
for a in range( len( arrays ) ):
  var   = np.nanmean( np.nanmean( arrays[a][0,:,:,:], axis=1), axis=1 )
  ncvar = dsout.createVariable( 'init_atmosphere_{}'.format( names[a] ), 'f4', ( zs[a], ), 
                                fill_value=-9999.0 )
  ncvar[:] = var
  ncvar.units = units[a]
  ncvar.source = "MEPS analysis for 20170609"
  ncvar.long_name = lnames[a]
  ncvar.standard_name = ""
  ncvar.lod = 1

# -------------------------------------------------------------------------------------------------#
# forcing:
# -------------------------------------------------------------------------------------------------#

# surface_forcing_surface_pressure (time)
var   = np.nanmean( np.nanmean( p[:,0,:,:], axis=1 ), axis=1 )
ncvar = dsout.createVariable( 'surface_forcing_surface_pressure', 'f4', ( 'time', ), 
                              fill_value=-9999.0 )
ncvar[:] = var
ncvar.units = "Pa"
ncvar.source = "MEPS analysis for 20170609"
ncvar.long_name = "surface pressure"
ncvar.standard_name = ""
ncvar.lod = 1


# ls_forcing_ug/vg (time,z) 
lnames = ["geostrophic wind (u component)", "geostrophic wind (v component)"]

for a in range(2):
  var   = np.zeros( [ntimesteps, nz], dtype=float )
  ncvar = dsout.createVariable( 'ls_forcing_{}g'.format( names[a] ), 'f4', ( 'time', zs[a] ), 
                                fill_value=-9999.0 )
  ncvar[:] = var
  ncvar.units = units[a]
  ncvar.source = "MEPS analysis for 20170609"
  ncvar.long_name = lnames[a]
  ncvar.standard_name = ""
  ncvar.lod = 2


# ls_forcing_left_XX (time,z/zw,y/yv)
lnames = ["large-scale forcing for left model boundary for the wind component in x direction",\
          "large-scale forcing for left model boundary for the wind component in y direction",\
          "large-scale forcing for left model boundary for the wind component in z direction",\
          "large-scale forcing for left model boundary for the potential temperature",\
          "large-scale forcing for left model boundary for the specific humidity"]

for a in range( len( arrays ) ):
  var   = arrays[a][:,:,:,0]
  ncvar = dsout.createVariable( 'ls_forcing_left_{}'.format( names[a] ), 'f4', ( 'time', zs[a], ys[a] ), 
                                fill_value=-9999.0 )
  ncvar[:] = var
  ncvar.units = units[a]
  ncvar.source = "MEPS analysis for 20170609"
  ncvar.long_name = lnames[a]
  ncvar.standard_name = ""
  ncvar.lod = 2

  
# ls_forcing_right_XX (time,z/zw,y/yv)
lnames = ["large-scale forcing for right model boundary for the wind component in x direction",\
          "large-scale forcing for right model boundary for the wind component in y direction",\
          "large-scale forcing for right model boundary for the wind component in z direction",\
          "large-scale forcing for right model boundary for the potential temperature",\
          "large-scale forcing for right model boundary for the specific humidity"]

for a in range( len( arrays ) ):
  var   = arrays[a][:,:,:,-1]
  ncvar = dsout.createVariable( 'ls_forcing_right_{}'.format( names[a] ), 'f4', ( 'time', zs[a], ys[a] ), 
                                fill_value=-9999.0 )
  ncvar[:] = var
  ncvar.units = units[a]
  ncvar.source = "MEPS analysis for 20170609"
  ncvar.long_name = lnames[a]
  ncvar.standard_name = ""
  ncvar.lod = 2


# ls_forcing_north_XX (time,z/zw,x/xu)
lnames = ["large-scale forcing for north model boundary for the wind component in x direction",\
          "large-scale forcing for north model boundary for the wind component in y direction",\
          "large-scale forcing for north model boundary for the wind component in z direction",\
          "large-scale forcing for north model boundary for the potential temperature",\
          "large-scale forcing for north model boundary for the specific humidity"]

for a in range( len( arrays ) ):
  var   = arrays[a][:,:,-1,:]
  ncvar = dsout.createVariable( 'ls_forcing_north_{}'.format( names[a] ), 'f4', ( 'time', zs[a], xs[a] ), 
                                fill_value=-9999.0 )
  ncvar[:] = var
  ncvar.units = units[a]
  ncvar.source = "MEPS analysis for 20170609"
  ncvar.long_name = lnames[a]
  ncvar.standard_name = ""
  ncvar.lod = 2


# ls_forcing_south_XX (time,z/zw,x/xu)
lnames = ["large-scale forcing for south model boundary for the wind component in x direction",\
          "large-scale forcing for south model boundary for the wind component in y direction",\
          "large-scale forcing for south model boundary for the wind component in z direction",\
          "large-scale forcing for south model boundary for the potential temperature",\
          "large-scale forcing for south model boundary for the specific humidity"]

for a in range( len( arrays ) ):
  var   = arrays[a][:,:,0,:]
  ncvar = dsout.createVariable( 'ls_forcing_south_{}'.format( names[a] ), 'f4', ( 'time', zs[a], xs[a] ), 
                                fill_value=-9999.0 )
  ncvar[:] = var
  ncvar.units = units[a]
  ncvar.source = "MEPS analysis for 20170609"
  ncvar.long_name = lnames[a]
  ncvar.standard_name = ""
  ncvar.lod = 2


# ls_forcing_top_XX (time,y/yv,x/xu)
lnames = ["large-scale forcing for top model boundary for the wind component in x direction",\
          "large-scale forcing for top model boundary for the wind component in y direction",\
          "large-scale forcing for top model boundary for the wind component in z direction",\
          "large-scale forcing for top model boundary for the potential temperature",\
          "large-scale forcing for top model boundary for the specific humidity"]

for a in range( len( arrays ) ):
  var   = arrays[a][:,-1,:,:]
  ncvar = dsout.createVariable( 'ls_forcing_top_{}'.format( names[a] ), 'f4', ( 'time', ys[a], xs[a] ), 
                                fill_value=-9999.0 )
  ncvar[:] = var
  ncvar.units = units[a]
  ncvar.source = "MEPS analysis for 20170609"
  ncvar.long_name = lnames[a]
  ncvar.standard_name = ""
  ncvar.lod = 2

del arrays

# init_soil_t
ncvar = dsout.createVariable( 'init_soil_t', 'f4', ( 'zsoil', 'y', 'x' ), 
                              fill_value=-9999.0 )
ncvar[:] = init_soil_t
ncvar.units = "K"
ncvar.source = "MEPS analysis for 20170609"
ncvar.long_name = "initial soil temperature"
ncvar.standard_name = ""
ncvar.lod = 2

# init_soil_m
ncvar = dsout.createVariable( 'init_soil_m', 'f4', ( 'zsoil', 'y', 'x' ), 
                              fill_value=-9999.0 )
ncvar[:] = init_soil_m
ncvar.units = "m^3/m^3"
ncvar.source = "MEPS analysis for 20170609"
ncvar.long_name = "initial soil moisture"
ncvar.standard_name = ""
ncvar.lod = 2


#%% Air pollutants for the child domain

# Provide this information:
maxspec = 7
ls_names = ['left','right','north','south']
ls_dimsname = ['y','y','x','x']

# -------------------------------------------------------------------------------------------------#
def integrate_profile( old_z, old_prof, new_z, new_prof ):

  kk = 0
  for k in np.arange(0,len(new_z),1):

    if ( kk < len(old_z) ):
      while ( old_z[kk+1] <= new_z[k] ):
        kk = kk +1
        if ( kk == len(old_z)-1 ):
          break
    
    if ( kk < (len(old_z)-1) ):
      new_prof[k] = old_prof[kk] + ( new_z[k] - old_z[kk] ) / ( old_z[kk+1] - old_z[kk] ) * \
                    ( old_prof[kk+1] - old_prof[kk] )
    else:
      new_prof[k] = old_prof[kk]
  
  return new_prof

# Read in the file
adchem = adchem = hdf5storage.loadmat( fname_adchem )

dmid = np.squeeze( adchem['dmid'] )
prof_z = np.squeeze( adchem['z'] )

smonth = int( adchem['smonth'][0][0] )
emonth = int( adchem['emonth'][0][0] )
sday   = int( adchem['sday'][0][0] )
eday   = int( adchem['eday'][0][0] )
shour  = int( adchem['shour'][0][0] )
ehour  = int( adchem['ehour'][0][0] )

# Select correct times:
adchem_timestring = np.arange( '{:04d}-{:02d}-{:02d}T{:02d}:00'.format( sim_year, smonth, sday, shour ),
                               '{:04d}-{:02d}-{:02d}T{:02d}:00'.format( sim_year, emonth, eday, ehour ),
                               dtype='datetime64[h]')
for t in range( len( time ) ):
  for tt in range( len( adchem_timestring ) ):
    if timestring[t] == adchem_timestring[tt]:
      if ( timestring[t].hour + plusUTC ) == start_hour:  ### CHECK plusutc FROM PONTUS!!!
        ss = tt
      if ( timestring[t].hour + plusUTC ) == end_hour:  ### CHECK plusutc FROM PONTUS!!! 
        ee = tt

ncc = len( adchem['composition_name'][0] )
composition_index = np.linspace( 1, ncc, ncc ) 
composition_name = []
for cc in range ( ncc ):
  composition_name.append('{:25}'.format( adchem['composition_name'][0][cc][0][0] ) )
max_string_length = np.linspace( 1, len( composition_name[0] ), len( composition_name[0] ) )

# PALM grid:

# PALM model extents:
lx_N03 = (nx_N03+1) * dx_N03
ly_N03 = (ny_N03+1) * dy_N03

z_column = np.zeros( nz_N03+1, dtype=float )
z_column[0] = dz_N03 * 0.5
for k in range( 1, len(z_column), 1):
  z_column[k] = z_column[k-1] + dz_N03
z_N03 = z_column[0:-1]
  
x_N03 = np.linspace( 0.5*dx_N03, lx_N03 - 0.5*dx_N03, nx_N03+1)
y_N03 = np.linspace( 0.5*dy_N03, ly_N03 - 0.5*dy_N03, ny_N03+1)


# Save dimensions to the dynamic input file:

dims = ['x','y','z','time','Dmid','composition_index','max_string_length']
arrays = [x_N03, y_N03, z_N03, time_palm, dmid, composition_index, max_string_length]
std_names = ['x coordinate of cell centers', 'y coordinate of cell centers', 
             'z coordinate of cell centers', 'time', 'aerosol_geometric_mean_diameter','','']
units = ['m','m','m',
         'seconds since {:04d}{:02d}{:02d} {}:00 UTC'.format( sim_year, sim_month, sim_day, start_hour),
         'm','']
types = ['f4','f4','f4','f4','f4','i4','i4']

for d in range( len( dims ) ):
  dsout_N03.createDimension( dims[d], len( arrays[d] ) )
  v = dsout_N03.createVariable( dims[d], 'f4', ( dims[d],) )
  v[:] = arrays[d]
  v.standard_name = std_names[d]
  if dims[d]=='time':
    v.long_name = "time"
    
# -----------------------------------------------------------------------------------------------#   
# Initial profiles: 

# aerosol chemical composition:                      
init_adchem_mf_a = adchem['mass_fracs'][ss,:,:]
init_mf_a = np.zeros( [ len( z_N03 ), ncc ], dtype=float  )
for c in range( ncc ):
  init_mf_a[:,c] = integrate_profile( prof_z, init_adchem_mf_a[:,c], z_N03, init_mf_a[:,c] )
init_mf_b = 0.0 * init_mf_a

# aerosol size distribution:
init_adchem_psd = adchem['psd'][ss,:,:]
nbins = np.shape( init_adchem_psd )[1]
init_psd = np.zeros( [ len( z_N03 ), nbins ], dtype=float  )
for b in range( nbins ):
  init_psd[:,b] = integrate_profile( prof_z, init_adchem_psd[:,b], z_N03, init_psd[:,b] )

# -----------------------------------------------------------------------------------------------#   
# Forcing:
  
ls_dims = [ y_N03, y_N03, x_N03, x_N03]

# aerosol chemical composition:  
lsf_adchem_mf_a = adchem['mass_fracs'][ss:ee+1,:,:]
lsf_mf_a = np.zeros( [ ee-ss+1, len( z_N03 ), ncc ], dtype=float  )
for t in range( ee-ss+1 ):
  for c in range( ncc ):
    lsf_mf_a[t,:,c] = integrate_profile( prof_z, lsf_adchem_mf_a[t,:,c], z_N03, lsf_mf_a[t,:,c] )
lsf_mf_b = 0.0 * lsf_mf_a

# aerosol size distribution:
lsf_adchem_psd = adchem['psd'][ss:ee+1,:,:]
lsf_psd = np.zeros( [ ee-ss+1, len( z_N03 ), nbins ], dtype=float  )
for t in range( ee-ss+1 ):
  for b in range( nbins ):
    lsf_psd[t,:,b] = integrate_profile( prof_z, lsf_adchem_psd[t,:,b], z_N03, lsf_psd[t,:,b] )

# -----------------------------------------------------------------------------------------------# 
# SAVE VARIABLES TO PIDS_DYNAMIC:

# namelist of chemical components:                            
cn = dsout_N03.createVariable( 'composition_name', 'S1', ('composition_index','max_string_length',) )
cn[:] = list( map( lambda x : list(x), composition_name ) )
cn.long_name = 'aerosol composition name'
cn.standard_name = 'composition_name' 
                         
# aerosol chemical composition:  
init_mf_a_v = createNetcdfVariable( dsout_N03, init_mf_a, 'init_atmosphere_mass_fracs_a', 0,  
                                    '', 'f4', ('z','composition_index',), variable, fill_value=-9999.0 )
init_mf_a_v.long_name = "initial mass fraction profile: a bins" 
                                       
init_mf_b_v = createNetcdfVariable( dsout_N03, init_mf_b, 'init_atmosphere_mass_fracs_b', 0, 
                                    '', 'f4', ('z','composition_index',), variable, fill_value=-9999.0  )                                          
init_mf_b_v.long_name = "initial mass fraction profile: b bins"

# forcing for aerosol chemical composition:
nvi = 0
for nv in ls_names:
  dummy = np.zeros( [ ee-ss+1, len(z_N03), len(ls_dims[nvi]), len(composition_index) ] )
  for i in range( len( ls_dims[nvi] ) ):
    dummy[:,:,i,:] = lsf_mf_a
  namev = 'ls_forcing_{}_mass_fracs_a'.format(nv)
  lsf_mf_a_v = createNetcdfVariable( dsout_N03, dummy, namev, 0, '', 'f4', 
                                     ('time','z',ls_dimsname[nvi], 'composition_index',), 
                                     variable, fill_value=-9999.0 )
  lsf_mf_a_v.long_name = "boundary conditions of mass fraction profile: a bins" 
      
  namev = 'ls_forcing_{}_mass_fracs_b'.format(nv)                                 
  lsf_mf_b_v = createNetcdfVariable( dsout_N03, dummy*0, namev, 0, '', 'f4', 
                                     ('time','z',ls_dimsname[nvi],'composition_index',), 
                                     variable, fill_value=-9999.0  )                                          
  lsf_mf_b_v.long_name = "boundary conditions of mass fraction profile: b bins"  
  
  nvi += 1

dummy = np.zeros( [ ee-ss+1, len( y_N03 ), len( x_N03 ), len(composition_index) ] )
for j in range( len( y_N03 ) ):
  for i in range( len( x_N03 ) ):
    dummy[:,j,i,:] = lsf_mf_a[:,-1,:]
namev = 'ls_forcing_top_mass_fracs_a'
lsf_mf_a_v = createNetcdfVariable( dsout_N03, dummy, namev, 0, '', 'f4', 
                                  ('time','y','x','composition_index',), 
                                  variable, fill_value=-9999.0 )
lsf_mf_a_v.long_name = "boundary conditions of mass fraction profile: a bins" 
    
namev = 'ls_forcing_top_mass_fracs_b'                               
lsf_mf_b_v = createNetcdfVariable( dsout_N03, dummy*0, namev, 0, '', 'f4', 
                                  ('time','y','x','composition_index',), 
                                  variable, fill_value=-9999.0  )
lsf_mf_b_v.long_name = "boundary conditions of mass fraction profile: b bins"  

                            
# aerosol size distribution                                           
init_psdv = createNetcdfVariable( dsout_N03, init_psd, 'init_atmosphere_aerosol', 0, '#/m3', 
                                  'f4', ('z','Dmid',), variable, fill_value=-9999.0 )   
init_psdv.long_name = 'initial vertical profile of aerosol concentration' 
init_psdv.lod = 1

nvi = 0
for nv in ls_names:
  dummy = np.zeros( [ ee-ss+1, len(z_N03), len(ls_dims[nvi]), nbins ] )
  for i in range( len( ls_dims[nvi] ) ):
    dummy[:,:,i,:] = lsf_psd
  namev = 'ls_forcing_{}_aerosol'.format( nv )
  lsf_psdv = createNetcdfVariable( dsout_N03, dummy, namev, 0, '', 'f4', 
                                   ('time','z',ls_dimsname[nvi],'Dmid',), 
                                   variable, fill_value=-9999.0 )   
  lsf_psdv.long_name = 'boundary condition of aerosol concentration' 
  lsf_psdv.lod = 1
  nvi += 1
  
dummy =  np.zeros( [ ee-ss+1, len( y_N03 ), len( x_N03 ), nbins ], dtype=float ) 
for j in range( len( y_N03 ) ):
  for i in range( len( x_N03 ) ):
    dummy[:,j,i,:] = lsf_psd[:,-1,:]
namev = 'ls_forcing_top_aerosol'
lsf_psdv = createNetcdfVariable( dsout_N03, dummy, namev, 0, '', 'f4', 
                                 ('time','y','x','Dmid',), variable, fill_value=-9999.0 )   
lsf_psdv.long_name = 'boundary condition of aerosol concentration' 
lsf_psdv.lod = 1

# -----------------------------------------------------------------------------------------------# 
# Initial profiles:  

# H2SO4:
init_achem_h2so4 = adchem['H2SO4'][ss,:] * 1e6 # from 1/cm3 to 1/m3
init_h2so4 = np.zeros( len(z_N03), dtype=float )  
init_h2so4 = integrate_profile( prof_z, init_achem_h2so4, z_N03, init_h2so4 ) 

# HNO3:
init_achem_hno3 = adchem['HNO3'][ss,:] * 1e6 # from 1/cm3 to 1/m3
init_hno3 = np.zeros( len(z_N03), dtype=float )  
init_hno3 = integrate_profile( prof_z, init_achem_hno3, z_N03, init_hno3 ) 

# NH3:
init_achem_nh3 = adchem['NH3'][ss,:] * 1e6 # from 1/cm3 to 1/m3
init_nh3 = np.zeros( len(z_N03), dtype=float )  
init_nh3 = integrate_profile( prof_z, init_achem_nh3, z_N03, init_nh3 )

# OCSV:
init_achem_ocsv = adchem['OCSV'][ss,:] * 1e6 # from 1/cm3 to 1/m3
init_ocsv = np.zeros( len(z_N03), dtype=float )  
init_ocsv = integrate_profile( prof_z, init_achem_ocsv, z_N03, init_ocsv )

# OCNV:
init_achem_ocnv = adchem['OCNV'][ss,:] * 1e6 # from 1/cm3 to 1/m3
init_ocnv = np.zeros( len(z_N03), dtype=float )  
init_ocnv = integrate_profile( prof_z, init_achem_ocnv, z_N03, init_ocnv )

# -----------------------------------------------------------------------------------------------# 
# Forcing:

# H2SO4:
lsf_adchem_h2so4 = adchem['H2SO4'][ss:ee+1,:] * 1e6 # from 1/cm3 to 1/m3
lsf_h2so4 = np.zeros( [ ee-ss+1, len(z_N03) ], dtype=float  )
for t in range( ee-ss+1 ):
  lsf_h2so4[t,:]  = integrate_profile( prof_z, lsf_adchem_h2so4[t,:], z_N03, lsf_h2so4[t,:] )

# HNO3:
lsf_adchem_hno3 = adchem['HNO3'][ss:ee+1,:] * 1e6 # from 1/cm3 to 1/m3
lsf_hno3 = np.zeros( [ ee-ss+1, len(z_N03) ], dtype=float  )
for t in range( ee-ss+1 ):
  lsf_hno3[t,:]  = integrate_profile( prof_z, lsf_adchem_hno3[t,:], z_N03, lsf_hno3[t,:] )

# NH3:
lsf_adchem_nh3 = adchem['NH3'][ss:ee+1,:] * 1e6 # from 1/cm3 to 1/m3
lsf_nh3 = np.zeros( [ ee-ss+1, len(z_N03) ], dtype=float  )
for t in range( ee-ss+1 ):
  lsf_nh3[t,:]  = integrate_profile( prof_z, lsf_adchem_nh3[t,:], z_N03, lsf_nh3[t,:] )

# OCSV:
lsf_adchem_ocsv = adchem['OCSV'][ss:ee+1,:] * 1e6 # from 1/cm3 to 1/m3
lsf_ocsv = np.zeros( [ ee-ss+1, len(z_N03) ], dtype=float  )
for t in range( ee-ss+1 ):
  lsf_ocsv[t,:]  = integrate_profile( prof_z, lsf_adchem_ocsv[t,:], z_N03, lsf_ocsv[t,:] )

# OCNV:
lsf_adchem_ocnv = adchem['OCNV'][ss:ee+1,:] * 1e6 # from 1/cm3 to 1/m3
lsf_ocnv = np.zeros( [ ee-ss+1, len(z_N03) ], dtype=float  )
for t in range( ee-ss+1 ):
  lsf_ocnv[t,:]  = integrate_profile( prof_z, lsf_adchem_ocnv[t,:], z_N03, lsf_ocnv[t,:] )

# -----------------------------------------------------------------------------------------------# 
# SAVE VARIABLES TO PIDS_DYNAMIC:  

# H2SO4: 
h2so4_adchem = createNetcdfVariable( dsout_N03, init_h2so4, 'init_atmosphere_h2so4', len(z_N03), 
                                     '#/m3', 'f4', ('z',), variable, fill_value=-9999.0  )   
h2so4_adchem.long_name = "initial vertical profile of H2SO4" 

nvi = 0
for nv in ls_names:
  dummy = np.zeros( [ ee-ss+1, len(z_N03), len(ls_dims[nvi])] )
  for i in range( len( ls_dims[nvi] ) ):
    dummy[:,:,i] = lsf_h2so4
  namev = 'ls_forcing_{}_h2so4'.format( nv )
  lsf_h2so4v = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4', 
                                     ('time','z',ls_dimsname[nvi],), 
                                     variable, fill_value=-9999.0  )   
  lsf_h2so4v.long_name = "background concentration of H2SO4"   
  nvi += 1
  
dummy = np.zeros( [ ee-ss+1, len( y_N03 ), len( x_N03 ) ] )  
for j in range( len( y_N03 ) ):
  for i in range( len( x_N03 ) ):
    dummy[:,j,i] = lsf_h2so4[:,-1]
namev = 'ls_forcing_top_h2so4'
lsf_h2so4v = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4', ('time','y','x',), 
                                   variable, fill_value=-9999.0  )   
lsf_h2so4v.long_name = "background concentration of H2SO4"


# HNO3:
hno3_adchem = createNetcdfVariable( dsout_N03, init_hno3, 'init_atmosphere_hno3', len(z_N03), 
                                    '#/m3', 'f4', ('z',), variable, fill_value=-9999.0  )   
hno3_adchem.long_name = "initial vertical profile of HNO3"    

nvi = 0
for nv in ls_names:
  dummy = np.zeros( [ ee-ss+1, len(z_N03), len(ls_dims[nvi])] )
  for i in range( len( ls_dims[nvi] ) ):
    dummy[:,:,i] = lsf_hno3
  namev = 'ls_forcing_{}_hno3'.format( nv )
  lsf_hno3v = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4', 
                                   ('time','z',ls_dimsname[nvi],), 
                                    variable, fill_value=-9999.0  )   
  lsf_hno3v.long_name = "background concentration of HNO3" 
  nvi += 1
  
dummy = np.zeros( [ ee-ss+1, len( y_N03 ), len( x_N03 ) ] )  
for j in range( len( y_N03 ) ):
  for i in range( len( x_N03 ) ):
    dummy[:,j,i] = lsf_hno3[:,-1]  
namev = 'ls_forcing_top_hno3'
lsf_hno3v = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4', ('time','y','x',),
                                  variable, fill_value=-9999.0  )   
lsf_hno3v.long_name = "background concentration of HNO3"


# NH3:
nh3_adchem = createNetcdfVariable( dsout_N03, init_nh3, 'init_atmosphere_nh3', len(z_N03), 
                                   '#/m3', 'f4', ('z',), variable, fill_value=-9999.0  )   
nh3_adchem.long_name = "initial vertical profile of NH3"

nvi = 0
for nv in ls_names:
  dummy = np.zeros( [ ee-ss+1, len(z_N03), len(ls_dims[nvi])] )
  for i in range( len( ls_dims[nvi] ) ):
    dummy[:,:,i] = lsf_nh3
  namev = 'ls_forcing_{}_nh3'.format( nv )
  lsf_nh3v = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4', 
                                  ('time','z',ls_dimsname[nvi],), 
                                   variable, fill_value=-9999.0  )   
  lsf_nh3v.long_name = "background concentration of NH3"
  nvi += 1
  
dummy = np.zeros( [ ee-ss+1, len( y_N03 ), len( x_N03 ) ] )  
for j in range( len( y_N03 ) ):
  for i in range( len( x_N03 ) ):
    dummy[:,j,i] = lsf_nh3[:,-1]
namev = 'ls_forcing_top_nh3'
lsf_nh3v = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4', ('time','y','x',), 
                                 variable, fill_value=-9999.0  )   
lsf_nh3v.long_name = "background concentration of NH3"


# OCSV:
ocsv_adchem = createNetcdfVariable( dsout_N03, init_ocsv, 'init_atmosphere_ocsv', len(z_N03), 
                                    '#/m3', 'f4', ('z',), variable, fill_value=-9999.0  )   
ocsv_adchem.long_name = "initial vertical profile of OCSV"

nvi = 0
for nv in ls_names:
  dummy = np.zeros( [ ee-ss+1, len(z_N03), len(ls_dims[nvi])] )
  for i in range( len( ls_dims[nvi] ) ):
    dummy[:,:,i] = lsf_ocsv
  namev = 'ls_forcing_{}_ocsv'.format( nv )
  lsf_ocsvv = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4',
                                   ('time','z',ls_dimsname[nvi],), 
                                    variable, fill_value=-9999.0  )   
  lsf_ocsvv.long_name = "background concentration of OCSV"   
  nvi += 1
  
dummy = np.zeros( [ ee-ss+1, len( y_N03 ), len( x_N03 ) ] )  
for j in range( len( y_N03 ) ):
  for i in range( len( x_N03 ) ):
    dummy[:,j,i] = lsf_ocsv[:,-1]
namev = 'ls_forcing_top_ocsv'
lsf_ocsvv = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4', ('time','y','x',), 
                                  variable, fill_value=-9999.0  )   
lsf_ocsvv.long_name = "background concentration of OCSV"


# OCNV:
ocnv_adchem = createNetcdfVariable( dsout_N03, init_ocnv, 'init_atmosphere_ocnv', len(z_N03), 
                                    '#/m3', 'f4', ('z',), variable, fill_value=-9999.0  )   
ocnv_adchem.long_name = "initial vertical profile of OCNV"  

nvi = 0
for nv in ls_names:
  dummy = np.zeros( [ ee-ss+1, len(z_N03), len(ls_dims[nvi])] )
  for i in range( len( ls_dims[nvi] ) ):
    dummy[:,:,i] = lsf_ocnv
  namev = 'ls_forcing_{}_ocnv'.format( nv )
  lsf_ocnvv = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4', 
                                   ('time','z',ls_dimsname[nvi],),
                                    variable, fill_value=-9999.0  )   
  lsf_ocnvv.long_name = "background concentration of OCNV"  
  nvi += 1
  
dummy = np.zeros( [ ee-ss+1, len( y_N03 ), len( x_N03 ) ] )  
for j in range( len( y_N03 ) ):
  for i in range( len( x_N03 ) ):
    dummy[:,j,i] = lsf_ocnv[:,-1]
namev = 'ls_forcing_top_ocnv'
lsf_ocnvv = createNetcdfVariable( dsout_N03, dummy, namev, 0, '#/m3', 'f4', ('time','y','x',), 
                                  variable, fill_value=-9999.0  )   
lsf_ocnvv.long_name = "background concentration of OCNV"


# gases:
# NO --> cnc_NO_gas
# NO2 --> cnc_NO2_gas
# O3 --> cnc_O3_gas
# HO2 --> cnc_HO2_gas
# OH --> cnc_OH_gas

# alkanes (RH) !!!!
# alkyl peroxy radical (RO2) !!!!
# carbonyl (RCHO)  !!!!

# HNO3 (g) --> cnc_HNO3_gas
# NH3 (g) --> cnc_NH3_gas
# H2SO4 (g) --> cnc_H2SO4_gas
# OCNV !!!!
# OCSV !!!!

# aerosol:
# SO4 --> cnc_SO4_m_20 and cnc_SO4_m_70
# OC !!!!
# BC --> cnc_EC_m_05, m_50, m1_8, m6_0 (note: EC, not BC)
# NH --> cnc_NH415SO4_m_20, cnc_NH415SO4_m_70, cnc_NH4NO3_m_70
# NO --> cnc_NO3_c_m3_0
# SS --> cnc_sslt_m_05, m_50, m3_0, m9_0
# DU --> cnc_mineral_m_05, m_50, m1_8, m6_0


#%% Close files

dsout.close()
dsout_N03.close()
