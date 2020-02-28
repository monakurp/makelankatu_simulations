import numpy as np
import netCDF4 as nc
import os
import datetime as dtime
from psdLib import define_bins
from scipy.interpolate import interp1d, CubicSpline
import matplotlib.pyplot as plt
import matplotlib as mpl
import warnings
warnings.filterwarnings( 'ignore', category=RuntimeWarning )

os.chdir(os.path.expanduser("~"))
os.chdir(('makelankatu_simulations'))
np.set_printoptions( precision=10 )

parameter = True
variable  = False

def_fs = 9

mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.labelsize'] = def_fs
mpl.rcParams['ytick.labelsize'] = def_fs
mpl.rcParams['xtick.labelsize'] = def_fs 
mpl.rcParams['grid.linewidth'] = 1
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['legend.fontsize'] = def_fs-1
mpl.rcParams['axes.titlesize'] = def_fs
mpl.rc( 'font', size=def_fs )

# ------------------------------------------------------------------------------------------------#
# PROVIDE THESE:
# ------------------------------------------------------------------------------------------------#

# Simulation time: start of drone measurements - ~15min --> end of drone measurements
sim_year    = 2017
sim_month   = 6 # 12
sim_day     = 9 # 7
sim_time    = 'morning'
datestr     = '{}{:02d}{:02d}'.format( sim_year, sim_month, sim_day )

precursor   = True

# ------------------------------------------------------------------------------------------------#

if ( datestr=='20170609' and sim_time=='morning' ):
  start_hour = 7
  start_min  = 0
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
  
if precursor:
  end_hour = start_hour
  start_hour -= 1
  
#%% Filenames

pids_folder = 'input_data_to_palm/cases/{}_{}'.format( datestr, sim_time )
harmonie_folder = 'source_data/cases/{}_{}/'.format( datestr, sim_time )

if precursor:
  suffix = '_precursor'
else:
  suffix = ''

fname_in = '{}/PIDS_DYNAMIC{}'.format( pids_folder, suffix )
fname_in_N03 = '{}/PIDS_DYNAMIC_N03_FRES'.format( pids_folder)

fname_out_wd  = '{}/PIDS_DYNAMIC_wd_kivenlahti{}'.format( pids_folder, suffix  ) 
fname_out_all = '{}/PIDS_DYNAMIC_all_kivenlahti{}'.format( pids_folder, suffix  )
fname_out_N03 = '{}/PIDS_DYNAMIC_N03_smear'.format( pids_folder )

fname_kivenlahti = 'source_data/kivenlahti/kivenlahti_mast_{}{}{}.txt'.format( sim_day, sim_month,
                                                                               sim_year )
fname_smear_met = 'source_data/smear/smeardata_{}{:02d}.txt'.format( sim_year, sim_month )
fname_smear_dmps = 'source_data/smear/dmps/dp{}{:02d}{:02d}.sum'.format( sim_year-2000, sim_month,
                                                                         sim_day )

fname_full_backup = '{}/meps_mbr0_full_backup_2_5km_{}T00Z.nc'.format( harmonie_folder, datestr )
fname_subset      = '{}/meps_subset_2_5km_{}T00Z.nc'.format( harmonie_folder, datestr )

#%% Copy the original PIDS_DYNAMIC files
  
# Input files:
dsin = nc.Dataset( fname_in )
if not precursor:
  dsin_N03 = nc.Dataset( fname_in_N03 )

# Output files:
dsout_wd  = nc.Dataset( fname_out_wd, "w" )
dsout_all = nc.Dataset( fname_out_all, "w" )
if not precursor:
  dsout_N03 = nc.Dataset( fname_out_N03, "w" )
 
d = 0
if not precursor:
  datas = [dsout_wd, dsout_all, dsout_N03]
  namesout = [fname_out_wd, fname_out_all, fname_out_N03]
else:
  datas = [dsout_wd, dsout_all]
  namesout = [fname_out_wd, fname_out_all]

for dst in datas:
  if d<2:
    src = dsin
    namein = fname_in
  else:
    src = dsin_N03
    namein = fname_in_N03
  print( 'Copy {} to {}'.format( namein, namesout[d] ) )

  # Copy attributes
  dst.setncatts( src.__dict__ )

  # Copy dimensions
  for dname, the_dim in src.dimensions.items():
    dst.createDimension( dname, len( the_dim ) if not the_dim.isunlimited() else None ) 

  # Copy variables
  for v_name, varin in src.variables.items():
      outVar = dst.createVariable( v_name, varin.datatype, varin.dimensions )
      # Copy variable attributes
      outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
      outVar[:] = varin[:]  

  d += 1 

# close the input files
dsin.close()
if not precursor:
  dsin_N03.close()
  dsout_N03.comment = 'air quality data from https://avaa.tdata.fi/web/smart/smearfrom' 

# ------------------------------------------------------------------------------------------------#

# Set global attributes:
dsout_wd.comment = 'Wind direction data retrieved from https://en.ilmatieteenlaitos.fi/open-data'
dsout_all.comment = 'Meteorological data retrieved from https://en.ilmatieteenlaitos.fi/open-data'

for dsouti in datas:
  dsouti.title = "PALM dynamic input file for scenario Makelankatu, {} {:02d}/{:02d}/{}".format(\
                 sim_time, sim_day, sim_month, sim_year)

ntimesteps = len( dsout_wd['time'] )

del datas 

#%% Read in Harmonie data

orig = [60.1663312, 24.873065935]
radius_earth = 6371229.0 # m
npoints_interpolation = 4

lats = np.zeros( [npoints_interpolation], dtype=float ); lats[0] = orig[0]
lons = np.zeros( [npoints_interpolation], dtype=float ); lons[0] = orig[1]

full_backup = nc.Dataset( fname_full_backup )
subset = nc.Dataset( fname_subset )

# Check also that the files are for the same time and coordinates!

time = np.array( full_backup['time'] ) # "seconds since 1970-01-01 00:00:00 +00:00"
timestring = dict()
for t in range( len( time ) ):
  timestring[t] = dtime.datetime.utcfromtimestamp( time[t] )
  if ( timestring[t].hour + plusUTC ) == start_hour:
    ss = t
  if ( timestring[t].hour + plusUTC ) == end_hour:
    ee = t
if end_min > 0:
  ee += 1
    
longitude = np.array( full_backup['longitude'] )
longitude = np.reshape( longitude, ( len(lats),len(lons) ) )

latitude = np.array( full_backup['latitude'] )
latitude = np.reshape( latitude, ( len(lats),len(lons) ) )

ensemble = np.array( full_backup['ensemble_member'] )

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

# Define the pressure and height levels

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

# Potential temperature

pt = np.zeros( [ ntimesteps, len(hybrid), len(lats), len(lons) ] ) + np.nan
for t in range( ntimesteps ):
  for h in range( len( hybrid ) ):
    pt[t,h,:,:] = T[t,h,:,:] * ( 100000.0 / p[t,h,:,:] ) ** 0.286

#%% Read in SMEARIII met data

# 0 = Year
# 1 = Month
# 2 = Day
# 3 = Hour
# 4 = Minute
# 5 = Second
# 6 = KUM_META.t
# 8 = KUM_META.rh
# 9 = KUM_META.p
# 10 = KUM_META.ws
# 11 = KUM_META.wdir
# 12 = KUM_META.Tower_T_4m
# 13 = KUM_META.Tower_T_8m
# 14 = KUM_META.Tower_T_16m
# 15 = KUM_META.Tower_T_32m
# 16 = KUM_META.Tower_WS_32m
# 17 = KUM_META.Tower_WS_16m
# 18 = KUM_META.Tower_WS_8m
# 19 = KUM_META.Tower_WS_4m
# 20 = KUM_META.Tower_WDIR_4m
# 21 = KUM_META.Tower_WDIR_8m
# 22 = KUM_META.Tower_WDIR_16m
# 23 = KUM_META.Tower_WDIR_32m
# 24 = KUM_META.Tower_VT_32m
# 25 = KUM_META.Tower_VT_16m
# 26 = KUM_META.Tower_VT_8m
# 27 = KUM_META.Tower_VT_4m

kumpula = np.loadtxt( fname_smear_met )
if sim_month==6:
  kumpula[:,3] += 1
  kumpula[kumpula[:,3]>23,3] = 0

pressure_smear = np.zeros( ntimesteps )

h = start_hour
t = 0
while t < ntimesteps:
  ind = ( kumpula[:,2]==sim_day ) & ( ( ( kumpula[:,3] == h ) & ( kumpula[:,4] < 10 ) ) | \
                                      ( ( kumpula[:,3] == h-1 ) & ( kumpula[:,4] > 49 ) ) )
  pressure_smear[t] = np.nanmean( kumpula[ind,9] )
  h += 1
  t += 1

#%% Read in Kivenlahti data

# 0 = year
# 1 = month
# 2 = day
# 3 = hour(utc+3)
# 4 = minute
# 5 = T_2
# 6 = TD_2
# 7 = U_2
# 8 = WD_2
# 9 = T_26
# 10 = TD_26
# 11 = U_26
# 12 = WD_26
# 13 = T_49
# 14 = TD_49
# 15 = U_49
# 16 = WD_49
# 17 = T_92
# 18 = TD_92
# 19 = U_92
# 20 = WD_92
# 21 = T_141
# 22 = TD_141
# 23 = U_141
# 24 = WD_141
# 25 = T_217
# 26 = TD_217
# 27 = U_217
# 28 = WD_217
# 29 = T_265
# 30 = TD_265
# 31 = U_265
# 32 = WD_265
# 33 = T_327
# 34 = TD_327
# 35 = U_327
# 36 = WD_327

kivenlahti = np.loadtxt( fname_kivenlahti, skiprows=1 )
kivenlahti[kivenlahti==-999] = np.nan

uk = np.zeros([ntimesteps,8])
vk = np.copy( uk )
WDk = np.copy( uk )
Uk = np.copy( WDk )
Tk = np.copy( WDk )
TDk = np.copy( WDk )
zk = np.array( [2, 26, 49, 92, 141, 217, 265, 327] ) 
cols = [ 7, 11, 15, 19, 23, 27, 31, 35 ] 

h = start_hour
t = 0
while t < ntimesteps:
  ind = ( kivenlahti[:,3]==h ) & ( ( kivenlahti[:,4]==0 ) )
  for zi in range( 8 ):
    uk[t,zi] = -kivenlahti[ind,cols[zi]] * np.sin( kivenlahti[ind,cols[zi]+1] * np.pi/180.0 )
    vk[t,zi] = -kivenlahti[ind,cols[zi]] * np.cos( kivenlahti[ind,cols[zi]+1] * np.pi/180.0 )
    Uk[t,zi] = np.nanmean( kivenlahti[ind,cols[zi]] )
    WDk[t,zi] = 270.0 - ( np.arctan2( vk[t,zi], uk[t,zi] ) * 180.0 / np.pi )
    Tk[t,zi] = np.nanmean( kivenlahti[ind,cols[zi]-2] )
    TDk[t,zi] = np.nanmean( kivenlahti[ind,cols[zi]-1] )

  t += 1
  h += 1
  
# Check that the dew point temperature does not exceed the temperature
for t in range( ntimesteps ):
  for zi in range( 8 ):
    if Tk[t,zi] < TDk[t,zi]:
      print( 'too high dew point temperature: T={:2.1f}C, TD={:2.1f}C --> T={:2.1f}C'\
             ' (0.98*T)'.format( Tk[t,zi], TDk[t,zi], Tk[t,zi] * 0.98 ) )
      TDk[t,zi] = 0.98*Tk[t,zi]
  
# Dew point temperature to q: assume that p=pressure_smear
def specific_humidity_from_dew( dew, p ):

  # Give temperatures in Kelvin and p in hPa!
  
  e0 = 6.113 # saturation vapor pressure in hPa
  c_water = 5423 # L/R for water in Kelvin
  T0 = 273.15 # Kelvin
  
  # saturation vapor not required, uncomment to calculate it (units in hPa becuase of e0)
  #sat_vapor = e0 * np.exp((self.c * (temp -T0))/(temp * T0)) 
  
  #calculating specific humidity, q directly from dew point temperature
  #using equation 4.24, Pg 96 Practical Meteorolgy (Roland Stull)
  q = ( 622 * e0 * np.exp( c_water * ( dew - T0 ) / ( dew * T0 ) ) ) / p # g/kg 
  # 622 is the ratio of Rd/Rv in g/kg
  return q/1000.0 # in kg/kg
  
qk = np.copy( uk )
for t in range( len( pressure_smear ) ):
  qk[t,:] = specific_humidity_from_dew( TDk[t,:]+273.15, pressure_smear[t] )
  
  
#%% Read SMEARIII DMPS data
  
if not precursor:
  
  dmps = np.genfromtxt( fname_smear_dmps )
  dmid_dmps = dmps[0,2::]
  day_dmps = dmps[:,0]
  ntot_dmps = dmps[:,1] # 1/cm3
  
  # Change the days from the beginning of the year (day_dmps) to datetime (time_dmps)
  for i in range( len( day_dmps ) ):
    
    # days from the beginning of the year to datetime
    timecol = dtime.datetime( sim_year, 1, 1 ) + dtime.timedelta( day_dmps[i] - 1 ) 
    
    # separate to year, month, day, hour, minute, second
    timecol = np.array( [timecol.year, timecol.month, timecol.day, timecol.hour, timecol.minute, 
                             timecol.second] ).T
    if i==0:
      time_dmps = timecol
    else:
      time_dmps = np.vstack(( time_dmps, timecol ))
  
  # Select all measurement times during the simulation time
  ind = (time_dmps[:,1]==sim_month) & (time_dmps[:,2]==sim_day) & \
        ( ( (time_dmps[:,3]>=start_hour) ) & ( (time_dmps[:,3]<=end_hour) ) )
  ind[(time_dmps[:,3]==end_hour) & (time_dmps[:,4]>end_min)] = False
  
  # Calculate the background concentration as the average psd over the simulation time
  dNdlogD_dmps = np.nanmean( dmps[ind,2::], axis=0 )
  
  # Change dNdlogD to N
  dlogdi = np.log10( dmid_dmps[1::]/dmid_dmps[0:-1] )[-1]
  dmidi = np.append( dmid_dmps, dmid_dmps[-1]*10**dlogdi )
  psd_dmps = dNdlogD_dmps * np.log10( dmidi[1::]/dmidi[0:-1] )
  
  # Adjust to the PALM-SALSA grid
  dmid_SALSA = dsout_N03['Dmid'][:].data
  
  _, bin_limits = define_bins( [3, 7], [2.5e-9, 15.0e-9, 1e-6] )
  
  psd_to_palm = np.zeros( len( dmid_SALSA ) )
  for ib in range( len( dmid_SALSA ) ):
    ind = ( dmid_dmps > bin_limits[ib] ) & ( dmid_dmps <= bin_limits[ib+1] );
    psd_to_palm[ib] = np.sum( psd_dmps[ind] )
    
  psd_to_palm *= 1.0e6 # from 1/cm3 to 1/m3
  
#%% Interpolate Kivenlahti measurements to PALM grid

plt.close('all')  

z = dsout_all['z'][:].data
zw = dsout_all['zw'][:].data
u_prof  = np.zeros( [ntimesteps, len( z )] )
v_prof  = np.copy( u_prof )
w_prof  = np.zeros( [ntimesteps, len( zw )] )
pt_prof = np.copy( u_prof )
qv_prof = np.copy( u_prof )

arrays_out = [u_prof, v_prof, w_prof, pt_prof, qv_prof]

iv = 0
for var in [uk, vk, uk*0, Tk, qk]:
  if iv == 2:
    zi = zw
  else:
    zi = z
  
  for t in range( ntimesteps ):
    ind = ~np.isnan( var[t,:] )
    fnearest = interp1d( zk[ind], var[t,ind], kind='nearest', fill_value="extrapolate" )
    flinear = interp1d( zk[ind], var[t,ind], kind='linear', fill_value="extrapolate" )
    fcubic  = CubicSpline( zk[ind], var[t,ind] )
    
    ind2 = zi<=np.max( zk[ind] )
    arrays_out[iv][t,ind2] = flinear( zi[ind2] )
    ind2 = zi>np.max( zk[ind] )
    arrays_out[iv][t,ind2] = fnearest( zi[ind2] )
    
  iv += 1
  
u_prof  = arrays_out[0]
v_prof  = arrays_out[1]
w_prof  = arrays_out[2]
pt_prof = arrays_out[3]
qv_prof = arrays_out[4]


# Plot:
txts = ['$U~\mathrm{(m~s^{-1})}$','$\mathrm{WD~(^\circ)}$']
fig1 = plt.figure(figsize=(12/2.54, 8/2.54), dpi=100  )
fig2 = plt.figure(figsize=(12/2.54, 8/2.54), dpi=100  )
f = 0
for figi in [fig1, fig2]:
  figi.subplots_adjust( hspace=0.4, wspace=0.1, left=0.11, right=0.99, bottom=0.15, top=0.9 ) 
  figi.text( 0.5, 0.02, txts[f], ha='center', fontsize=def_fs )
  f += 1
for t in range(ntimesteps):  
  zi = np.nanmean( np.nanmean( Z[t,:,:,:], axis=-1 ), axis=-1 )
  
  ax1 = fig1.add_subplot(1,ntimesteps,t+1)
  U = np.sqrt( u[t,:,:,:]**2+v[t,:,:,:]**2 )
  minlim = np.percentile( np.percentile( U, 1, axis=-1 ), 1, axis=-1 )
  maxlim = np.percentile( np.percentile( U, 99, axis=-1 ), 99, axis=-1 )
  ax1.fill_betweenx( zi, minlim, maxlim, color='k', alpha=0.3 )
  ax1.plot( np.sqrt( u_prof[t,:]**2 + v_prof[t,:]**2 ), z, '--r*', markevery=10 )
  
  ax2 = fig2.add_subplot(1,ntimesteps,t+1)
  WD = 270.0 - ( np.arctan2( v[t,:,:,:], u[t,:,:,:] ) * 180.0 / np.pi )
 # WD[WD>360] = WD[WD>360]-360
  minlim = np.percentile( np.percentile( WD, 1, axis=-1 ), 1, axis=-1 )
  maxlim = np.percentile( np.percentile( WD, 99, axis=-1 ), 99, axis=-1 )
  ax2.fill_betweenx( zi, minlim, maxlim, color='k', alpha=0.3 )
  WD = 270.0 - ( np.arctan2( v_prof[t,:], u_prof[t,:] ) * 180.0 / np.pi )
  WD[WD>360] = WD[WD>360]-360
  ax2.plot( WD, z, '--r*', markevery=10 )
  
  for axi in [ax1, ax2]:
    axi.set_ylim( 0, 700 )
    axi.text( 0.15, 0.85, '{}:00'.format( start_hour+t ), transform=axi.transAxes )
    if t>0:
      axi.set_yticklabels([])
    else:
      axi.set_ylabel('$z$ (m)', labelpad=0)
  ax1.set_xlim( -1, 5 )
  ax1.set_xticks( np.arange( 0,5,2 ) )
  ax2.set_xlim( 0, 415 )
  ax2.set_xticks( np.arange( 0, 361, 90 ) )
  for figi in [fig1, fig2]:
    figi.legend(['Measurements','Harmonie'], loc='upper center', ncol=2)

#%% Create arrays
    
nx = len( dsout_all['x'][:].data )-1
ny = len( dsout_all['y'][:].data )-1
nz = len( dsout_all['z'][:].data )
    
u_palm  = np.zeros( [ ntimesteps, nz,   ny+1, nx  ], dtype=float )
v_palm  = np.zeros( [ ntimesteps, nz,   ny,   nx+1], dtype=float )
w_palm  = np.zeros( [ ntimesteps, nz-1, ny+1, nx+1], dtype=float )
pt_palm = np.zeros( [ ntimesteps, nz,   ny+1, nx+1], dtype=float )
qv_palm = np.zeros( [ ntimesteps, nz,   ny+1, nx+1], dtype=float )
u2_palm  = np.copy( u_palm )
v2_palm  = np.copy( v_palm )

arrays_out = [u_palm, v_palm, pt_palm, qv_palm]
arrays_in  = [u_prof, v_prof, pt_prof, qv_prof]
zs         = [z,      z,      z,       z      ]
nys        = [ny+1,   ny,     ny+1,    ny+1   ]
nxs        = [nx,     nx+1,   nx+1,    nx+1   ]
names      = ['u',    'v',    'pt',    'qv'   ]

print('Both wind speed and direction:')
for a in range(  len( arrays_out ) ):
  print('   {}'.format(names[a]) )  
  for t in range( ntimesteps ):
    for j in range( nys[a] ):
      for i in range( nxs[a] ):
        arrays_out[a][t,:,j,i] = arrays_in[a][t,:]

u_palm  = arrays_out[0]
v_palm  = arrays_out[1]
pt_palm = arrays_out[2]
qv_palm = arrays_out[3]


arrays_out = [u2_palm, v2_palm]
arrays_in  = [u_prof, v_prof  ]
zs         = [z,      z       ]
nys        = [ny+1,   ny      ]
nxs        = [nx,     nx+1    ]
names      = ['u',    'v'     ]

wdk = 270.0 - ( np.arctan2( v_prof, u_prof ) * 180.0 / np.pi )

# Copy data from PIDS_DYNAMIC
for a in range(  len( arrays_out ) ):
  arrays_out[a][:,:,:,0]  = dsout_wd['ls_forcing_left_{}'.format( names[a] )][:].data
  arrays_out[a][:,:,:,-1] = dsout_wd['ls_forcing_right_{}'.format( names[a] )][:].data
  arrays_out[a][:,:,:,-2] = dsout_wd['ls_forcing_right_{}'.format( names[a] )][:].data
  arrays_out[a][:,:,-1,:] = dsout_wd['ls_forcing_north_{}'.format( names[a] )][:].data
  arrays_out[a][:,:,-2,:] = dsout_wd['ls_forcing_north_{}'.format( names[a] )][:].data
  arrays_out[a][:,:,0,:]  = dsout_wd['ls_forcing_south_{}'.format( names[a] )][:].data
  arrays_out[a][:,-1,:,:] = dsout_wd['ls_forcing_top_{}'.format( names[a] )][:].data

print('Only wind direction: both u & v')
for t in range( ntimesteps ):
  print('   t={}'.format( t ) )
  wdh = wdk[t,:]
  for j in range( ny ):
    for i in range( nx ):
      Uh = np.sqrt( arrays_out[0][t,:,j,i]**2 + arrays_out[1][t,:,j,i]**2 )
      uh2 = -Uh * np.sin( wdh * np.pi/180.0 )
      vh2 = -Uh * np.cos( wdh * np.pi/180.0 )
      arrays_out[0][t,:,j,i] = uh2
      arrays_out[1][t,:,j,i] = vh2

arrays_out[0][:,:,-1,:] = arrays_out[0][:,:,-2,:]
arrays_out[0][:,:,-2,:] = 0
arrays_out[0][:,:,:,-1] = arrays_out[0][:,:,:,-2]
arrays_out[0][:,:,:,-2] = 0
    
u2_palm  = arrays_out[0]
v2_palm  = arrays_out[1]

del arrays_out

#%% Replace meteorological data in PIDS_DYNAMIC with Kivenlahti measurements

names    = ['u',    'v',     'w',    'pt',    'qv'    ]
arrays      = [ u_palm, v_palm,  w_palm, pt_palm, qv_palm]
arrays2     = [ u2_palm, v2_palm ]
arrays_prof = [ u_prof, v_prof,  w_prof, pt_prof, qv_prof]

# init_atmosphere_XX (z/zw)
for a in range( len( arrays ) ):
  dsout_all['init_atmosphere_{}'.format( names[a] )][:] = arrays_prof[a][0,:]
  if a<2:
    dsout_wd['init_atmosphere_{}'.format( names[a] )][:] = arrays_prof[a][0,:]

# ls_forcing_left_XX (time,z/zw,y/yv)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_left_{}'.format(names[a] )][:] = arrays[a][:,:,:,0]
  if a<2:
    dsout_wd['ls_forcing_left_{}'.format(names[a] )][:] = arrays2[a][:,:,:,0]
  
# ls_forcing_right_XX (time,z/zw,y/yv)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_right_{}'.format(names[a] )][:] = arrays[a][:,:,:,-1]
  if a<2:
    dsout_wd['ls_forcing_right_{}'.format(names[a] )][:] = arrays2[a][:,:,:,-1]

# ls_forcing_north_XX (time,z/zw,x/xu)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_north_{}'.format(names[a] )][:] = arrays[a][:,:,-1,:]
  if a<2:
    dsout_wd['ls_forcing_north_{}'.format(names[a] )][:] = arrays2[a][:,:,-1,:]

# ls_forcing_south_XX (time,z/zw,x/xu)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_south_{}'.format(names[a] )][:] = arrays[a][:,:,0,:]
  if a<2:
    dsout_wd['ls_forcing_south_{}'.format(names[a] )][:] = arrays2[a][:,:,0,:]

# ls_forcing_top_XX (time,y/yv,x/xu)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_top_{}'.format(names[a] )][:] = arrays[a][:,-1,:,:]
  if a<2:
    dsout_wd['ls_forcing_top_{}'.format(names[a] )][:] = arrays2[a][:,-1,:,:]

del arrays, arrays2

  
#%% Replace size distribution data in PIDS_DYNAMIC_N03 with SMEARIII DMPS data

if not precursor:
  
  # PS: use the chemical composition from ADCHEM

  for k in range( len( dsout_N03['z'][:] ) ):
    dsout_N03['init_atmosphere_aerosol'][k,:] = psd_to_palm

  bounds = ['north','south','left','right']  
  for bound in bounds:
    shp = np.shape( dsout_N03['ls_forcing_{}_aerosol'.format( bound )] )
    for t in range( ntimesteps ):
      for k in range( shp[-3] ):
        for i in range( shp[-2] ):
          dsout_N03['ls_forcing_{}_aerosol'.format( bound )][t,k,i,:] = psd_to_palm

  bound = 'top'
  shp = np.shape( dsout_N03['ls_forcing_{}_aerosol'.format( bound )] )
  for t in range( ntimesteps ):
    for j in range( shp[-2] ):
      for i in range( shp[-1] ):
        dsout_N03['ls_forcing_{}_aerosol'.format( bound )][t,j,i,:] = psd_to_palm


#%% Close files

dsout_all.close()
dsout_wd.close()
if not precursor:
  dsout_N03.close()