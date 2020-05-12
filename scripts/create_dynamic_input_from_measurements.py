import numpy as np
import netCDF4 as nc
import os
import datetime as dtime
from psdLib import define_bins
from scipy.interpolate import interp1d, CubicSpline
import matplotlib.pyplot as plt
import warnings
import metpy.calc as mpcalc
from metpy.units import units

'''
Description:
Create the dynamic input file for offline nesting in PALM.
For the root:
- using only observed meteorology (PIDS_DYNAMIC_all_kivenlahti) and
- modifying MEPS data with observations (PIDS_DYNAMIC_wd_kivenlahti and PIDS_DYNAMIC_wd_kumpula)
For the child:
- using the measured aerosol size distribution, but all the rest from ADCHEM.

Author: Mona Kurppa
        mona.kurppa@helsinki.fi
        Univeristy of Helsinki
'''

# Ignore runtime warning:
warnings.filterwarnings('ignore', category=RuntimeWarning)

# File path:
os.chdir(os.path.expanduser("~"))
os.chdir(('makelankatu_simulations'))

# Printing precision
np.set_printoptions(precision=10)

parameter = True
variable  = False

# ------------------------------------------------------------------------------------------------#
# PROVIDE THESE:
# ------------------------------------------------------------------------------------------------#

# Simulation time: start of drone measurements - ~15min --> end of drone measurements
sim_year    = 2017
sim_month   = 6
sim_day     = 9
sim_time    = 'morning'
datestr     = '{}{:02d}{:02d}'.format(sim_year, sim_month, sim_day)

precursor   = False

# ------------------------------------------------------------------------------------------------#

if ( datestr=='20170609' and sim_time=='morning' ):
  start_hour = 7
  start_min  = 0
  end_hour   = 9
  end_min    = 15
  plusUTC    = 3
elif ( datestr=='20170609' and sim_time=='evening' ):
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

#%% Constants

rd = 287 # gas constant for dry air, J/(kg*K)
rv = 461.51 # gas constant for moist air, J/(kg*K)

#%% Inline functions

# ------------------------------------------------------------------------------------------------#
# Convert the height to pressure:
def barometric_formula( z, t0, p0 ):

  # Give height in metres, temperature in Kelvin and p in Pa!

  g = 9.81 # m/s2
  cp = 1005.0 # J/kg
  rd = 287.0 # J/(kg*K)
  cp_d_rd = cp/rd

  # Hydrostatic pressure
  hyp = p0 * ( ( t0 - ( g/cp ) * z ) / t0 )**( cp_d_rd )

  return hyp

# ------------------------------------------------------------------------------------------------#
# Exner formula to convert potential temperature to temperature and back
def exner_function( p ):

  p0 = 101325.0 # Pa
  rd = 287.0
  cp = 1005.0
  rd_d_cp = rd/cp
  exner = ( p / p0 )**( rd_d_cp )

  return exner
# ------------------------------------------------------------------------------------------------#

#%% Filenames

# Folders:
pids_folder = 'input_data_to_palm/cases/{}_{}'.format(datestr, sim_time)
harmonie_folder = 'source_data/cases/{}_{}/'.format(datestr, sim_time)

if precursor:
  suffix = '_precursor'
else:
  suffix = ''

# Filenames:

# Existing PIDS_DYNAMIC files created from model data:
fname_in = '{}/PIDS_DYNAMIC{}'.format(pids_folder, suffix)
fname_in_N03 = '{}/PIDS_DYNAMIC_N03_FRES'.format(pids_folder)

# Output files:
fname_out_wd  = '{}/PIDS_DYNAMIC_wd_kivenlahti{}'.format(pids_folder, suffix)
fname_out_all = '{}/PIDS_DYNAMIC_all_kivenlahti{}'.format(pids_folder, suffix)
fname_out_wdk = '{}/PIDS_DYNAMIC_wd_kumpula{}'.format(pids_folder, suffix)
fname_out_N03 = '{}/PIDS_DYNAMIC_N03_smear'.format(pids_folder)

# Kivenlahti meteorological observations:
fname_kivenlahti = 'source_data/kivenlahti/kivenlahti_mast_{}{}{}.txt'.format(sim_day, sim_month,
                                                                              sim_year)
# SMEARIII meteorological observations:
fname_smear_met = 'source_data/smear/smeardata_{}{:02d}.txt'.format(sim_year, sim_month)

# SMEARIII DMPS (aerosol size distribution) data
fname_smear_dmps = 'source_data/smear/dmps/dp{}{:02d}{:02d}.sum'.format(sim_year-2000, sim_month,
                                                                        sim_day)

# MEPS data: full data and a subset for the ensemble member 0
fname_full_backup = '{}/meps_mbr0_full_backup_2_5km_{}T00Z.nc'.format(harmonie_folder, datestr)
fname_subset      = '{}/meps_subset_2_5km_{}T00Z.nc'.format(harmonie_folder, datestr)

#%% Copy the original PIDS_DYNAMIC files

# Input files:
dsin = nc.Dataset(fname_in)
if not precursor:
  dsin_N03 = nc.Dataset(fname_in_N03)

# Output files:
dsout_wd  = nc.Dataset(fname_out_wd, "w") # kivenlahti
dsout_wdk = nc.Dataset(fname_out_wdk, "w") # smeariii (kumpula)
dsout_all = nc.Dataset(fname_out_all, "w")
if not precursor:
  dsout_N03 = nc.Dataset(fname_out_N03, "w")

d = 0
if not precursor:
  datas = [dsout_wd, dsout_all, dsout_wdk, dsout_N03]
  namesout = [fname_out_wd, fname_out_all, fname_out_wdk, fname_out_N03]
else:
  datas = [dsout_wd, dsout_all, dsout_wdk]
  namesout = [fname_out_wd, fname_out_all, fname_out_wdk]

for dst in datas:
  if d<3:
    src = dsin
    namein = fname_in
  else:
    src = dsin_N03
    namein = fname_in_N03
  print('Copy {} to {}'.format(namein, namesout[d]))

  # Copy attributes
  dst.setncatts(src.__dict__)

  # Copy dimensions
  for dname, the_dim in src.dimensions.items():
    dst.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

  # Copy variables
  for v_name, varin in src.variables.items():
      outVar = dst.createVariable(v_name, varin.datatype, varin.dimensions)
      # Copy variable attributes
      outVar.setncatts({k: varin.getncattr(k) for k in varin.ncattrs()})
      outVar[:] = varin[:]

  d += 1 

# close the input files
dsin.close()
if not precursor:
  dsin_N03.close()
  dsout_N03.comment = 'air quality data from https://avaa.tdata.fi/web/smart/smear'

# ------------------------------------------------------------------------------------------------#

# Set global attributes:
dsout_wd.comment = 'Wind direction data retrieved from https://en.ilmatieteenlaitos.fi/open-data'
dsout_all.comment = 'Meteorological data retrieved from https://en.ilmatieteenlaitos.fi/open-data'

for dsouti in datas:
  dsouti.title = "PALM dynamic input file for scenario Makelankatu, {} {:02d}/{:02d}/{}".format(\
                 sim_time, sim_day, sim_month, sim_year)

ntimesteps = len(dsout_wd['time'])

del datas

#%% Read in Harmonie data

orig = [60.1663312, 24.873065935]
radius_earth = 6371229.0 # m
npoints_interpolation = 4

lats = np.zeros( [npoints_interpolation], dtype=float ); lats[0] = orig[0]
lons = np.zeros( [npoints_interpolation], dtype=float ); lons[0] = orig[1]

full_backup = nc.Dataset(fname_full_backup)
subset = nc.Dataset(fname_subset)

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
T_avg = 0.5 * (T[:,0:-1,:,:] + T[:,1::,:,:])

# geopotential height
Z = np.zeros([ntimesteps, len(hybrid), len(lats), len(lons)]) + np.nan

for t in range( ntimesteps ):
  for h in range( len( hybrid ) - 1 ):
    if h==0:
      Z[t,h,:,:] = Z_0[t,:,:]/g0
    else:
      Z[t,h,:,:] = Z[t,h-1,:,:] + T_avg[t,h-1,:,:] * R/g0 * np.log(p[t,h-1,:,:]/p[t,h,:,:])
Z[:,-1,:,:]= Z[:,-2,:,:]

# Potential temperature

pt = np.zeros([ntimesteps, len(hybrid), len(lats), len(lons)]) + np.nan
for t in range(ntimesteps):
  for h in range(len(hybrid)):
    pt[t,h,:,:] = T[t,h,:,:] * (100000.0 / p[t,h,:,:]) ** 0.286

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

kumpula = np.loadtxt(fname_smear_met)
if sim_month==6:
  kumpula[:,3] += 1
  kumpula[kumpula[:,3]>23,3] = 0

pressure_smear = np.zeros(ntimesteps)
wd_smear = np.zeros(ntimesteps)

h = start_hour
t = 0
while t < ntimesteps:
  # 10 min before and after each hour
  ind = (kumpula[:,2]==sim_day) & (((kumpula[:,3] == h) & (kumpula[:,4] < 10)) | \
                                   ((kumpula[:,3] == h-1) & (kumpula[:,4] > 49)))
  pressure_smear[t] = np.nanmean(kumpula[ind,9])
  Ui = kumpula[ind,16]
  WDi = kumpula[ind,23]
  ui = np.nanmean(-Ui * np.sin(WDi * np.pi/180.0))
  vi = np.nanmean(-Ui * np.cos(WDi * np.pi/180.0))
  wd_smear[t] = 270.0 - (np.arctan2(vi, ui) * 180.0 / np.pi)
  h += 1
  t += 1

del Ui, WDi, ui, vi

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

kivenlahti = np.loadtxt(fname_kivenlahti, skiprows=1)
kivenlahti[kivenlahti==-999] = np.nan

uk = np.zeros([ntimesteps,8])
vk = np.copy(uk)
WDk = np.copy(uk)
Uk = np.copy(WDk)
Tk = np.copy(WDk)
TDk = np.copy(WDk)
zk = np.array( [2, 26, 49, 92, 141, 217, 265, 327]) # measurement heights
cols = [7, 11, 15, 19, 23, 27, 31, 35] # data columns including wind speed

h = start_hour
t = 0
while t < ntimesteps:
  ind = (kivenlahti[:,3]==h) & ((kivenlahti[:,4]==0))
  for zi in range(8):
    # Separate u and v from U and WD:
    uk[t,zi] = -kivenlahti[ind,cols[zi]] * np.sin(kivenlahti[ind,cols[zi]+1] * np.pi/180.0)
    vk[t,zi] = -kivenlahti[ind,cols[zi]] * np.cos(kivenlahti[ind,cols[zi]+1] * np.pi/180.0)
    # Wind speed
    Uk[t,zi] = np.nanmean(kivenlahti[ind,cols[zi]])
    # Wind direction
    WDk[t,zi] = 270.0 - (np.arctan2( vk[t,zi], uk[t,zi]) * 180.0 / np.pi)
    # air temperature
    Tk[t,zi] = np.nanmean( kivenlahti[ind,cols[zi]-2] )
    # dew-point temperature
    TDk[t,zi] = np.nanmean( kivenlahti[ind,cols[zi]-1] )

  t += 1
  h += 1

# Check that the dew point temperature does not exceed the temperature
for t in range( ntimesteps ):
  for zi in range( 8 ):
    if 0.98*Tk[t,zi] < TDk[t,zi]:
      print( 'too high dew point temperature: T={:2.1f}C, TD={:2.1f}C --> T={:2.1f}C'\
             ' (0.98*T)'.format( Tk[t,zi], TDk[t,zi], Tk[t,zi] * 0.98 ) )
      TDk[t,zi] = 0.98*Tk[t,zi]

# Change the temperature unit from degrees Celsius to Kelvin
Tk += 273.15
TDk += 273.15

# Derive pressure levels from the barometric equation (i.e. assume neutral stratification with a
# lapse rate of g/cp)
pk = np.zeros(np.shape(Tk))
for t in range(ntimesteps):
  pk[t,:] = barometric_formula(zk, 300, 101325)

# Derive the potential temperature
ptk = np.zeros(np.shape(Tk))
for t in range(ntimesteps):
  ptk[t,:] = Tk[t,:] / exner_function(pk[t,:])

# Calculate the specific humidity from the dew point temperature and pressure
qk = np.copy( uk )
for t in range( ntimesteps ):
  qk[t,:] =  mpcalc.specific_humidity_from_dewpoint(TDk[t,:]*units.kelvin, pk[t,:]*units.pascal)

#%% Read SMEARIII DMPS data

if not precursor:

  dmps = np.genfromtxt(fname_smear_dmps)
  dmid_dmps = dmps[0,2::]
  day_dmps = dmps[:,0]
  ntot_dmps = dmps[:,1] # 1/cm3

  # Change the days from the beginning of the year (day_dmps) to datetime (time_dmps)
  for i in range(len(day_dmps)):

    # days from the beginning of the year to datetime
    timecol = dtime.datetime(sim_year, 1, 1) + dtime.timedelta(day_dmps[i] - 1)

    # separate to year, month, day, hour, minute, second
    timecol = np.array([timecol.year, timecol.month, timecol.day, timecol.hour, timecol.minute,
                        timecol.second]).T
    if i==0:
      time_dmps = timecol
    else:
      time_dmps = np.vstack(( time_dmps, timecol ))

  # Select all measurement times during the simulation time
  ind = (time_dmps[:,1]==sim_month) & (time_dmps[:,2]==sim_day) & \
        (((time_dmps[:,3]>=start_hour)) & ((time_dmps[:,3]<=end_hour)))
  ind[(time_dmps[:,3]==end_hour) & (time_dmps[:,4]>end_min)] = False

  # Calculate the background concentration as the average psd over the simulation time
  dNdlogD_dmps = np.nanmean(dmps[ind,2::], axis=0)

  # Change dNdlogD to N
  dlogdi = np.log10(dmid_dmps[1::]/dmid_dmps[0:-1])[-1]
  dmidi = np.append(dmid_dmps, dmid_dmps[-1]*10**dlogdi)
  psd_dmps = dNdlogD_dmps * np.log10( dmidi[1::]/dmidi[0:-1])

  # Adjust to the PALM-SALSA grid:
  # Read the geometric mean diameter per bin from the existing PIDS_DYNAMIC_N03
  dmid_SALSA = dsout_N03['Dmid'][:].data

  # Get the bin limits
  _, bin_limits = define_bins( [3, 7], [2.5e-9, 15.0e-9, 1e-6] )

  psd_to_palm = np.zeros(len(dmid_SALSA))
  for ib in range(len(dmid_SALSA)):
    ind = (dmid_dmps > bin_limits[ib]) & (dmid_dmps <= bin_limits[ib+1]);
    psd_to_palm[ib] = np.sum( psd_dmps[ind] )

  psd_to_palm *= 1.0e6 # from 1/cm3 to 1/m3

#%% Interpolate Kivenlahti measurements to PALM grid

plt.close('all')

z = dsout_all['z'][:].data
zw = dsout_all['zw'][:].data
u_prof  = np.zeros([ntimesteps, len(z)])
v_prof  = np.copy(u_prof)
w_prof  = np.zeros([ntimesteps, len(zw)])
pt_prof = np.copy(u_prof)
qv_prof = np.copy(u_prof)
p_prof = np.copy(u_prof)
T_prof = np.copy(u_prof)

arrays_out = [u_prof, v_prof, w_prof, pt_prof, qv_prof]

iv = 0
for var in [uk, vk, uk*0, ptk, qk]:
  if iv == 2:
    zi = zw
  else:
    zi = z

  # profiles to interpolate
  zin = np.copy(zk)
  if iv == 4:
    zin = np.append(zin, 600.0)

  # Loop through time steps:
  for t in range(ntimesteps):

    # profiles to interpolate
    varin = np.copy(var[t,:])
    if iv == 4:
      varin = np.append(varin, 0.0)

    ind = ~np.isnan(varin)
    fnearest = interp1d( zin[ind], varin[ind], kind='nearest', fill_value="extrapolate" )
    flinear = interp1d( zin[ind], varin[ind], kind='linear', fill_value="extrapolate" )
    fcubic  = CubicSpline( zin[ind], varin[ind] )

    # Use linear interpolation between the observation points
    ind2 = zi<=np.max(zin[ind])
    arrays_out[iv][t,ind2] = flinear(zi[ind2])

    # Use the nearest value above the highest measurement point
    ind2 = zi>np.max(zin[ind])
    arrays_out[iv][t,ind2] = fnearest(zi[ind2])

  iv += 1

u_prof  = arrays_out[0]
v_prof  = arrays_out[1]
w_prof  = arrays_out[2]
pt_prof = arrays_out[3]
qv_prof = arrays_out[4]

#%% Create arrays

nx = len(dsout_all['x'][:].data)-1
ny = len(dsout_all['y'][:].data)-1
nz = len(dsout_all['z'][:].data)

u_palm  = np.zeros([ntimesteps, nz,   ny+1, nx  ], dtype=float)
v_palm  = np.zeros([ntimesteps, nz,   ny,   nx+1], dtype=float)
w_palm  = np.zeros([ntimesteps, nz-1, ny+1, nx+1], dtype=float)
pt_palm = np.zeros([ntimesteps, nz,   ny+1, nx+1], dtype=float)
qv_palm = np.zeros([ntimesteps, nz,   ny+1, nx+1], dtype=float)
u_palm_kivenlahti_wd = np.copy(u_palm)
v_palm_kivenlahti_wd = np.copy(v_palm)
u_palm_kumpula_wd = np.copy(u_palm)
v_palm_kumpula_wd = np.copy(v_palm)


# ------------------------------------------------------------------------------------------------#

# All meteorological data:

arrays_out = [u_palm, v_palm, pt_palm, qv_palm]
arrays_in  = [u_prof, v_prof, pt_prof, qv_prof]
zs         = [z,      z,      z,       z      ]
nys        = [ny+1,   ny,     ny+1,    ny+1   ]
nxs        = [nx,     nx+1,   nx+1,    nx+1   ]
names      = ['u',    'v',    'pt',    'qv'   ]

print('Both wind speed and direction:')
for a in range(len(arrays_out)):
  print('   {}'.format(names[a]))
  for t in range(ntimesteps):
    for j in range(nys[a]):
      for i in range(nxs[a]):
        arrays_out[a][t,:,j,i] = arrays_in[a][t,:]

u_palm  = arrays_out[0]
v_palm  = arrays_out[1]
pt_palm = arrays_out[2]
qv_palm = arrays_out[3]

# ------------------------------------------------------------------------------------------------#

# Only wind direction (use Kivenlahti data):

arrays_out = [u_palm_kivenlahti_wd, v_palm_kivenlahti_wd]
zs         = [z,      z       ]
nys        = [ny+1,   ny      ]
nxs        = [nx,     nx+1    ]
names      = ['u',    'v'     ]

wdk = 270.0 - (np.arctan2(v_prof, u_prof) * 180.0 / np.pi)

# Copy data from PIDS_DYNAMIC
for a in range(  len( arrays_out ) ):
  arrays_out[a][:,:,:,0]  = dsout_wd['ls_forcing_left_{}'.format(names[a])][:].data
  arrays_out[a][:,:,:,-1] = dsout_wd['ls_forcing_right_{}'.format(names[a])][:].data
  arrays_out[a][:,:,:,-2] = dsout_wd['ls_forcing_right_{}'.format(names[a])][:].data
  arrays_out[a][:,:,-1,:] = dsout_wd['ls_forcing_north_{}'.format(names[a])][:].data
  arrays_out[a][:,:,-2,:] = dsout_wd['ls_forcing_north_{}'.format(names[a])][:].data
  arrays_out[a][:,:,0,:]  = dsout_wd['ls_forcing_south_{}'.format(names[a])][:].data
  arrays_out[a][:,-1,:,:] = dsout_wd['ls_forcing_top_{}'.format(names[a])][:].data

print('Only wind direction (Kivenlahti): both u & v')
for t in range(ntimesteps):
  print('   t={}'.format(t))
  wdh = wdk[t,:]
  for j in range(ny):
    for i in range(nx):
      Uh = np.sqrt(arrays_out[0][t,:,j,i]**2 + arrays_out[1][t,:,j,i]**2)
      uh2 = -Uh * np.sin(wdh * np.pi/180.0)
      vh2 = -Uh * np.cos(wdh * np.pi/180.0)
      wdi = 270.0 - (np.arctan2(vh2, uh2) * 180.0 / np.pi)
      arrays_out[0][t,:,j,i] = uh2
      arrays_out[1][t,:,j,i] = vh2
      arrays_out[0][t,-1,j,i] = uh2[-1]
      arrays_out[1][t,-1,j,i] = vh2[-1]

arrays_out[0][:,:,-1,:] = arrays_out[0][:,:,-2,:]
arrays_out[1][:,:,:,-1] = arrays_out[1][:,:,:,-2]

u_palm_kivenlahti_wd  = arrays_out[0]
v_palm_kivenlahti_wd  = arrays_out[1]

del arrays_out, uh2, vh2

# ------------------------------------------------------------------------------------------------#

# Only wind direction (use Kumpula data):

arrays_out = [u_palm_kumpula_wd, v_palm_kumpula_wd]
zs         = [z,      z       ]
nys        = [ny+1,   ny      ]
nxs        = [nx,     nx+1    ]
names      = ['u',    'v'     ]

# Copy data from PIDS_DYNAMIC
for a in range(len(arrays_out)):
  arrays_out[a][:,:,:,0]  = dsout_wdk['ls_forcing_left_{}'.format(names[a])][:].data
  arrays_out[a][:,:,:,-1] = dsout_wdk['ls_forcing_right_{}'.format(names[a])][:].data
  arrays_out[a][:,:,:,-2] = dsout_wdk['ls_forcing_right_{}'.format(names[a])][:].data
  arrays_out[a][:,:,-1,:] = dsout_wdk['ls_forcing_north_{}'.format(names[a])][:].data
  arrays_out[a][:,:,-2,:] = dsout_wdk['ls_forcing_north_{}'.format(names[a])][:].data
  arrays_out[a][:,:,0,:]  = dsout_wdk['ls_forcing_south_{}'.format(names[a])][:].data
  arrays_out[a][:,-1,:,:] = dsout_wdk['ls_forcing_top_{}'.format(names[a])][:].data

print('Only wind direction (Kumpula): both u & v')
for t in range(ntimesteps):
  print('   t={}'.format(t))
  wdk = u_prof[t,:]*0.0 + wd_smear[t]
  wdh = wdk
  for j in range(ny):
    for i in range(nx):
      Uh = np.sqrt(arrays_out[0][t,:,j,i]**2 + arrays_out[1][t,:,j,i]**2)
      uh2 = -Uh * np.sin(wdh * np.pi/180.0)
      vh2 = -Uh * np.cos(wdh * np.pi/180.0)
      wdi = 270.0 - (np.arctan2(vh2, uh2) * 180.0 / np.pi)
      arrays_out[0][t,:,j,i] = uh2
      arrays_out[1][t,:,j,i] = vh2
      arrays_out[0][t,-1,j,i] = uh2[-1]
      arrays_out[1][t,-1,j,i] = vh2[-1]

arrays_out[0][:,:,-1,:] = arrays_out[0][:,:,-2,:]
arrays_out[1][:,:,:,-1] = arrays_out[1][:,:,:,-2]

u_palm_kumpula_wd  = arrays_out[0]
v_palm_kumpula_wd  = arrays_out[1]

del arrays_out, uh2, vh2

#%% Replace meteorological data in PIDS_DYNAMIC with Kivenlahti measurements

names             = ['u',    'v',     'w',    'pt',    'qv'    ]
arrays            = [ u_palm, v_palm,  w_palm, pt_palm, qv_palm]
arrays_kivenlahti = [ u_palm_kivenlahti_wd, v_palm_kivenlahti_wd ]
arrays_kumpula    = [ u_palm_kumpula_wd, v_palm_kumpula_wd ]
arrays_prof       = [ u_prof, v_prof,  w_prof, pt_prof, qv_prof]

# init_atmosphere_XX (z/zw)
for a in range( len( arrays ) ):
  dsout_all['init_atmosphere_{}'.format(names[a])][:] = arrays_prof[a][0,:]
  if a<2:
    init_kivenlahti = np.nanmean(arrays_kivenlahti[a][0,:,:,0],  axis=-1)
    dsout_wd['init_atmosphere_{}'.format(names[a])][:] = init_kivenlahti
    init_kumpula = np.nanmean(arrays_kumpula[a][0,:,:,0],  axis=-1 )
    dsout_wdk['init_atmosphere_{}'.format(names[a])][:] = init_kumpula

# ls_forcing_left_XX (time,z/zw,y/yv)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_left_{}'.format(names[a])][:] = arrays[a][:,:,:,0]
  if a<2:
    dsout_wd['ls_forcing_left_{}'.format(names[a])][:] = arrays_kivenlahti[a][:,:,:,0]
    dsout_wdk['ls_forcing_left_{}'.format(names[a])][:] = arrays_kumpula[a][:,:,:,0]
  
# ls_forcing_right_XX (time,z/zw,y/yv)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_right_{}'.format(names[a])][:] = arrays[a][:,:,:,-1]
  if a<2:
    dsout_wd['ls_forcing_right_{}'.format(names[a])][:] = arrays_kivenlahti[a][:,:,:,-1]
    dsout_wdk['ls_forcing_right_{}'.format(names[a])][:] = arrays_kumpula[a][:,:,:,-1]

# ls_forcing_north_XX (time,z/zw,x/xu)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_north_{}'.format(names[a])][:] = arrays[a][:,:,-1,:]
  if a<2:
    dsout_wd['ls_forcing_north_{}'.format(names[a])][:] = arrays_kivenlahti[a][:,:,-1,:]
    dsout_wd['ls_forcing_north_{}'.format(names[a])][:] = arrays_kumpula[a][:,:,-1,:]

# ls_forcing_south_XX (time,z/zw,x/xu)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_south_{}'.format(names[a])][:] = arrays[a][:,:,0,:]
  if a<2:
    dsout_wd['ls_forcing_south_{}'.format(names[a])][:] = arrays_kivenlahti[a][:,:,0,:]
    dsout_wdk['ls_forcing_south_{}'.format(names[a])][:] = arrays_kumpula[a][:,:,0,:]

# ls_forcing_top_XX (time,y/yv,x/xu)
for a in range( len( arrays ) ):
  dsout_all['ls_forcing_top_{}'.format(names[a])][:] = arrays[a][:,-1,:,:]
  if a<2:
    dsout_wd['ls_forcing_top_{}'.format(names[a])][:] = arrays_kivenlahti[a][:,-1,:,:]
    dsout_wdk['ls_forcing_top_{}'.format(names[a])][:] = arrays_kumpula[a][:,-1,:,:]

del arrays, arrays_kivenlahti


#%% Replace size distribution data in PIDS_DYNAMIC_N03 with SMEARIII DMPS data

if not precursor:

  # PS: use the chemical composition from ADCHEM

  for k in range( len( dsout_N03['z'][:] ) ):
    dsout_N03['init_atmosphere_aerosol'][k,:] = psd_to_palm

  bounds = ['north','south','left','right']
  for bound in bounds:
    shp = np.shape(dsout_N03['ls_forcing_{}_aerosol'.format(bound)])
    for t in range(ntimesteps):
      for k in range(shp[-3]):
        for i in range(shp[-2]):
          dsout_N03['ls_forcing_{}_aerosol'.format(bound)][t,k,i,:] = psd_to_palm

  bound = 'top'
  shp = np.shape(dsout_N03['ls_forcing_{}_aerosol'.format(bound)])
  for t in range(ntimesteps):
    for j in range(shp[-2]):
      for i in range(shp[-1]):
        dsout_N03['ls_forcing_{}_aerosol'.format(bound)][t,j,i,:] = psd_to_palm


#%% Close files

dsout_all.close()
dsout_wd.close()
dsout_wdk.close()
if not precursor:
  dsout_N03.close()
