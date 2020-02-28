
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import os
import pandas as pd
from scipy.ndimage import binary_erosion, binary_dilation
from pylab import cm
from psdLib import define_bins, interpolate_psd_from_figure, ntot_from_psd_and_pm
import netCDF4 as nc

plt.close('all')

folder = '/home/monakurp/makelankatu_simulations'
year = 2017
month = 12
day = 7
tod  = 'morning'

# Aerosol size distribution:
nbins = 10  # number of size bins applied in SALSA
nbin = [3, 7]  # number of bins per subrange
reglim = [2.5e-9, 15.0e-9, 1.0e-6]  # subrange diameter limits

fill_value = np.nan#-9999

hourly = True # calculate hourly emissions (if False, calculate every 15 minutes)
biobus = False # assume that all busses use bioful (like HSL buses)
PM_HSY = False # use emission data provided by Anu Kousa from HSY
hietikko = True # use Hietikko et al. to define the number emission instead EEA data
VTT = False # use VTT fuel consumption data instead of EEA

# coordinates of the HSY supersite
#easting_HSY  = 25497337.38
#northing_HSY = 6675959.94
#y_HSY = origo[0] - northing_HSY
#x_HSY = easting_HSY - origo[1]

#%% Exact hours

date = '{}{:02d}{:02d}'.format(year,month,day)

if ( date=='20170609' and tod=='morning' ):
  start_hour = 7
  start_min  = 0
  end_hour   = 9
  end_min    = 15
  plusUTC    = 3

if ( date=='20170609' and tod=='evening' ):
  start_hour = 20
  start_min  = 55
  end_hour   = 21
  end_min    = 15
  plusUTC    = 3

elif ( date=='20170614' and tod=='morning' ):
  start_hour = 6
  start_min  = 55
  end_hour   = 9
  end_min    = 0
  plusUTC    = 3
  
elif ( date=='20171207' and tod=='morning' ):
  start_hour = 7
  start_min  = 5
  end_hour   = 9
  end_min    = 15
  plusUTC    = 2

#%% Different street types around Makelankatu:

st_index = np.arange(1,7.1,1)
st_names = ['asuntokatu','kokoojakatu','paakatu_keskustaan', 'paakatu_keskustasta','kumpulantie',\
            'makelanrinne_keskustaan','makelanrinne_keskustasta']
st_traffic_relative = np.array([200., 5400., 28100., 28100., 8700., 43200., 43200.])
st_traffic_relative = st_traffic_relative / st_traffic_relative[2]
st_traffic_relative[-2::] = 1.0 + 0.3/0.7
st_width = np.array([5., 6., 10., 10., 10., 10., 10.])

#%% Colors:
cmap_init_st = cm.get_cmap('PiYG', 3) 
cmap_new_st  = cm.get_cmap('rainbow', len(st_names) )

#%% Open files:

os.chdir( ( folder ) )

file_child_topo                = 'input_data_to_palm/topo_child.npz'
file_child_lanes_and_buildings = 'input_data_to_palm/lanes_child.npz'
file_child_st                  = 'input_data_to_palm/street_types_child.npz'
file_child_oro                 = 'input_data_to_palm/oro_child.npz'

if not hourly:
  suffix = '_15min'
else:
  suffix = ''

file_emissions = 'source_data/EEA_calculated_EF{}.csv'.format( suffix )
file_psd_shape = 'source_data/hietikko2018_fig7.csv'
file_HDD = 'source_data/HDD_Kaisaniemi.csv'
file_emissions_HSY = 'source_data/HSY/emissions_PM_makelankatu_2017.xlsx'

pids_folder = 'input_data_to_palm/cases/{}_{}'.format( date, tod )
pids_salsa_out = '{}/PIDS_SALSA_N03'.format( pids_folder )
if not hourly:
  pids_salsa_out += '_15min'
if PM_HSY:
  pids_salsa_out += '_HSY'
  if biobus:
    pids_salsa_out += '_biobus' 
else:
  if not hietikko:
    pids_salsa_out += '_EEA'
if hietikko:
  pids_salsa_out += '_hietikko'
  if VTT:
    pids_salsa_out += '_VTT'
  else:
    pids_salsa_out += '_EEA'
  
pids_chem_out = '{}/PIDS_CHEM_N03{}'.format( pids_folder, suffix )

#%% Read in input files

# Maps
child_topo                = np.load( file_child_topo )
child_lanes_and_buildings = np.load( file_child_lanes_and_buildings )
child_st                  = np.load( file_child_st )
child_oro                 = np.load( file_child_oro )

# Emission data in g/m/veh: 
emission = pd.read_csv( file_emissions, header=1 )
emission.rename( columns=lambda x: x.strip() )
emission['start_date'] = pd.to_datetime( emission['start time'] )
mask = ( emission['start_date'] >= '{}-{}-{} {}:00:00'.format( date[0:4],date[4:6],date[6::],start_hour ) ) & \
       ( emission['start_date'] <= '{}-{}-{} {}:00:00'.format( date[0:4],date[4:6],date[6::],end_hour+1 ) )
emission = emission.loc[mask]  # select only hours needed

# Heating degree days
HDD_data = np.genfromtxt( file_HDD, delimiter=',' ) # year, month, day, T, HDD

# HSY emission data in g/km (multiplied by the traffic rate already)
def dateparse_fleetdistr(date_string):
    dt = pd.datetime.strptime(date_string, '%Y %m %d %H')
    return dt
emission_HSY = pd.read_excel( file_emissions_HSY, 'Taul1', delimiter=',', header=4, 
                              date_parser=dateparse_fleetdistr, 
                              parse_dates={'datetime' : ['Year','Month','Day','Hour']})
mask = ( emission_HSY['datetime'] >= '{}-{}-{} {}:00:00'.format( date[0:4],date[4:6],date[6::],start_hour ) ) & \
       ( emission_HSY['datetime'] <= '{}-{}-{} {}:00:00'.format( date[0:4],date[4:6],date[6::],end_hour+1 ) )
emission_HSY = emission_HSY.loc[mask]  # select only hours needed

#%% Manipulate lane and street maps

# Erode the street type data to make street areas slightly narrower (needed for filtering lane data)
st = child_st['R']
st = binary_erosion( binary_erosion( st ) )
street_types = np.copy(child_st['R'])
street_types[st==False] = np.nan

# Filter the lane map to create a street type map
street_type_map = np.copy( child_lanes_and_buildings['R'] )
street_type_map[ street_type_map==6 ] = fill_value # exclude buildings
street_type_map[ street_type_map==2 ] = fill_value # set surface to nan
street_type_map[ ( ( street_type_map==4 ) & ( street_types==4 ) ) ] = 4 # main streets
street_type_map[ ( ( street_type_map==3 ) & ( street_types==4 ) ) ] = 4 # main streets
street_type_map[ ( ( street_type_map==3 ) & ( street_types==5 ) ) ] = 5 # main streets

# Increase the width of main street lanes
main_streets = binary_dilation( binary_dilation( street_type_map==4 ) )
street_type_map[ main_streets ] = 4

# Normalise street type values to start from 1
street_type_map = street_type_map - np.nanmin( street_type_map ) + 1

#%% Further separate lanes based on traffic rates

plt.close('all')

from PIL import Image, ImageDraw

# polygon points for separating different traffic rate areas
x1 = [  0,   0,   0,   0,   0,   0]
x2 = [554, 321, 214,   0, 495,  64]
x3 = [  0,   0,   0, 245,   0,   0]
y1 = [123, 343, 248,   0, 158,  95]
y2 = [576,   0,   0, 290, 576, 176]
y3 = [576,   0,   0,   0, 576, 248]

nx = np.shape( street_type_map )[0]
ny = np.shape( street_type_map )[1]

# emission map
st_map = np.zeros( np.shape( street_type_map ) ) + np.nan

# mask for separating different traffic rate areas
mask = np.zeros( np.shape( street_type_map ) )

# Create the mask
for i in range( len( x1 ) ):
  
  polygon = [(x1[i],y1[i]),(x2[i],y2[i]),(x3[i],y3[i])]
  
  maskIm = Image.new( 'L', ( st_map.shape[1], st_map.shape[0]), 0 )
  ImageDraw.Draw( maskIm ).polygon( polygon, outline=1, fill=1 )
  mask += np.array( maskIm ) * (i+1)
 
# Separate different areas / lane types:

# asuntokatu
st_map[ street_type_map==1 ] = st_index[0]

# kokoojakatu
st_map[ street_type_map==3 ] = st_index[1]

# Makelankatu to the centre:
st_map[ ( ( mask==1 ) & ( street_type_map==2 ) ) ] = st_index[2]

# Makelankatu from the city centre
st_map[ ( ( mask==0 ) & ( street_type_map==2 ) ) ] = st_index[3]

# Kumpulantie
st_map[ ( ( mask==12 ) & ( street_type_map==2 ) ) ] = st_index[4]

# Makelankadun pohjoisosa:
# Etelaan:
st_map[ ( ( ( mask==3 ) | ( mask==7 ) | ( mask==8 ) | (mask==13 ) | (mask==15 ) | (mask==16 )) & \
          ( street_type_map==2 ) ) ] = st_index[5]
# Pohjoiseen:
st_map[ ( ( ( mask==2 ) | ( mask==6 ) | (mask==9 ) ) & ( street_type_map==2 ) ) ] = st_index[6]


# Plot maps:
gs   = gridspec.GridSpec( 1, 2, width_ratios=[20,1], left=0.05, right=0.69, wspace=0.02  )

# Plot relative traffic amount map
fig2   = plt.figure()
ax2    = fig2.add_subplot(gs[0])
im21   = ax2.imshow( child_topo['R'], cmap='Greys', interpolation='None' )
im22   = ax2.imshow( st_map, interpolation='None', cmap=cmap_new_st )
cax2   = fig2.add_subplot(gs[1])
cbar2  = fig2.colorbar( im22, cax=cax2, cmap=cmap_new_st, 
                        ticks=np.linspace(1.5,6.5,len(st_names) ) )
cbar2.ax.set_yticklabels(st_names, fontsize=10 )

#%% Heat emissions from buildings and traffic

# Anthropogenic heat from buildings
a0 = 6.8e6 * 3600.0 # W*h/capita * J/W*h = J/capita
a2 = 0.0149 * 1.0e4 # W/K/capita

HDD = HDD_data[ ( ( HDD_data[:,0]==year ) & ( HDD_data[:,1]==month ) & \
                  ( HDD_data[:,2]==day ) ), -1 ]

QF_per_capita = ( a0 + a2*HDD ) 

#%% Calculate pollutant emissions

emission_name    = ['NO                       ','NO2                      ',
                    'OCSV                     ','RH                       ',
                    'H2SO4                    ','NH3                      ']  
composition_name = ['BC                       ','OC                       ', 
                    'H2SO4                    '] 

nspecies = len( emission_name )
ncc = len( composition_name )
ncat = np.arange( 1, 1.1, 1 )

emission_index = np.linspace( 1, nspecies, nspecies )
composition_index = np.linspace( 1, ncc, ncc )

dmid, bin_limits = define_bins( nbin, reglim )

# Traffic rates
tr = np.zeros([len(emission),len(st_traffic_relative)])
# asuntokaduilla, kokoojakaduilla ja Kumpulantiella Makelankadun summaan verrannollinen
for i in [0,1,4]:
  tr[:,i] = 0.5* ( emission[' traffic to city (veh/h)'] + \
            emission[' traffic from city (veh/h)'] ) * st_traffic_relative[i]
# Makelankadulla mitattu
tr[:,2] = emission[' traffic to city (veh/h)']
tr[:,3] = emission[' traffic from city (veh/h)']
# Makelankadun pohjoisosassa 30% enemman
tr[:,5] = emission[' traffic to city (veh/h)'] * st_traffic_relative[5]
tr[:,6] = emission[' traffic from city (veh/h)'] * st_traffic_relative[6]


# EF [g/m/veh] * TR[veh/h] * 1/3600 h/s * 1 / street_width --> g/m2/s

emission_values = np.zeros( [len( emission ), 1, ny, nx, nspecies] ) + 0.0
aerosol_emission_values = np.zeros( [len( emission ), ny, nx, len( ncat ) ] ) + 0.0
heat_traffic = np.zeros( [len( emission ), ny, nx] ) + 0.0

# Define gas emissions
si = 0
for s in emission_name:
  sname = emission_name[si].strip()
  if hourly:
    factor = 1.0 / 3600.0
  else:
    factor = 1.0 / 900.0
  if sname == 'RH':
    sname = 'alkanes'
  elif sname == 'H2SO4':
    sname = 'SO2'
    factor = factor * 0.1 # EF_[H2SO4] = 0.1 * EF_[SO2]

  for t in range( len( emission ) ): # time
    EF = np.array( emission[' {}'.format(sname)] )[t] # g/m/veh

    for i in range( len( st_traffic_relative ) ): # street type
      TR = tr[t,i]  # veh/s
      is_street = ( st_map==st_index[i] )
      emission_values[t,0,is_street,si] = EF * TR * 1./st_width[i] * factor

  si += 1

# Define aerosol and heat emissions  
if hourly:
  factor = 1.0 / 3600.0
else:
  factor = 1.0 / 900.0

# emission factor EF in g/m/veh
if PM_HSY:
  TR_HSY = np.sum( emission_HSY.iloc[:,1:5+1].values, axis=1 ) # veh/h
  if biobus:
    EF = emission_HSY.iloc[:,-1].values # g/km/H
  else:
    EF = emission_HSY.iloc[:,-3].values
  EF = EF / 1000.0 / TR_HSY
else:
  EF = np.array( emission[' PM'] )  

do_plot = True

# Size distribution for Hietikko et al. 
_, rpb_hietikko = interpolate_psd_from_figure( file_psd_shape, bin_limits, do_plot )

for t in range( len( emission ) ): # time
  
  for i in range( len( st_traffic_relative ) ): # street type
    TR = tr[t,i]  # veh/s
    is_street = ( st_map==st_index[i] )
    
    if hietikko:
      # Hietikko et al. based on the fuel consumption
      if VTT:
        ntot = np.array( emission[' FC_VTT'] )[t] * 4.22E15 * 1e-3
      else:
        ntot = np.array( emission[' FC_EEA'] )[t] * 4.22E15 * 1e-3
    else:
      # Mass emission (PM2.5) to number
      ntot = ntot_from_psd_and_pm( EF[t], rpb_hietikko, dmid )
    
    aerosol_emission_values[t,is_street,:] = ntot * TR / st_width[i] * factor
    heat_traffic[t,is_street] = np.array( emission[' E (J/m)'] )[t] * TR / st_width[i] \
                                * factor # J/(m*veh) * veh/s * 1/m = J/m2s = W/m2
    do_plot = False                             

number_fracs = np.zeros([ len( ncat ), nbins ])
number_fracs[0,:] += rpb_hietikko

emission_mass_fracs = np.zeros( [len( ncat ), ncc], dtype = float )
emission_mass_fracs[0,0] = np.mean( emission[' BC'] / emission[' PM'] )
emission_mass_fracs[0,1] = np.mean( emission[' OC'] / emission[' PM'] )
emission_mass_fracs[0,2] = 1.0 - np.sum( emission_mass_fracs[0:-1] )

#dmid = np.sqrt( bin_limits[0:-1]*bin_limits[1::] )

emission_values[emission_values==0] = -9999.0
aerosol_emission_values[aerosol_emission_values==0] = -9999.0
  
#%% Save into a file: aerosol and gas emissions

pids_static = nc.Dataset( '{}/PIDS_STATIC_N03'.format( pids_folder, tod ), 'r')
pids_salsa  = nc.Dataset( pids_salsa_out, 'w', format='NETCDF4' )
pids_chem   = nc.Dataset( pids_chem_out, 'w', format='NETCDF4' )

dims = ['x','y']
max_string_length = np.linspace( 1, 25, 25 )

seconds_in_hour = 3600.0
time_start = ( start_hour - plusUTC ) * seconds_in_hour
if hourly:
  dt = 3600.0
else:
  dt = 900.0
time_emission = np.arange( time_start, time_start + ( len( emission ) - 1 )*dt+1, dt ) 
time_emission[0] += start_min*60
time_emission -= time_emission[0]

for dsout in [pids_salsa, pids_chem]:
  
  # Global attributes:
  dsout.setncatts( pids_static.__dict__ )
  dsout.author = 'Mona Kurppa (mona.kurppa@helsinki.fi)'
  dsout.contact_person = 'Mona Kurppa (mona.kurppa@helsinki.fi)'
  dsout.source = 'http://lipasto.vtt.fi, HSY (personal communication),'\
    'https://www.eea.europa.eu/publications/emep-eea-guidebook-2016/part-b-sectoral-guidance-chapters/1-energy/1-a-combustion/1-a-3-b-i/view'
 
  # Dimensions:
  for name, dimension in pids_static.dimensions.items():
    if name in dims:
      dsout.createDimension( name, (len( dimension ) if not dimension.isunlimited() else None))
  for name, variable in pids_static.variables.items():
    if name in dims:
      x = dsout.createVariable( name, variable.datatype, variable.dimensions )
      dsout[name][:] = pids_static[name][:]
      dsout[name].setncatts( pids_static[name].__dict__ ) # attributes
  
  # Time: 
  dsout.createDimension( 'time', len( time_emission ) )
  time_emissionn = dsout.createVariable( 'time', 'f4', ('time',) )
  time_emissionn[:] = time_emission
  time_emissionn.long_name = 'emission data time step (seconds in time utc)'
  time_emissionn.standard_name = 'emission_time_step'  
  
  dsout.createDimension( 'max_string_length', len( max_string_length ) )
  cai = dsout.createVariable( 'max_string_length', 'i4', ('max_string_length',) )
  cai[:] = max_string_length


# Only for chemistry:

# z:
pids_chem.createDimension( 'z', 1 )
zi = pids_chem.createVariable( 'z', 'f4', ('z',),)
zi[:] = 0.0
zi.units = 'm'
                               
# number of emission species:   
pids_chem.createDimension( 'nspecies', len( emission_index ) )
nspeciesi = pids_chem.createVariable( 'nspecies', 'i4', ('nspecies',),)
nspeciesi[:] = emission_index                         

# indices of emitted species:                                  
ei = pids_chem.createVariable( 'emission_index', 'u2', ('nspecies',) )
ei[:] = emission_index 
ei.long_name = 'emission species index'
ei.standard_name = 'emission_index'  

# namelist of all emitted species:
en = pids_chem.createVariable( 'emission_name', 'S1', ('nspecies','max_string_length',) )
en[:] = list( map(lambda x : list(x), emission_name ) )
en.long_name = 'emission species name'
en.standard_name = 'emission_name'

# gas emissions
ev = pids_chem.createVariable( 'emission_values', 'f4', ('time','z','y','x','nspecies',), fill_value=-9999.0 )
evmod = np.copy( emission_values[:,:,::-1,:,:] )
ev[:] = evmod
ev.units = 'g/m2/s'                      
ev.long_name = 'emission values'
ev.standard_name = 'emission_values'
ev.lod = 2


# Only for salsa:    

# category number
pids_salsa.createDimension( 'ncat', len( ncat ) )
ncati = pids_salsa.createVariable( 'ncat', 'i4', ('ncat',) )
ncati[:] = ncat
 
# category name
emission_category_name = ['traffic exhaust          ']
ecn = pids_salsa.createVariable( 'emission_category_name', 'S1', ('ncat','max_string_length',) )
ecn[:] = list( map( lambda x : list(x), emission_category_name ) )
ecn.long_name = 'emission category name'
ecn.standard_name = 'emission_category_name'

# indices of chemical components:          
pids_salsa.createDimension( 'composition_index', len( composition_index ) )                  
cai = pids_salsa.createVariable('composition_index', 'i4', ('composition_index',) ) 
cai[:] = composition_index  
                          
# namelist of chemical components:                            
cn = pids_salsa.createVariable( 'composition_name', 'S1', ('composition_index','max_string_length',) )
cn[:] = list( map(lambda x : list(x), composition_name ) )
cn.long_name = 'aerosol composition name'
cn.standard_name = 'composition_name'  
         
# emission mass fractions                 
ca = pids_salsa.createVariable( 'emission_mass_fracs', 'f4', ('ncat','composition_index',),
                                fill_value=-9999.0 )  
ca.units = ''
ca[:] = emission_mass_fracs
ca.long_name = 'mass fractions of chemical components in aerosol emissions'
ca.standard_name = 'emission_mass_fractions'    

# mean diameters of aerosol size bins:  
pids_salsa.createDimension( 'Dmid', len( dmid ) )                      
dmidi = pids_salsa.createVariable( 'Dmid', 'f4', ('Dmid',))  
dmidi[:] = dmid
dmidi.units = 'm'

# Number fractions per bin in aerosol emissions:                    
enf = pids_salsa.createVariable( 'emission_number_fracs', 'f4', ('ncat','Dmid',), fill_value=-9999.0 )
enf.units = ''
enf[:] = number_fracs
enf.long_name = 'number fractions of aerosol size bins in aerosol emissions'
enf.standard_name = 'emission_number_fractions' 

# aerosol emissions
aev = pids_salsa.createVariable( 'aerosol_emission_values', 'f4', ('time','y','x','ncat',), 
                                 fill_value=-9999.0 )
aevmod = np.copy( aerosol_emission_values[:,::-1,:,:] )
aev[:] = aevmod
aev.units = '#/m2/s'
aev.long_name = 'aerosol emission values'
aev.standard_name = 'aerosol_emission_values'
aev.lod = 2

pids_salsa.close()
pids_chem.close()
pids_static.close()

print('Saved {} successfully!'. format( pids_salsa_out ) )
print('Saved {} successfully!'. format( pids_chem_out ) )

#%% Save into a file: anthropogenic heat
  
topo_top = np.ceil( child_oro['R'] ).astype(int)
  
# ANTHROPOGENIC_HEAT and ANTHROPOGENIC_HEAT_PROFILE (which is set to 1)
f1 = open('input_data_to_palm/cases/{}_{}/ANTHROPOGENIC_HEAT'.format( date, tod ), 'w+')
f2 = open('input_data_to_palm/cases/{}_{}/ANTHROPOGENIC_HEAT_PROFILE'.format( date, tod ), 'w+')
  
for i in range( len( street_types ) ):
  for j in range( len( street_types ) ):
    if heat_traffic[0,j,i]>0.0:
      f1.write( '{},{},{},{}\n'.format(i,ny-1-j,int(topo_top[j,i]),np.nanmean( heat_traffic[:,j,i] ) ) )
for t in range(24):
  f2.write( '{},{},{}\n'.format(t,int(topo_top[j,i]),1.0))
f1.close()
f2.close()

print('naheatlayers: {}'.format( np.max( topo_top[heat_traffic[0,:,:]>0] ) +1 ) )
  
#%% Create a test set for anthropogenic heat

if month==6 and day==9 and tod == 'morning':

  nx_test = 20
  ny_test = 20
  
  topo_top_test = np.zeros( [ny_test, nx_test] ) + 4.0
  heat_traffic_test = np.zeros( [ny_test, nx_test] )
  
  heat_traffic_test[1:3,:]   = np.nanpercentile( heat_traffic[heat_traffic>0], 10 )
  heat_traffic_test[5:7,:]   = np.nanpercentile( heat_traffic[heat_traffic>0], 30 )
  heat_traffic_test[9:11,:]  = np.nanpercentile( heat_traffic[heat_traffic>0], 50 )
  heat_traffic_test[13:15,:] = np.nanpercentile( heat_traffic[heat_traffic>0], 70 )
  heat_traffic_test[17:19,:] = np.nanpercentile( heat_traffic[heat_traffic>0], 90 )
  
  f1 = open('input_data_to_palm/cases/{}_{}/ANTHROPOGENIC_HEAT_test'.format( date, tod ), 'w+')
  f2 = open('input_data_to_palm/cases/{}_{}/ANTHROPOGENIC_HEAT_PROFILE_test'.format( date, tod ), 'w+')
      
  for i in range( nx_test ):
    for j in range( ny_test ):
      f1.write( '{},{},{},{}\n'.format(i,j,int(topo_top_test[j,i]+1),heat_traffic_test[j,i]) )
  for t in range(24):
    f2.write( '{},{},{}\n'.format(t,int(topo_top_test[j,i]+1),1.0))
  f1.close()
  f2.close()