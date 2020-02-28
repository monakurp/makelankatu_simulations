import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from psdLib import define_bins, interpolate_psd_from_figure, ntot_from_psd_and_pm

year = 2017
month = 6
day = 9
tod  = 'morning'

biobus = False

nbins = 10  # number of size bins applied in SALSA
nbin = [2, 8]  # number of bins per subrange
reglim = [2.5e-9, 10.0e-9, 1.0e-6]  # subrange diameter limits

folder = '/home/monakurp/makelankatu_simulations'

file_psd_shape = 'source_data/hietikko2018_fig7.csv'
file_emissions = 'source_data/EEA_calculated_EF.csv'
file_emissions_HSY = 'source_data/HSY/emissions_PM_makelankatu_2017.xlsx'

#%% Date and times

date = '{}{:02d}{:02d}'.format(year,month,day)

if ( date=='20170609' and tod=='morning' ):
  start_hour = 7
  start_min  = 0
  end_hour   = 9
  end_min    = 15
  plusUTC    = 3
  
#%% Download data  

os.chdir( ( folder ) )

# Bin diameters and limits
dmid, bin_limits = define_bins( nbin, reglim )
dmid = np.sqrt( bin_limits[0:-1] * bin_limits[1::] )

# Measured number size distribution
psd_hietikko = np.genfromtxt( file_psd_shape, delimiter=',' )

# EEA emission data in g/m/veh: 
emission = pd.read_csv( file_emissions, header=1 )
emission.rename( columns=lambda x: x.strip() )
emission['start_date'] = pd.to_datetime( emission['start time'] )
mask = ( emission['start_date'] >= '{}-{}-{} {}:00:00'.format( date[0:4],date[4:6],date[6::],start_hour ) ) & \
       ( emission['start_date'] <= '{}-{}-{} {}:00:00'.format( date[0:4],date[4:6],date[6::],end_hour+1 ) )
emission = emission.loc[mask]  # select only hours needed

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

#%% Calculate emission

plt.close('all')
do_plot = True

# Define aerosol and heat emissions  
factor = 1.0 / 3600.0
km_in_m = 1000.0

# Emission factor EF in g/m/veh

# EEA:
EF_EEA = np.array( emission[' PM'] )  

# HSY:
TR_HSY = np.sum( emission_HSY.iloc[:,1:5+1].values, axis=1 ) # veh/h
if biobus:
  EF_HSY = emission_HSY.iloc[:,-1].values # g/km/H
else:
  EF_HSY = emission_HSY.iloc[:,-3].values
EF_HSY = EF_HSY / ( km_in_m * TR_HSY )

# Size distribution for Hietikko et al. 
_, rpb_hietikko = interpolate_psd_from_figure( file_psd_shape, bin_limits, do_plot )

for t in range( len( emission ) ): # time

  # Mass emission (PM2.5) to number
  ntot_EEA = ntot_from_psd_and_pm( EF_EEA[t], rpb_hietikko, dmid )
  ntot_HSY = ntot_from_psd_and_pm( EF_HSY[t], rpb_hietikko, dmid )
  
  # Hietikko et al. based on the fuel consumption
  ntot_hietikko_VTT = np.array( emission[' FC_VTT'] )[t] * 4.22E15 * 1e-3
  ntot_hietikko_EEA = np.array( emission[' FC_EEA'] )[t] * 4.22E15 * 1e-3
  
  # Print results:
  print( 'Start time: {}'.format( emission['start time'].iloc[t] ) )
  print( 'Fuel consumption: {}/{}  g/km/veh (VTT, EEA)'.format(
    np.array( emission[' FC_VTT'] )[t]*1e3, np.array( emission[' FC_EEA'] )[t]*1e3 ) )
  
  print( 'Number emission: {:2.3g}, {:2.3g}, {:2.3g}, {:2.3g} #/m/veh'.format(
    ntot_EEA, ntot_HSY, ntot_hietikko_VTT, ntot_hietikko_EEA ) )
  print( '                (EEA, HSY, hietikko VTT, hietikko EEA) \n' ) 
  
  do_plot = False
  
#%% Estimate sensitivity of number emission on the PM emission value
  
#EF = np.array( emission[' PM'] )[0]
EF = 2e-5

# Size distribution for Hietikko et al. 
_, rpb_hietikko = interpolate_psd_from_figure( file_psd_shape, bin_limits, False )

multiplied = np.arange(1, 2, 0.1)
multiplied = np.append( multiplied, np.arange( 2, 10.1, 1 ) )

ntots = ntot_from_psd_and_pm( EF*multiplied, rpb_hietikko, dmid )
  
deriv = 5e-5 * ntots[0] / EF
print('{:2.3g}'.format( deriv ) )
  
fig, ax = plt.subplots()
ax.plot( multiplied, ntots/ntots[0], '*')
