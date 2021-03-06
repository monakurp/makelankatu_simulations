# -*- coding: utf-8 -*-
import numpy as np
import os
import pandas as pd

'''
Description:
Calculate the unit traffic emissions for Makelankatu in Helsinki

Author: Mona Kurppa
        mona.kurppa@helsinki.fi
        Univeristy of Helsinki
'''

# File path:
folder =  '/home/monakurp/makelankatu_simulations/source_data'
os.chdir((folder))

# Calculate hourly emission values (if False, calculate every 15 minutes)
hourly = True

#%% Inline functions

def dateparse_traffic(d,t):
  '''For parsing traffic data'''
  dt = d + " " + t
  return pd.datetime.strptime(dt, '%d.%m.%Y %H:%M:%S')

def dateparse_fleetdistr(date_string):
  '''For parsing fleet distribution data'''
  dt = pd.datetime.strptime(date_string, '%Y %m %d %H')
  return dt


#%% Read in traffic data

# Filenames
filename = ['traffic/TrafficDataSummer.csv',
            'traffic/TrafficDataWinter.csv']

# Read in data
for i in range(len(filename)):
  data = pd.read_csv(filename[i], delimiter=',', header=0, date_parser=dateparse_traffic,
                      parse_dates={'datetime' : ['pvm', 'aika']})
  if i==0:
    traffic = data
  else:
    traffic = pd.concat([traffic,data])

del data
traffic = traffic.set_index(traffic['datetime'])
traffic = traffic.drop(['datetime'], axis = 1)

#%% Read in vehicle fleet distribution

# Filename
filename = 'HSY/emissions_PM_makelankatu_2017.xlsx'

fleet_distr = pd.read_excel(filename, 'Taul1', delimiter=',', header=4,
                            date_parser=dateparse_fleetdistr,
                            parse_dates={'datetime' : ['Year','Month','Day','Hour']})

# Driving speed
speed = fleet_distr[['Nopeus']]
speed = speed.set_index(fleet_distr['datetime'])

# Fleet distribution
fleet_distr = fleet_distr.iloc[:,0:6] # select only fleet distribution
fleet_distr[['HA','PA','KAIP','KAP','LA']] = \
  fleet_distr[['Ha','pa','ka','ra','la']].div( fleet_distr.iloc[:,1:6].sum( axis=1 ), axis=0 )
fleet_distr = fleet_distr.set_index(fleet_distr['datetime'])
fleet_distr = fleet_distr.drop(['datetime','Ha','pa','ka','ra','la'], axis=1)

# HA = henkiloauto
# PA = pakettiauto
# LA = linja-auto
# KAIP = ka = kuorma-auto
# KAP  = ra = rekka-auto

#%% Vehicle performance distributions

# Filename
filename = 'VTT_suoritejakaumat.xlsx'

# Read in data
suorite = pd.read_excel(filename, skiprows=4, skipfooter=4, usecols=np.arange(0,9,1))

# drop uncecessary data rows
suorite = suorite.drop(suorite.index[[9,10,20,21,26,27,33,34]])
suorite = suorite.set_index(suorite[2017])
suorite = suorite.drop([2017], axis=1)
suorite.index.name = 'luokka'

# Use HSL information for buses
suorite.loc['LA kaasu']=0.0
suorite['Euro 6'].loc['LA sahko BEV'] = 0.004
suorite['Euro 0'].loc['LA diesel'] = 0.0
suorite['Euro 1'].loc['LA diesel'] = 0.0
suorite['Euro 2'].loc['LA diesel'] = 0.001
suorite['Euro 3'].loc['LA diesel'] = 0.033
suorite['Euro 4'].loc['LA diesel'] = 0.004
suorite['Euro 5'].loc['LA diesel'] = 0.415 # includes also EEV
suorite['Euro 6'].loc['LA diesel'] = 0.543

#%% Fuel consumption from VTT

# Filename
filename = 'VTT_fuel_consumption.xlsx'

# Read in data
FC_VTT = pd.read_excel(filename, skiprows=1, skipfooter=4, index_col=[0])
FC_VTT.index.name = 'luokka'
FC_VTT = FC_VTT.fillna(0)

#%% Read in European Environmental Agency (EEA) emission data

# Filename
filename = 'EEA_road_transport_hot_emission_factor_2017_annex.xlsx'

# Read in data
EEA_data = pd.read_excel(filename, 'HOT_EMISSIONS_PARAMETERS final')

# Replace NaN by 0
EEA_data = EEA_data.fillna(0)

## Unify column names:
EEA_data['Euro Standard'].loc[EEA_data['Euro Standard']=='Euro 1'] = 'Euro I'
EEA_data['Euro Standard'].loc[EEA_data['Euro Standard']=='Euro 2'] = 'Euro II'
EEA_data['Euro Standard'].loc[EEA_data['Euro Standard']=='Euro 3'] = 'Euro III'
EEA_data['Euro Standard'].loc[EEA_data['Euro Standard']=='Euro 4'] = 'Euro IV'
EEA_data['Euro Standard'].loc[EEA_data['Euro Standard']=='Euro 5'] = 'Euro V'
EEA_data['Euro Standard'].loc[EEA_data['Euro Standard']=='Euro 6'] = 'Euro VI'
EEA_data['Euro Standard'].loc[EEA_data['Euro Standard']=='Euro 6 up to 2016'] = 'Euro VI'
EEA_data['Euro Standard'].loc[EEA_data['Euro Standard']=='Euro 6 up to 2017'] = 'Euro VI'

# Biobuses

bb = EEA_data[((EEA_data['Category']=='Buses') & (EEA_data['Fuel']=='Biodiesel') &
               (EEA_data['Road Slope']==0) & (EEA_data['Load']==1))]

# - Use “passenger car medium” in EEA for “HA” in VTT
# - Use “ECE-15.04” or “Conventional” in EEA “Euro 0” in VTT
# - Use the same emissions of PM, NOx and VOC for FFV as for gasoline
# - Use Euro 6 up to 2016
# - Use “CNG Bifuel – CNG” in EEA for “kaasu” in VTT
# - Use “Light commercial vehicle N1-II” in EEA for “PA” in VTT
# - Use “Large-SUV-Executive” in EEA for “PA kaasu” in VTT
# - Use “Urban Buses Standard 15-18 t” with road slope=0 and load=1 in EEA for “LA” in VTT
# - For “KAIP” in VTT, use an average of “Rigid <= 7,5 t” and “Rigid 14-20 t” with a road slope=0
#     and load=0.5
# - For “KAP” in VTT, use an average of “Articulated 34-40 t” and “Articulated 50-60 t” with a road
#     slope=0 and load=0.5
# - For “KAIP” in VTT, use "pieni jakelukuorma-auto" with full load.
# - For “KAP” in VTT, use "puoliperavaunuyhdistelma" with full load.

# Details on NH3 and N2O emissions
partition = pd.read_excel('EEA_data_collected_details.xls')

# Fuel consumption
FC = pd.read_excel('EEA_FC.xls')

#%% Calculate emission factors

year = 2017

slope    = 0
mode     = 'Urban Peak'
pol      = 'PM Exhaust'#'NOx'#'PM Exhaust'
x_biobus = 0.3  # share of biobuses of all diesel buses

# Vehicle categories:
# EEA:
EEA_cat = ['Passenger Cars', 'Light Commercial Vehicles', \
           'Heavy Duty Trucks', 'Heavy Duty Trucks', 'Buses']
# VTT:
VTT_cat = ['HA', 'PA', 'KAIP', 'KAP', 'LA']

# Fuel categories:
# EEA:
EEA_fuel = ['Petrol',   'Petrol', 'Diesel', 'CNG Bifuel ~ Petrol',
            'Petrol',              'Diesel',            'None',           'None' ]
# VTT
VTT_fuel = ['bensiini', 'FFV',    'diesel', 'kaasu',
            'sahko PHEV bensiini', 'sahko PHEV diesel', 'sahko BEV', 'vety' ]

# Euro categories:
# EEA:
EEA_euros = ['Euro 0', 'Euro I', 'Euro II', 'Euro III','Euro IV', 'Euro V', 'Euro VI']
#VTT
VTT_euros = ['Euro 0', 'Euro 1', 'Euro 2', 'Euro 3','Euro 4', 'Euro 5', 'Euro 6'] 

# Pollutants included:
pols = ['PM Exhaust', 'NOx', 'VOC', 'SO2', 'N2O',  'NH3', 'E_heat', 'FC']


# Output filename:
if hourly:
  fileout = 'EEA_calculated_EF.csv'
else:
  fileout = 'EEA_calculated_EF_15min.csv'

# Open a new output file:
f1 = open(fileout, 'w+')

# Write the header:
header = 'start time, end time, traffic to city (veh/h), traffic from city (veh/h), PM, BC, OC, '\
         'NOx, NO, NO2, VOC, OCSV, alkanes, SO2, N2O, NH3, E (J/m), FC_VTT, FC_EEA'
print("# EF in g/m/veh assuming a fleet distribution similar to Makelankatu", file=f1)
print(header, file=f1)


# Start the loops:
ti = 1

# Months:
for month in [6, 12]:

  if month==6:
    days = [9, 14]
  elif month==12:
    days = [5, 7]

  # Days:
  for day in days:

    # Hours:
    for hour in range(24):

      if hourly:
        minutes = [0]
      else:
        minutes = [0, 15, 30, 45]

      # Time strings:
      timestr  = '{}-{:02d}-{:02d} {:02d}:00:00'.format( year, month, day, hour )
      timestr2 = '{}-{:02d}-{:02d} {:02d}:59:00'.format( year, month, day, hour )

      # Driving speed:
      V = speed['Nopeus'].loc[timestr]

      # Traffic rate: 30% less traffic at the supersite than at the traffic measurement location
      T =  0.7 * ( traffic['autot'].loc[timestr:timestr2].sum() + \
                   traffic['autot.1'].loc[timestr:timestr2].sum() )
      to_city   = traffic['autot'].loc[timestr:timestr2].sum() * 0.7
      from_city = traffic['autot.1'].loc[timestr:timestr2].sum() * 0.7

      # Initialise the data row:
      datarow = '{}, {}, {}, {}'.format(timestr, timestr2, to_city, from_city)

      # Print loop information on the console:
      print( '{}-{} ({}/{})'.format(timestr, timestr2, ti, 4*24))

      # Minutes (if hourly = False)
      for minute in minutes:
        timestr_min  = '{}-{:02d}-{:02d} {:02d}:{:02d}:00'.format(year, month, day, hour, minute)
        timestr2_min = '{}-{:02d}-{:02d} {:02d}:{:02d}:00'.format(year, month, day, hour, minute+14)

        if not hourly:  # use 15 min traffic data if hourly==False
          # Traffic rate: 30% less traffic at the supersite than at the traffic measurement location
          T =  0.7 * (traffic['autot'].loc[timestr_min] + traffic['autot.1'].loc[timestr_min])

          to_city   = 0.7 * traffic['autot'].loc[timestr_min]
          from_city = 0.7 * traffic['autot.1'].loc[timestr_min]

          datarow = '{}, {}, {}, {}'.format( timestr_min, timestr2_min, to_city, from_city )

        if minute == 0:

          # Initialise EF data row:
          EF_row = ''

          # Loop through different pollutants:
          for pol in pols:

            EF  = 0
            BC  = 0
            OC  = 0
            NO  = 0
            NO2 = 0
            FCi = 0

            ci = 0 # category index

            # Loop through vehicle categories: hereafter seg
            for cat in VTT_cat:

              if cat=='HA':
                seg = ['Medium']
                load = 0.0
              elif cat=='PA':
                seg = ['N1-II']
                load = 0.0
              elif cat=='KAIP':
                seg = ['Rigid <=7,5 t', 'Rigid 14 - 20 t']
                load = 0.5
              elif cat=='KAP':
                seg = ['Articulated 34 - 40 t', 'Articulated 50 - 60 t']
                load = 0.5 
              elif cat=='LA':
                seg = ['Urban Buses Standard 15 - 18 t']
                load = 1.0 

              fi = 0 # fuel index

              # Loop through fuel types:
              for fuel in VTT_fuel:

                # Loop through EURO classes:
                for ei in range(len(EEA_euros)):

                  EEA_euro = EEA_euros[ei]
                  cati = EEA_cat[ci]
                  fueli = EEA_fuel[fi]
                  E = 0

                  # Redefine some special cases:
                  if EEA_euro=='Euro 0':
                    if cat=='HA' and EEA_fuel[fi]=='Petrol':
                      EEA_euro = 'ECE 15/04' # petrol
                    else:
                      EEA_euro = 'Conventional' # diesel

                  if ((fuel=='kaasu')  & (cat=='PA')): # Light Commercial Vehicle running with gas
                    cati = 'Passenger Cars'
                    seg = ['Large-SUV-Executive']

                  if ((fuel=='kaasu')  & (cat=='LA')): # bus running with gas
                    fueli = 'CNG'
                    seg = ['Urban CNG Buses']
                    slope = 0.0
                    load = 0.0
                    if ei==5:
                      EEA_euro = 'EEV'

                  # Define the performance fraction ("suoriteosuus"):
                  if '{} {}'.format( cat, fuel ) in suorite.index:
                    s = suorite.loc['{} {}'.format(cat, fuel)][VTT_euros[ei]]
                  else:
                    continue

                  # For PM, NOx and VOC, find matching rows:
                  if pol in str(['PM Exhaust', 'NOx', 'VOC']):
                    r = EEA_data[((EEA_data['Category']                    == cati    ) &
                                  (EEA_data['Fuel']                        == fueli   ) &
                                  (EEA_data['Euro Standard']               == EEA_euro) &
                                  (EEA_data['Pollutant/Energy consumption']== pol     ) &
                                  (EEA_data['Road Slope']                  == slope   ) &
                                  (EEA_data['Load']                        == load    ))]

                    # Find matching segments:
                    if len(seg)==1:
                      r = r[r['Segment']==seg[0]]
                    else:
                      r = r[((r['Segment']==seg[0]) | (r['Segment']==seg[1]))]

                    # Find matching driving mode:
                    rr = r[r['Mode']==0]
                    if len(rr)==0:
                      rr =  r[r['Mode']=='Urban Peak']

                    # If data exists, calculate the emission factor
                    if len(rr)>0:
                      E = ( rr['Alpha'] * V**2 + rr['Beta'] * V + rr['Gamma'] + rr['Delta'] / V)\
                          / ( rr['Epsilon'] * V**2 + rr['Zita'] * V + rr['Hta']) * \
                          ( 1 - rr['Reduction Factor [%]'] )
                      E = np.nanmean(E.values)

                    # For hybrid vehicles, emission factor ~80% of a regular vehicle (VTT)
                    if (fuel=='sähkö PHEV bensiini' or fuel=='sähkö PHEV bensiini'):
                      E = 0.8 * E

                    # Include biobuses
                    if ( cat=='LA' and fuel==VTT_fuel[-1]):

                      fueli = 'Biodiesel'
                      bbb = bb[(bb['Pollutant/Energy consumption']== pol     ) &
                               (bb['Euro Standard']               == EEA_euro)]
                      Eb = (bbb['Alpha'] * V**2 + bbb['Beta'] * V + bbb['Gamma'] + bbb['Delta'] / V)\
                           / (bbb['Epsilon'] * V**2 + bbb['Zita'] * V + bbb['Hta']) * \
                           (1 - bbb['Reduction Factor [%]'])
                      E += np.nanmean(Eb.values)

                      s = suorite.loc['LA diesel'][VTT_euros[ei]] * x_biobus # performance fraction

                    if ( cat=='LA' and fuel=='diesel'):
                      s = suorite.loc['LA diesel'][VTT_euros[ei]] * (1 - x_biobus) # performance

                    # Calculate emission factor per vehicle category, fuel and Euro standard
                    EF += E * fleet_distr[cat].loc[timestr:timestr2].sum() * s

                    # Calculate black and organic carbon from the PM emission:
                    if pol == 'PM Exhaust':
                      pp = partition[((partition['Category']         == cati) &
                                      (partition['Fuel']             == fueli) &
                                      (partition['Emission standard']== EEA_euro))]['BC/PM2.5']
                      if ( cat=='LA' and fuel==VTT_fuel[-1]):
                        pp = partition[((partition['Category']         == cati) &
                                        (partition['Fuel']             == 'Diesel') &
                                        (partition['Emission standard']== EEA_euro))]['BC/PM2.5']
                      if len(pp)>0:
                        BC += E * fleet_distr[cat].loc[timestr:timestr2].sum() * s * pp.values[0]

                      pp = partition[((partition['Category']         == cati) &
                                      (partition['Fuel']             == fueli) &
                                      (partition['Emission standard']== EEA_euro))]['OM/PM2.5']
                      if ( cat=='LA' and fuel==VTT_fuel[-1]):
                        pp = partition[((partition['Category']         == cati) &
                                        (partition['Fuel']             == 'Diesel') &
                                        (partition['Emission standard']== EEA_euro))]['OM/PM2.5']
                      if len(pp)>0:
                        OC += E * fleet_distr[cat].loc[timestr:timestr2].sum() * s * pp.values[0]

                    # NO and NO2 from NOx:
                    if pol == 'NOx':
                      pp = partition[((partition['Category']         == cati) &
                                      (partition['Fuel']             == fueli) &
                                      (partition['Emission standard']== EEA_euro))]['f-NO2']
                      if ( cat=='LA' and fuel==VTT_fuel[-1]):
                        pp = partition[((partition['Category']         == cati) &
                                        (partition['Fuel']             == 'Diesel') &
                                        (partition['Emission standard']== EEA_euro))]['f-NO2']
                      if len(pp)>0:
                        NO2 += E * fleet_distr[cat].loc[timestr:timestr2].sum() * s * pp.values[0]
                        NO += E * fleet_distr[cat].loc[timestr:timestr2].sum() * s * (1-pp.values[0])

                  # Sulphur dioxide:
                  if pol=='SO2':
                    r = FC[((FC['Category']==cati) & (FC['Fuel']==fueli))]
                    if len(seg)==1:
                      rr = r[r['Segment']==seg[0]]
                    else:
                      rr = r[((r['Segment']==seg[0]) | (r['Segment']==seg[1]))]
                    if len(rr)>0:
                      E = np.nanmean(2 * rr['sulphur_content'] * rr['FC (g/km)'])
                      EF += E * fleet_distr[cat].loc[timestr:timestr2].sum() * s

                  # Nitrous oxide and ammonia:
                  if (pol=='N2O' or pol=='NH3'):
                    r = FC[((FC['Category']==cati) & (FC['Fuel']==fueli))]
                    pp = partition[((partition['Category']         == cati) &
                                    (partition['Fuel']             == fueli) &
                                    (partition['Emission standard']== EEA_euro))].filter(regex=pol)
                    if len(seg)==1:
                      rr = r[ r['Segment']==seg[0] ]
                    else:
                      rr = r[ ( ( r['Segment']==seg[0] ) | ( r['Segment']==seg[1] ) ) ]
                    if len(rr)>0 and len(pp)>0:
                      ppi = pp
                      rri = rr
                      E = np.nanmean( rr['FC (g/km)'] * pp.values[0] )
                      EF += E * fleet_distr[cat].loc[timestr:timestr2].sum() * s

                  # Heat from traffic:
                  if pol=='E_heat': 
                    r = FC[((FC['Category']==cati) & (FC['Fuel']==fueli))]
                    if len(seg)==1:
                      rr = r[r['Segment']==seg[0]]
                    else:
                      rr = r[((r['Segment']==seg[0]) | (r['Segment']==seg[1]))]
                    if len(rr)>0:
                      E = np.nanmean(rr['Caloric value (J/kg)'] * rr['FC (g/km)'] * 1.0e-3) # J/km
                      EF += E * fleet_distr[cat].loc[timestr:timestr2].sum() * s # J/km

                  # Fuel consumption:
                  if pol=='FC': 
                    r = FC_VTT.loc['{} {}'.format(cat, fuel)][VTT_euros[ei]]
                    EF += np.nanmean(r) * fleet_distr[cat].loc[timestr:timestr2].sum() * s

                    r = FC[((FC['Category']==cati) & (FC['Fuel']==fueli))]
                    if len(seg)==1:
                      rr = r[r['Segment']==seg[0]]
                    else:
                      rr = r[((r['Segment']==seg[0]) | (r['Segment']==seg[1]))]
                    if len(rr)>0:
                      FCi += np.nanmean(rr['FC (g/km)']) * \
                             fleet_distr[cat].loc[timestr:timestr2].sum() * s

                fi += 1
              ci += 1

            # Print to screen and save to the output file
            print( pol, EF )
            EF_row += ', {:4.3e}'.format(EF / 1000.0)
            if pol=='PM Exhaust':
              EF_row += ', {:4.3e}'.format(BC / 1000.0)
              EF_row += ', {:4.3e}'.format(OC / 1000.0)
            if pol=='NOx':
              EF_row += ', {:4.3e}'.format(NO  / 1000.0)
              EF_row += ', {:4.3e}'.format(NO2 / 1000.0)
            if pol=='VOC':
              EF_row += ', {:4.3e}'.format(0.01 * EF / 1000.0)
              EF_row += ', {:4.3e}'.format(0.4  * EF / 1000.0)
            if pol=='FC':
              EF_row += ', {:4.3e}'.format(FCi / 1000.0)
              print(FCi / 1000.0)

          outrow = datarow+EF_row
          print(datarow+EF_row, file=f1)
          ti += 1
        else:
          if not hourly:
            outrow = datarow+EF_row
            print(datarow+EF_row, file=f1)
f1.close()
