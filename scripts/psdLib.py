#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

# - - - - - - - - - - - - - - - - - - - - - - - #
def define_bins(nbin, reglim):
  """
  Define the geometric mean diameters for the sectional aerosol size distribution

  Input:
    nbin = number of bins per subrange (e.g. [2,8])
    reglim = subrange diameter limits (m) (e.g. [2.5e-9, 10.0e-9, 2.5e-6])

  Output:
    dmid = geometric mean diameter per bin (m)
    bin_limits = bin limits of the sectional representation (m)
  """

  # Total number of bins
  nbins = np.sum(nbin) # = subrange 1 + subrange 2

  # Initialise arrays
  vlolim = np.zeros(nbins)
  vhilim = np.zeros(nbins)
  dmid   = np.zeros(nbins)
  bin_limits = np.zeros(nbins)

  # Define the geometric mean diameters
  # subrange 1:
  ratio_d = reglim[1] / reglim[0]
  for b in range(nbin[0]):
    vlolim[b] = np.pi / 6.0 * (reglim[0] * ratio_d **(float(b) / nbin[0]))**3
    vhilim[b] = np.pi / 6.0 * (reglim[0] * ratio_d **(float(b+1) / nbin[0]))**3
    dmid[b] = np.sqrt((6.0 * vhilim[b] / np.pi)**0.33333333 * \
                      (6.0 * vlolim[b] / np.pi )**0.33333333)
  # subrange 2:
  ratio_d = reglim[2] / reglim[1]
  for b in np.arange(nbin[0], nbins,1):
    c = b-nbin[0]
    vlolim[b] = np.pi / 6.0 * (reglim[1] * ratio_d ** (float(c) / nbin[1]))**3
    vhilim[b] = np.pi / 6.0 * (reglim[1] * ratio_d ** (float(c+1) / nbin[1]))**3
    dmid[b] = np.sqrt((6.0 * vhilim[b] / np.pi)**0.33333333 * \
                      (6.0 * vlolim[b] / np.pi)**0.33333333)

  # Bin limits
  bin_limits = (6.0 * vlolim / np.pi)**0.33333333
  bin_limits = np.append(bin_limits, reglim[-1])

  return dmid, bin_limits

# - - - - - - - - - - - - - - - - - - - - - - - #

def psd_from_data(input_file, EF, bin_limits, input_type, plot):
  """
  Plot size distribution data and a sectional size distribution.

  Input:
    input_file = filename for the aerosol number size distribution
    EF = emission factor in g/m/veh or #/m/veh
    bin_limits= bin limits of the sectional representation (m)
    input_type = 'mass' or 'number'
    plot = True or False

  Output:
    relative_per_bin = relative concentration per bin (sum equals 1)
    psd_sect = number emission per bin #/m/s
  """

  # Fit manually:
  dpg   = np.array([4.0 , 13.0, 75.0, 220.0]) * 1e-9 # in m
  sigma = np.array([1.45, 2.00, 1.60,  1.45])
  rho   = 2000.0 * 1e3 # g/m3

  if input_type == 'mass':
    m = np.array([0.00015*EF, 0.008*EF, 0.3*EF, 0.22*EF ])
    V = m / rho # V = rho/m (m3)
    n = 6.0 * V / (np.pi * dpg**3 * np.exp(9.0/2.0 * np.log(sigma)**2)) # 1/m3
  elif input_type == 'number':
    n = np.array([11000.0,  4000.0, 2500.0, 130.0]) # in 1/cm3
    n = n / np.sum(n) * EF * 1e6 # 1/m3

  # Read in data
  data = np.genfromtxt( input_file, delimiter=',' )

  # col 0: D (nm)
  # col 1: dN/dlogD (1/cm3)

  dmid = np.sqrt( bin_limits[0:-1] * bin_limits[1::] )

  # Calculate size distributions:

  # Log-normal
  d, psd_lognorm = create_lognorm( dpg, n, sigma )

  # Sectional
  [psd_sect_per_logD, psd_sect] = create_sectional( dpg, n, sigma, bin_limits )

  # Plot data and fit
  if plot:
    fig = plt.figure()
    ax  = fig.add_subplot(111)
    ax.loglog( data[:,0], data[:,1], 'r-', label='Measurements' )
    ax.loglog( d, psd_lognorm, 'b--', label='Log-normal' )
    ax.loglog( dmid*1e9, psd_sect_per_logD*1e-6, 'k:o', label='Sectional' )
    ax.set_xlim(1.0,2500.0)
    ax.set_ylim(1e2,1e5)
    ax.set_xlabel('D (nm)')
    ax.set_ylabel('$dN/dlogD\,\mathrm{(1/cm^3)}$')
    ax.grid(True)
    plt.legend()
    plt.show()

  relative_per_bin = psd_sect / np.sum( psd_sect ) 

  return relative_per_bin, psd_sect

# - - - - - - - - - - - - - - - - - - - - - - - #
  
def create_sectional(gmD, N, sigma, bin_limits):
  """
  Input:
    gmD = geometric mean diameter of a mode (m)
    N = total number concentration of a mode (1/m3)
    sigma = standard deviation of a mode
    bin_limits = bin limits of the sectional representation (m)
  
  Output:
    nsect = number concentration per bin (1/m3)
    nsect_per_dlogD = number concentration per bin (1/m3)(dN/dlogD)
  """

  nsect = np.zeros(len(bin_limits)-1,dtype=float)
  for l in range(1,len(bin_limits)): 
    d1 = bin_limits[l-1]
    d2 = bin_limits[l]
    delta_d = (d2-d1)/10.0
    for ib in range(1,11):
      d1 = bin_limits[l-1] + (ib-1)*delta_d
      d2 = d1+delta_d
      dmidi = (d1+d2)/2.0
      deltadp = np.log10(d2/d1)
      nsect[l-1] = nsect[l-1] + np.sum( N  * deltadp / ( np.sqrt( 
         2.0 * np.pi ) * np.log10( sigma ) ) * np.exp( -np.log10( dmidi 
         / gmD )**2.0 / ( 2.0 * np.log10( sigma ) ** 2.0 ) ) )
  nsect_per_dlogD = nsect / np.log10(bin_limits[1::]/bin_limits[0:-1])
  return nsect_per_dlogD, nsect
  
# - - - - - - - - - - - - - - - - - - - - - - - #
  
def create_lognorm(gmD, N, sigma):
  """
  Input:
    gmD = geometric mean diameter of a mode (m)
    N = total number concentration of a mode (1/m3)
    sigma = standard deviation of a mode
    
  Output: 
    x = diameter in nm
    n_N_sum = number concentration per bin (1/m3) (dN/dlogD)
  """
  
  x = np.linspace( 0.000001, 10.0, 100000 )
  n_N_sum = np.zeros((len(x)),dtype=float) 
  if np.shape(gmD)==():
    n_N = (N*1e-6) / ( np.sqrt(2*np.pi) * np.log10(sigma) ) * \
                np.exp( -( np.log10(x) - np.log10(gmD*1e6) )**2 / \
                ( 2.0 * np.log10(sigma)**2 ) ) 
    if np.sum(n_N>0.0):
      n_N_sum += n_N
  else:
    for i in range(0, gmD.size):
      n_N = (N[i]*1e-6) / ( np.sqrt(2*np.pi) * np.log10(sigma[i]) ) * \
                np.exp( -( np.log10(x) - np.log10(gmD[i]*1e6) )**2 / \
                ( 2.0 * np.log10(sigma[i])**2 ) ) 
      if np.sum(n_N>0.0):
        n_N_sum += n_N
  return x*1000., n_N_sum  

# - - - - - - - - - - - - - - - - - - - - - - - #
  
def interpolate_psd_from_figure( input_file, bin_limits, plot ):
  """
  Input:
    input_file = filename for the aerosol number size distribution
    bin_limits= bin limits of the sectional representation (m)
    plot = True or False
  
  Output:
    relative_per_bin = relative concentration per bin (sum equals 1)
    nsect = number concentration per bin (1/m3)
  """
  
  from scipy.interpolate import interp1d
  
  # Read in data 
  data = np.genfromtxt( input_file, delimiter=',' ) # columns: diamater, dNdlogD
  
  # Geometric and arithemic mean diameters
  dmid = np.sqrt( bin_limits[0:-1] * bin_limits[1::] )
  dmean = 0.5 * ( bin_limits[0:-1] + bin_limits[1::] )
  
  # Define the function to interpolate
  x = data[:,0]*1e-9
  y = data[:,1]
  f = interp1d( x, y )
  
  # Select diameters that fall within the range of the data
  within = ( dmid>=x[0] ) & ( dmid<=x[-1] )
  first_out = np.where( within==False )[0][0]
  xnew = dmid[ within ]
  ynew = f( xnew )
  
  # Calculate the absolute and relative number concentration per size bin:
  nsect = np.zeros( len( dmid ) )
  nsect[within] = ynew 
  nsect = nsect / np.log10( bin_limits[1::] / bin_limits[0:-1] ) * 1e6 # in #/m3
  relative_per_bin = nsect / np.sum( nsect )
  
  if plot:
    fig, ax = plt.subplots()
    ax.plot( x, y, '--')
    ax.plot( xnew, ynew, '*r')
    wbar = ( bin_limits[1:first_out+1]-bin_limits[0:first_out] )
    ax.bar( dmean[within], ynew, width=wbar, ec="k", alpha=0.3 ) 
    ax.legend( ['Data', 'Linear interpolation', 'Sectional representation'] )
    ax.set_xscale( 'log' )
    ax.set_yscale( 'log' )
    ax.set_xlabel( '$D_{mid}$ (m)' )
    ax.set_ylabel( '$dN/d\logD$ (cm$^{-3}$)' )
    
  return nsect, relative_per_bin
  
# - - - - - - - - - - - - - - - - - - - - - - - #
  
def ntot_from_psd_and_pm( pm, psd, dmid, rho=2.0e6 ):
  """
  Input:
    pm = particulate mass concentration or emission (1/m3 or e.g. #/m/veh) 
    psd = aerosol NUMBER size distribution (no need to be normalised)
    dmid = geometric mean diameter per bin (m)
    rho = aerosol particle density (g/m3)
    
  Output:
    ntot = total number concentration or emission (1/m3 or e.g. #/m/veh)
  """
  
  psd = psd / np.sum( psd )
  
  ntot = 6 * pm / ( np.pi * rho * np.sum( psd*dmid**3 ) )
  
  return ntot
  
  
