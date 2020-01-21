#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

def_fs = 11

mpl.rcParams['lines.linewidth'] = 2
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.labelsize'] = def_fs
mpl.rcParams['ytick.labelsize'] = def_fs
mpl.rcParams['xtick.labelsize'] = def_fs
mpl.rcParams['grid.linewidth'] = 1
mpl.rcParams['lines.markersize'] = 4.5
mpl.rcParams['legend.fontsize'] = def_fs
#mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
#mpl.rc('text', usetex=True)

# - - - - - - - - - - - - - - - - - - - - - - - #
def define_bins( nbin, reglim ):

  nbins = np.sum( nbin ) # = subrange 1 + subrange 2

  # Log-normal to sectional
  
  vlolim = np.zeros( nbins )
  vhilim = np.zeros( nbins )
  dmid   = np.zeros( nbins )
  bin_limits = np.zeros( nbins )

  # Sectional bin limits
  ratio_d = reglim[1] / reglim[0]
  for b in range( nbin[0] ):
    vlolim[b] = np.pi / 6.0 * ( reglim[0] * ratio_d **( float(b) / nbin[0] ) )**3
    vhilim[b] = np.pi / 6.0 * ( reglim[0] * ratio_d **( float(b+1) / nbin[0] ) )**3
    dmid[b] = np.sqrt( ( 6.0 * vhilim[b] / np.pi )**0.33333333 * ( 6.0 * vlolim[b] / np.pi )**0.33333333 )                 

  ratio_d = reglim[2] / reglim[1]
  for b in np.arange( nbin[0], np.sum( nbin ),1 ):
    c = b-nbin[0]
    vlolim[b] = np.pi / 6.0 * ( reglim[1] * ratio_d ** ( float(c) / nbin[1] ) )**3
    vhilim[b] = np.pi / 6.0 * ( reglim[1] * ratio_d ** ( float(c+1) / nbin[1] ) ) ** 3
    dmid[b] = np.sqrt( ( 6.0 * vhilim[b] / np.pi )**0.33333333 * ( 6.0 * vlolim[b] / np.pi )**0.33333333 )

  bin_limits = ( 6.0 * vlolim / np.pi )**0.33333333
  bin_limits = np.append( bin_limits, reglim[-1] )

  return dmid, bin_limits

# - - - - - - - - - - - - - - - - - - - - - - - #

def psd_from_data( input_file, EF, bin_limits, input_type, plot ):
  '''
  psd_from_data( input_file, EF, nbins, input_type )
   input_file (filename)
   EF (emission factor in g/m2/s or #/m2/s)
   nbins (number of size bins)
   input_type ('mass' or 'number')
  '''
  
  # Fit manually:
  dpg   = np.array([    4.0,    13.0,   75.0, 220.0]) # in nm
  sigma = np.array([   1.45,    2.00,   1.60,  1.45])
  rho   = 2000.0
  
  if input_type == 'mass':
    V    = np.array([0.00015*EF, 0.008*EF, 0.3*EF, 0.22*EF ])
    V    = V / ( rho * 1e3 ) *1e12 
    n = 6.0 * V / (np.pi * (1.0e-3*dpg)**3 * np.exp(9.0/2.0 * np.log(sigma)**2) )
  elif input_type == 'number':
    n = np.array([11000.0,  4000.0, 2500.0, 130.0]) # in 1/cm3
    n = n / np.sum( n ) * EF
  
    
  # Read in data
    
  data = np.genfromtxt( input_file, delimiter=',' )
  
  # col 0: D (nm)
  # col 1: dN/dlogD (1/cm3)
                        
  dmid = np.sqrt( bin_limits[0:-1] * bin_limits[1::] )    
  
  # Calculate size distributions:
  
  # Log-normal 
  d, psd_lognorm = create_lognorm( dpg / 1000., n, sigma )
  
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
  
  return bin_limits, relative_per_bin, psd_sect

# - - - - - - - - - - - - - - - - - - - - - - - #
  
def create_sectional(gmD, N, sigma, bin_lims):
# gmD = geometric mean diameter of a mode (nm)
# N = total number concentration of a mode (1/cm3)
# sigma = standard deviation of a mode
# bin_lims = bin limits ofr section representation
# nsect = number concentration per bin (1/cm3)

 # if ( gmD.size!=N.size ) or ( N.size!=sigma.size ) or ( gmD!=sigma.size ): 
 #  exit('Error.')
  nsect = np.zeros(len(bin_lims)-1,dtype=float)
  for l in range(1,len(bin_lims)): 
    d1 = bin_lims[l-1]
    d2 = bin_lims[l]
    delta_d = (d2-d1)/10.0
    for ib in range(1,11):
      d1 = bin_lims[l-1] + (ib-1)*delta_d
      d2 = d1+delta_d
      dmidi = (d1+d2)/2.0
      deltadp = np.log10(d2/d1)
      nsect[l-1] = nsect[l-1] + np.sum( 1.0e6 * N  * deltadp / ( np.sqrt( 
         2.0 * np.pi ) * np.log10( sigma ) ) * np.exp( -np.log10( dmidi 
         / ( 1.0E-9 * gmD ) )**2.0 / ( 2.0 * np.log10( sigma ) ** 2.0 ) ) )
  nsect_per_dlogD = nsect / np.log10(bin_lims[1::]/bin_lims[0:-1])
  return nsect_per_dlogD, nsect
  
# - - - - - - - - - - - - - - - - - - - - - - - #  
  
def create_lognorm(gmD, N, sigma):
  x = np.linspace( 0.000001, 10.0, 100000 )
  n_N_sum = np.zeros((len(x)),dtype=float) 
  if np.shape(gmD)==():
    n_N = N / ( np.sqrt(2*np.pi) * np.log10(sigma) ) * \
                np.exp( -( np.log10(x) - np.log10(gmD) )**2 / \
                ( 2.0 * np.log10(sigma)**2 ) ) 
    if np.sum(n_N>0.0):
      n_N_sum += n_N
  else:
    for i in range(0, gmD.size):
      n_N = N[i] / ( np.sqrt(2*np.pi) * np.log10(sigma[i]) ) * \
                np.exp( -( np.log10(x) - np.log10(gmD[i]) )**2 / \
                ( 2.0 * np.log10(sigma[i])**2 ) ) 
      if np.sum(n_N>0.0):
        n_N_sum += n_N
  return x*1000., n_N_sum  
 
