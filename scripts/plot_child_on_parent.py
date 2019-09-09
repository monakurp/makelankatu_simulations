import numpy as np
import matplotlib.pyplot as plt
from scipy import ndimage
import os
import matplotlib as mpl
import os
from mpl_toolkits.axes_grid1 import make_axes_locatable

def_fs = 11
zoom_in = True

mpl.rcParams['lines.linewidth'] = 1.5
mpl.rcParams['axes.linewidth'] = 1
mpl.rcParams['axes.labelsize'] = def_fs
mpl.rcParams['ytick.labelsize'] = def_fs
mpl.rcParams['xtick.labelsize'] = def_fs 
mpl.rcParams['grid.linewidth'] = 1
mpl.rcParams['lines.markersize'] = 6
mpl.rcParams['legend.fontsize'] = def_fs
mpl.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
mpl.rc('text', usetex=True)

folder = '/home/monakurp/makelankatu_simulations/'

os.chdir( ( folder ) )

plt.close('all')

file_root = 'input_data_to_palm/topo_root.npz'
file_parent = 'input_data_to_palm/topo_parent.npz'
file_child = 'input_data_to_palm/topo_child.npz'
fname_plot = 'input_data_to_palm/topos.pdf'

root = np.load(file_root)
parent = np.load(file_parent)
child = np.load(file_child)

for k in parent.files:
  print(k)
  
# root to parent resolution:  
ratio_rp = root['dPx'][0]/parent['dPx'][0]
dim_parent_in_root = [int(ratio_rp)*i for i in np.shape( parent['R'] ) ]
plot_parent = np.zeros( dim_parent_in_root, dtype=float )+ np.nan
shift = ( ( root['GlobOrig'] - parent['GlobOrig'] ) / parent['dPx'][0] ).astype(int)
dim_parent = np.shape( parent['R'] )
plot_parent[ shift[0]:shift[0]+dim_parent[0], -shift[1]:-shift[1]+dim_parent[1] ] = parent['R']
plot_root = ndimage.zoom( root['R'], ratio_rp )

# child to parent resolution
ratio = parent['dPx'][0]/child['dPx'][0]
dim_child_in_root = dim_parent_in_root
plot_child = np.zeros( dim_child_in_root, dtype=float )+ np.nan
shift2 = ( ( root['GlobOrig'] - child['GlobOrig'] ) / ratio ).astype(int)
dim_child = np.shape( child['R'] )
plot_childi = ndimage.zoom( child['R'], 1/ratio )
x1 = shift2[0]
x2 = np.ceil(shift2[0]+(dim_child[0]/ratio)).astype(int)
y1 = -shift2[1]
y2 = np.ceil(-shift2[1]+dim_child[1]/ratio).astype(int)
plot_child[x1:x2,y1:y2] = plot_childi

fig = plt.figure(figsize=(6,6), dpi=100)  
ax = fig.add_subplot(111)
ax.imshow(plot_root, interpolation='None', cmap='binary', vmin=0, vmax=50, alpha=0.7)
ax.imshow(plot_parent, interpolation='None', cmap='viridis_r', vmin=0, vmax=50, alpha=0.8)
ax.imshow(plot_child, interpolation='None', cmap='rainbow')
#plt.plot(dim_parent_in_root[0]/2,dim_parent_in_root[1]/2, 'r*', markersize=10)

rect = mpl.patches.Rectangle((0,0), 2500.0/parent['dPx'][0], 2500.0/parent['dPx'][0], linewidth=1, 
                              edgecolor='none', facecolor='red',alpha=0.5 )
ax.add_patch(rect)

ax.plot([2500.0/parent['dPx'][0], 2500.0/parent['dPx'][0]],[0,dim_parent_in_root[0]+1],'k--')
ax.plot([5000.0/parent['dPx'][0], 5000.0/parent['dPx'][0]],[0,dim_parent_in_root[0]+1],'k--')
ax.plot([0,dim_parent_in_root[0]+1],[2500.0/parent['dPx'][0], 2500.0/parent['dPx'][0]],'k--')
ax.plot([0,dim_parent_in_root[0]+1],[5000.0/parent['dPx'][0], 5000.0/parent['dPx'][0]],'k--')

xticksi = np.arange(0,8000.0,1250.0)/parent['dPx'][0]
ax.set_xticks(xticksi)
ax.set_yticks(xticksi)
xticks_labels = ['0 m','','2500 m','','5000 m','']
ax.set_xticklabels(xticks_labels)
ax.set_yticklabels(xticks_labels)
ax.axis([0,dim_parent_in_root[0],dim_parent_in_root[1],0])

ax.text(450.0/parent['dPx'][0], 1400.0/parent['dPx'][0], '  Harmonie grid \n 2.5 km x 2.5 km', color='w')  

fig.savefig( fname_plot, format='pdf', dpi=300 ) 
  
#ratio = parent['dPx'][0]/child['dPx'][0]
#dim_child_in_parent = [int(ratio)*i for i in np.shape( parent['R'] ) ]
#plot_child = np.zeros( dim_child_in_parent, dtype=float )+ np.nan
#shift = ( parent['GlobOrig'] - child['GlobOrig'] ).astype(int)
#dim_child = np.shape( child['R'] )
#
#plot_child[ shift[0]:shift[0]+dim_child[0], -shift[1]:-shift[1]+dim_child[1] ] = child['R']
#
#plot_parent = ndimage.zoom( parent['R'], ratio )

