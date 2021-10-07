'''
flayer_ur.py

Track the phi and theta averaged radial velocity within the F-layer to rule out numerically 'gravitational-like' instabilities.
'''

import os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import h5py

from paropy.data_utils import load_dimensionless
from paropy.routines import ur_avg

matplotlib.use('Agg')  # backend for no display

#%%--------------------------------------------------------------------------%%
# INPUT PARAMETERS
#----------------------------------------------------------------------------%%
# PARODY simulation tag
run_ID = ['d_0_55a', 'd_0_6a', 'd_0_65a', 'c-200a', 'd_0_75a', 'd_0_8a']
# run_ID = ['d_0_6a', 'd_0_65a', 'c-200a', 'd_0_75a']
# path containing runs
# directory = '/Volumes/NAS/ipgp/Work/{}'.format(run_ID)

fig_aspect = 1  # figure aspect ratio
yLim = 1

saveOn = 1  # save figures?
saveDir = '/home/wongj/Work/figures/flayer_ur'  # path to save files
# saveDir = '/Users/wongj/Documents/isterre/parody/figures/flayer_ur'

#%%----------------------------------------------------------------------------
w, h = plt.figaspect(fig_aspect)
fig1, ax1 = plt.subplots(1, 1, figsize=(1.5*w, h))

# Perform ur average in F-layer
for run in run_ID:
    directory = '/data/geodynamo/wongj/Work/{}'.format(run)
    _, _, _, _, _, fi, rf = load_dimensionless(run, directory)
    if not os.path.exists('{}/flayer_ur'.format(directory)):
        ur = ur_avg(run, directory, rf)  # timeavg
    else:  # load flayer ur data
        print('Loading {}/flayer_ur'.format(directory))
        f = h5py.File('{}/flayer_ur'.format(directory), 'r')
        for key in f.keys():
            globals()[key] = np.array(f[key])

    # Plot
    ax1.plot(ur, label = "{}".format(run))

ax1.set_ylabel(r"$\langle u_r \rangle_{\theta,\phi}$")
ax1.set_ylim([-yLim, yLim])
ax1.set_yscale("symlog")
ax1.legend()

# Save
if saveOn == 1:
    if not os.path.exists('{}'.format(saveDir)):
        os.makedirs('{}'.format(saveDir))
    fig1.savefig('{}/compare.png'.format(saveDir),
                 format='png', dpi=200, bbox_inches='tight')
    fig1.savefig('{}/compare.pdf'.format(saveDir),
                 format='pdf', dpi=200, bbox_inches='tight')
    print('Figures saved as {}/compare.*'.format(saveDir))
