
'''
Find antenna position offsets in the appropriate frame
between two baseline solutions.

Goal is  to create an output that CASA's gencal can interpret
to correct baselines within CASA

gencal(vis=msname, caltable=f'{msname[:-3]}.antpos', caltype='antpos', antenna="",
      parameter=[])
'''

import numpy as np
from pyuvdata import get_telescope
from pyuvdata.utils import ECEF_from_rotECEF

orig_ants = np.loadtxt("antennas_orig")

new_ants = np.loadtxt("antennas_230104")


lat, lon, alt = get_telescope("SMA")._telescope_location.lat_lon_alt()


orig_antenna_positions = ECEF_from_rotECEF(orig_ants[:, 1:], lon)
new_antenna_positions = ECEF_from_rotECEF(new_ants[:, 1:], lon)

diff_positions = new_antenna_positions - orig_antenna_positions

# DO I HAVE THIS BACKWARDS? TESTING MODE:
# AS OF 01/05/2023, yes I believe this was backwards. The 01/04 VEX tracks
# has noticeably fewer phase gain wraps with the -1.
# Leaving like this for now so I remember, but this means we should be defining the
# change as ORIG - NEW positions.
diff_positions *= -1.

antenna_string = ''
xyz_diff_list = []

# Iterate over antennas
for ii in range(orig_ants.shape[0]):
    # Skip over any missing antennas: Skip the antenna on the reference pad (6)
    if (orig_ants[ii, 1:] == 0.).all():
        continue

    if (diff_positions[ii] != 0.).any():
        antenna_string += f"{ii+1},"
        xyz_diff_list.extend(list(diff_positions[ii]))

if antenna_string[-1] == ",":
    antenna_string = antenna_string[:-1]

# Now feed into gencal in mode='antpos'
# antenna=antenna_string
# parameter=xyz_diff_list

# In [38]: antenna_string
# Out[38]: '1,3,4,5,7,8'

# In [39]: xyz_diff_list
# Out[39]:

# [7.802013754565351e-05,
#  1.6907339741578653e-05,
#  -0.00015199999999992997,
#  -0.0006467228456088492,
#  2.2617713653971805e-05,
#  0.00014499999999995072,
#  -0.0005035325668423241,
#  7.255311257381436e-05,
#  0.00033100000000274576,
#  8.905469954356704e-06,
#  2.3847276679589413e-05,
#  -6.099999998809835e-05,
#  -0.0001734889423925523,
#  1.6479892877896418e-05,
#  0.00021200000000476393,
#  0.00016393552614601958,
#  -0.00017142095340716423,
#  -0.00016300000000057935]

