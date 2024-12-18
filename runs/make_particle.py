from crystal_dda.polygons import BranchedPlanar
from crystal_dda.shapes import Crystal
import matplotlib.pyplot as plt
import numpy as np

# create branched planar crystal polygon
a_axis_max = 3.
a_axis = 1.5
a_axis_core = 0.2
frac_tip = 0.2
frac_bcov = 0.3
frac_gap = 0.7
nsub = 5
bp = BranchedPlanar(a_axis, a_axis_max, a_axis_core,
                    frac_tip, frac_bcov, frac_gap, nsub)

# create 3d crystal
c_axis = 0.2
cr = Crystal(bp, c_axis)

# fill with dipoles
dip_len = 0.1
cr.create_dipoles(dip_len)
#cr.write_dipoles('crystal.txt')
cr.write_points('crystal_points.txt', dip_len)
