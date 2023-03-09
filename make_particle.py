import numpy as np
from crystal_dda.crystal_dda import branched_planar_dda

# set values to create branched planar crystal with
a = 3.
amax = 3.
ac = 0.3

fb = 0.3
ft = 0.2
fg = 0.4

nsb = 3

ag = fg*amax+(1.-fg)*ac
nxp = 50
nzp = 3

fname, afrac = branched_planar_dda(a, amax, ac, ft, fb, fg, nsb, nxp, nzp)
print(fname, afrac)

# plot dipole locations from text file
dda_data = np.genfromtxt(fname, skip_header=3)
x3d = dda_data[:,0]
y3d = dda_data[:,1]
z3d = dda_data[:,2]

# get 2d slice of x y values
zval = z3d[0]
x2d = x3d[z3d==zval]
y2d = y3d[z3d==zval]
x2f = x2d.flatten()
y2f = y2d.flatten()

# write 2d positions to text file
np.savetxt('par_data2.txt', np.c_[x2f, y2f],
           fmt=('%d','%d'), comments='')
