import numpy as np
import matplotlib.pyplot as plt
from crystal_dda.crystal_dda import branched_planar_dda
from crystal_dda.polygons import make_branched_planar

# set values to create branched planar crystal with
a = 1.5
amax = 3.
ac = 0.2

fb = 0.3
ft = 0.2
fg = 0.7

nsb = 5

ag = fg*amax+(1.-fg)*ac
nxp = 30
nzp = 2

fname, afrac = branched_planar_dda(a, amax, ac, ft, fb, fg, nsb, nxp, nzp)
print(fname, afrac)

xp, yp = make_branched_planar(amax, ac, ft, fb, fg, nsb, 0.)

plt.plot(xp, yp, 'm-')
plt.savefig('poly.png')

'''
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
np.savetxt('par_data1.txt', np.c_[x2f, y2f],
           fmt=('%d','%d'), comments='')
'''
