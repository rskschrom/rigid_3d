import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

# read orientation data
orient = np.genfromtxt('orientation.txt')
nts = int(len(orient)/4)
orient.shape = (nts,4)

# get polar angle
orient = orient[:,[1,2,3,0]]
r = Rotation.from_quat(orient)
euler_angs = r.as_euler('ZYZ', degrees=True)
alpha = euler_angs[:,0]
beta = euler_angs[:,1]

# plot
ts = np.arange(nts)*0.01
sf = 5

print(alpha[::sf])

plt.plot(ts[::sf], np.cos(alpha[::sf]*np.pi/180.), 'k-')
plt.plot(ts[::sf], np.cos(beta[::sf]*np.pi/180.), 'r-')
plt.savefig('euler_ts.png')
