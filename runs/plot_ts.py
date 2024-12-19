import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation

# get euler angles
def get_euler(orient):
    orient = orient[:,[1,2,3,0]]
    r = Rotation.from_quat(orient)
    euler_angs = r.as_euler('ZYZ', degrees=True)
    alpha = euler_angs[:,0]
    beta = euler_angs[:,1]
    return alpha, beta

# read orientation data
orient_old = np.genfromtxt('orientation_old.txt')
orient = np.genfromtxt('orientation.txt')
nts = int(len(orient)/4)
orient.shape = (nts,4)
orient_old.shape = (nts,4)

# get polar angle
alp, bet = get_euler(orient)
alp_old, bet_old = get_euler(orient_old)

# plot
ts = np.arange(nts)*0.01
sf = 5

#plt.plot(ts[::sf], np.cos(alpha[::sf]*np.pi/180.), 'k-')
plt.plot(ts[::sf], np.cos(bet_old[::sf]*np.pi/180.), 'b-')
plt.plot(ts[::sf], np.cos(bet[::sf]*np.pi/180.), 'r-')
plt.savefig('euler_ts.png')
