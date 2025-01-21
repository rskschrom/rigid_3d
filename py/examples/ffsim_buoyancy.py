from rigidpy.freefallsimulation import FreeFallSimulation
from rigidpy.particle import Particle
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

# create particle and do simulation
dip_len = 0.015*1.e-3
par = Particle('../tests/data/crystal_points.txt', 920.*dip_len**3.)
dt = 1.e-3

# loop over different initial zenith angles
zens = np.arange(5)*10.+5.
omegaz = np.arange(5)*0.1+0.1
bets = [None]*5

for i in range(5):
    print(zens[i])
    par.set_omega_body([0.,0.,0.2])
    par.set_orient_zenith(zens[i])
    ffsim = FreeFallSimulation(par, 20000, dt, -10.)
    ffsim.evolve_motion_buoyancy()
    orient = ffsim.get_orient_history()
    alp, bet = get_euler(orient)
    bets[i] = bet

# plot
nts = len(bet)
sf = 5
ts = sf*np.arange(nts)*dt
print(ts.shape)
lns = ['k-','r-','c-','m-','y-']

for i in range(5):
    bet = bets[i]
    print(bet[::sf].shape, bet.shape, ts[::sf].shape)
    plt.plot(ts[::sf], np.cos(bet[::sf]*np.pi/180.), lns[i])
plt.savefig('euler_ts.png')

