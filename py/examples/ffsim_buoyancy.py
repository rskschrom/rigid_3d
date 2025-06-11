from rigidpy.freefallsimulation import FreeFallSimulation
from rigidpy.particle import Particle
from rigidpy.state import State
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
rho = 920.
par = Particle.from_file('../tests/data/crystal_remesh.stl', rho)
st = State(par.get_mat_inerm())
dt = 1.e-2
nstep = 20000

# loop over different initial zenith angles
zens = np.arange(5)*10.+5.
omegaz = np.arange(5)*0.1+0.1
bets = [None]*5

for i in range(5):
    print(zens[i])
    st.set_omega([0.,0.,0.2])
    st.set_orient_zenith(zens[i])
    ffsim = FreeFallSimulation(st, nstep, dt, -9.81, 10., rho, 1.e-5)
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

