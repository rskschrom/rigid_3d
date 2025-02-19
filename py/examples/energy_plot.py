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

# create particle
dip_len = 0.015*1.e-3
rhob = 920.
par = Particle('../tests/data/crystal_points.txt', rhob*dip_len**3.)
par.set_omega((0.,0.,0.02))
par.set_orient_zenith(5.)

# calculate the rotational energies for orientations
g = -9.81
rhof_grad = 10.
damp_frac = 1.e-5
# step through simulation and examine energies
nstep = 20000
dt = 1.e-2
ffsim = FreeFallSimulation(par, 1, dt, g, rhof_grad, rhob, damp_frac)
pe = np.empty([nstep])
ke = np.empty([nstep])

#ffsim.evolve_motion_buoyancy(1)

for i in range(nstep):
    ffsim.evolve_motion_buoyancy(1)
    pe[i] = ffsim.par.rotational_potential_energy(g, rhof_grad, rhob)
    ke[i] = ffsim.par.rotational_kinetic_energy()
    print(ffsim.get_sim_step(), pe[i], ke[i], pe[i]+ke[i])

# plot
ts = np.arange(nstep)*dt
plt.plot(ts, ke, 'm-', label='KE')
plt.plot(ts, pe, 'c-', label='PE')
plt.plot(ts, pe+ke, 'k-', label='PE+KE')
plt.legend()
plt.savefig('energy.png')
