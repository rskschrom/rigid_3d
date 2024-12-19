#from rigidpy.particle import Particle
from rigidpy.pybind11_lib._particle import Particle
from rigidpy.pybind11_lib._freefallsimulation import FreeFallSimulation

# test ffsim
par = Particle('data/crystal_points.txt', 920.*(0.1*1.e-3)**3.)
par.setOmega([0.5,0.,0.2])

ffsim = FreeFallSimulation(par, 2500, 1.e-3, -6.)
print(dir(ffsim))
ffsim.evolveMotionInertial(100)

orient = ffsim.orientHistory
print(orient)
#print(dir(par))
#print(par.get_rel_points())
