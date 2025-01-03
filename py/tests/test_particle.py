from rigidpy.particle import Particle
#from rigidpy.pybind11_lib._particle import Particle

# test particle
par = Particle('data/crystal_points.txt', 920.*(0.1*1.e-3)**3.)

print(dir(par))
print(par.get_rel_points())
