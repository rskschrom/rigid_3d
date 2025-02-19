from rigidpy.pybind11_lib._freefallsimulation import FreeFallSimulation
from rigidpy.pybind11_lib._state import State
from rigidpy.pybind11_lib._particle import Particle

import numpy as np

# test ffsim creation
def test_int_ffsim():
    
    par = Particle('data/crystal_points.txt', 920.*(0.015*1.e-3)**3.)
    st = State(par.getMatInerm())

    ffsim = FreeFallSimulation(st, 10, 1.e-2, -9.81, 10., 920., 1.e-5)
    ffsim.evolveMotionInertial(1)

    orient = ffsim.orientHistory
    assert orient[0]==1.