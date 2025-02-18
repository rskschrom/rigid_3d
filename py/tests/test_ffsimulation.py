from rigidpy.pybind11_lib._freefallsimulation import FreeFallSimulation
from rigidpy.pybind11_lib._state import State
import numpy as np

# test ffsim creation
def test_ffsim():
    st = State(np.eye(3))

    ffsim = FreeFallSimulation(st, 10, 1.e-2, -9.81, 10., 920., 1.e-5)
    ffsim.evolveMotionInertial(1)

    orient = ffsim.orientHistory
    assert orient[0]==1.