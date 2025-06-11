from rigidpy.freefallsimulation import FreeFallSimulation
from rigidpy.state import State
import numpy as np

# test ffsim creation
def test_ffsim():
    st = State(np.eye(3))

    ffsim = FreeFallSimulation(st, 10, 1.e-2, -9.81, 10., 920., 1.e-5)
    ffsim.evolve_motion_inertial(1)

    orient = ffsim.get_orient_history()
    assert orient[0][0]==1.