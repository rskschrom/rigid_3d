from rigidpy.pybind11_lib._state import State
import numpy as np

# test state creation
def test_state():
    inerm = np.eye(3)
    s1 = State(inerm)
    assert np.sum((s1.getMatIInerm()-inerm)**2.)==0
    
