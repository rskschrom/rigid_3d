from rigidpy.pybind11_lib._particle import Particle
import numpy as np

# test particle
def test_particle():
    # create points
    dx = 1.e-4
    xp = np.arange(51)*dx
    yp = np.zeros([51])
    zp = np.zeros([51])
    
    # subtract mean
    xp = xp-np.mean(xp)
    
    points = np.stack((xp,yp,zp), axis=1)
    points = points.flatten()
    
    par = Particle(points, 920.*(dx)**3.)
    inerm = par.inertiaMomentTensor()
    inerm_true = np.array([[2.77517007e-23,0.00000000e+00,0.00000000e+00],
                           [0.00000000e+00,1.01688386e-13,0.00000000e+00],
                           [0.00000000e+00,0.00000000e+00,1.01659993e-13]])
    assert np.sum((inerm-inerm_true)**2.)<1.e-16
