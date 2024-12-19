from .pybind11_lib._freefallsimulation import FreeFallSimulation as _FreeFallSimulation
from .particle import Particle
import numpy as np

class FreeFallSimulation():
    '''
    FreeFallSimulation class
    '''
    def __init__(self, par, nt, dt, g):
        self.par = par
        self.nt = nt
        self.dt = dt
        self.g = g
        self._FreeFallSimulation = _FreeFallSimulation(par._Particle, nt, dt, g)
        
    def get_pos_history(self):
        pos_history = np.array(self._FreeFallSimulation.posHistory)
        pos_history.shape = (int(pos_history.shape[0]/3),3)
        
        return pos_history
        
    def get_orient_history(self):
        orient_history = np.array(self._FreeFallSimulation.orientHistory)
        orient_history.shape = (int(orient_history.shape[0]/4),4)
        
        return orient_history
        
    def evolve_motion_inertial(self, nstep=None):
        if nstep is None:
            nstep = self.nt
        self._FreeFallSimulation.evolveMotionInertial(nstep)
        return
        
    def evolve_motion_buoyancy(self, nstep=None):
        if nstep is None:
            nstep = self.nt
        self._FreeFallSimulation.evolveMotionBuoyancy(nstep)
        return
    
