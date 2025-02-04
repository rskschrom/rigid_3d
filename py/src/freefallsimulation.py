from .pybind11_lib._freefallsimulation import FreeFallSimulation as _FreeFallSimulation
from .particle import Particle
import numpy as np

class FreeFallSimulation():
    '''
    FreeFallSimulation class
    '''
    def __init__(self, par, nt, dt, g, rhof_grad, rhob):
        self.par = par
        self.nt = nt
        self.dt = dt
        self.g = g
        self.rhof_grad = rhof_grad
        self.rhob = rhob
        self._FreeFallSimulation = _FreeFallSimulation(self.par._Particle, nt, dt, g, rhof_grad, rhob)
        
    def get_pos_history(self):
        '''
        Return position history during simulation.
        '''
        pos_history = np.array(self._FreeFallSimulation.posHistory)
        pos_history.shape = (int(pos_history.shape[0]/3),3)
        
        return pos_history
        
    def get_orient_history(self):
        orient_history = np.array(self._FreeFallSimulation.orientHistory)
        orient_history.shape = (int(orient_history.shape[0]/4),4)
        
        return orient_history
    
    def get_sim_step(self):
        sim_step = self._FreeFallSimulation.getSimStep()
        return sim_step
        
    def evolve_motion_inertial(self, nstep=None):
        if nstep is None:
            nstep = self.nt
        self._FreeFallSimulation.evolveMotionInertial(nstep)
        
        # update python particle object with wrapped pybind11 particle object
        self.par._Particle = self._FreeFallSimulation.par
        return
        
    def evolve_motion_buoyancy(self, nstep=None):
        if nstep is None:
            nstep = self.nt
        self._FreeFallSimulation.evolveMotionBuoyancy(nstep)
        
        # update python particle object with wrapped pybind11 particle object
        self.par._Particle = self._FreeFallSimulation.par
        return
    
