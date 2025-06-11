from .pybind11_lib._freefallsimulation import FreeFallSimulation as _FreeFallSimulation
from .pybind11_lib._state import State as _State
from .state import State
import numpy as np

class FreeFallSimulation():
    '''
    Class for simulating the kinematics of a body in freefall with a buoyancy gradient torque.
    
    Parameters
    ----------
    st : `State`
        The associated state of a body.
    nt : int
        The number of simulation time steps.
    dt : float
        The simulation timestep.
    rhof_grad : float
        The buoyancy gradient.
    rhob : float
        The density of the body.
    damp_frac : float
        The simulation damping fraction value.
        
    Returns
    -------
    The `FreeFallSimulation` object.
    '''
    def __init__(self, st, nt, dt, g, rhof_grad, rhob, damp_frac):
        self.st = st
        self.nt = nt
        self.dt = dt
        self.g = g
        self.rhof_grad = rhof_grad
        self.rhob = rhob
        self.damp_frac = damp_frac
        self._FreeFallSimulation = _FreeFallSimulation(self.st._State, nt, dt, g,
                                                       rhof_grad, rhob, damp_frac)
        
    def get_pos_history(self):
        '''
        Return position history during simulation.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        '''
        pos_history = np.array(self._FreeFallSimulation.posHistory)
        pos_history.shape = (int(pos_history.shape[0]/3),3)
        
        return pos_history
        
    def get_orient_history(self):
        '''
        Return orientation history during simulation.
        
        Parameters
        ----------
        None
        
        Returns
        -------
        None
        '''
        orient_history = np.array(self._FreeFallSimulation.orientHistory)
        orient_history.shape = (int(orient_history.shape[0]/4),4)
        
        return orient_history
    
    def get_sim_step(self):
        sim_step = self._FreeFallSimulation.getSimStep()
        return sim_step
        
    def evolve_motion_inertial(self, nstep=None):
        '''
        Evolve the body motion with no external forces.
        
        Parameters
        ----------
        nstep : int, optional
            The number of simulation steps.
        
        Returns
        -------
        None
        '''
        if nstep is None:
            nstep = self.nt
        self._FreeFallSimulation.evolveMotionInertial(nstep)
        
        # update python particle object with wrapped pybind11 particle object
        self.st._State = self._FreeFallSimulation.st
        return
        
    def evolve_motion_buoyancy(self, nstep=None):
        '''
        Evolve the body motion with an external buoyancy gradient.
        
        Parameters
        ----------
        nstep : int, optional
            The number of simulation steps.
        
        Returns
        -------
        None
        '''
        if nstep is None:
            nstep = self.nt
        self._FreeFallSimulation.evolveMotionBuoyancy(nstep)
        
        # update python particle object with wrapped pybind11 particle object
        self.par._State = self._FreeFallSimulation.st
        return
    
