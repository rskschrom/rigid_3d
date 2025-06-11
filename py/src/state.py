from .pybind11_lib._state import State as _State
import numpy as np

class State():
    '''
    Class for the body state.
    
    Parameters
    ----------
    mat_inerm : ndarray (3,3) of float
        The inertia momentum tensor of the body in the lab-relative reference frame.
    com_pos : tuple (3) of float, default=(0,0,0)
        The body center of mass position.
    com_vel : tuple (3) of float, default=(0,0,0)
        The body velocity.
    orient : tuple (4) of float, default=(1,0,0,0)
        The orientation quaternion of the body in the lab-relative reference frame.
    omega : tuple (3) of float, default=(0,0,0)
        The angular velocity of the body.
        
    Returns
    -------
    The `State` object.
    
    '''
    def __init__(self, mat_inerm, com_pos=(0.,0.,0.), com_vel=(0.,0.,0.),
                 orient=(1.,0.,0.,0.), omega=(0.,0.,0.)):
        self._State = _State(mat_inerm, com_pos, com_vel, orient, omega)     

    def get_mat_inerm(self):
        return self._State.getMatInerm()
        
    def get_com_pos(self):
        return self._State.getComPos()
        
    def get_com_vel(self):
        return self._State.getComVel()
        
    def get_orient(self):
        return self._State.getOrient()
    
    def get_omega(self):
        return self._State.getOmega()
    
    def set_mat_inerm(self, mat_inerm):
        self._State.setMatInerm(mat_inerm)
        
    def set_com_pos(self, com_pos):
        self._State.setComPos(com_pos)
        return
        
    def set_com_vel(self, com_vel):
        self._State.setComVel(com_vel)
        return
        
    def set_orient(self, orient):
        self._State.setOrient(orient)
        return
    
    def set_omega(self, omega):
        self._State.setOmega(omega)
        return
        
    def set_orient_zenith(self, theta):
        '''
        Helper function to set the body orientation according to a zenith angle.
        
        Parameters
        ----------
        theta : float
            The zenith angle in degrees.
            
        Returns
        -------
        None
        '''
        orient = [np.cos(theta*np.pi/180./2.),0.,
                  np.sin(theta*np.pi/180./2.),0.]
        self.set_orient(orient)
        return
        
    def rotational_potential_energy(self, g, rhof_grad, rhob):
        rotPE = self._State.rotationalPotentialEnergy(g, rhof_grad, rhob)
        return rotPE
    
    def rotational_kinetic_energy(self):
        rotKE = self._State.rotationalKineticEnergy()
        return rotKE
