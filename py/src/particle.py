from .pybind11_lib._particle import Particle as _Particle
from .pybind11_lib._quat import quat
import numpy as np
import os

class Particle():
    '''
    Particle class
    '''
    def __init__(self, point_file, point_masses):
        self.point_file = point_file
        self.point_mass = point_masses
        
        try:
            self._Particle = _Particle(point_file, point_masses)
        except:
            print("File not found.")       

    def set_com_pos(self, pos):
        self._Particle.setComPos(pos)
        return
        
    def set_com_vel(self, vel):
        self._Particle.setComVel(vel)
        return
        
    def set_orient(self, orient):
        self._Particle.setOrient(orient)
        return
        
    def set_orient_zenith(self, theta):
        orient = [np.cos(theta*np.pi/180./2.),0.,
                  np.sin(theta*np.pi/180./2.),0.]
        self.set_orient(orient)
        return
        
    def set_omega(self, omega):
        self._Particle.setOmega(omega)
        return
        
    def set_omega_body(self, omega):
        omega = quat.vecRotate(omega, self._Particle.orient)
        self._Particle.setOmega(omega)
        return
        
    def get_rel_points(self):
        rel_points = np.array(self._Particle.relPoints)
        rel_points.shape = (int(rel_points.shape[0]/3),3)
        return rel_points
    
    def rotational_potential_energy(self, g, rhof_grad, rhob):
        rotPE = self._Particle.rotationalPotentialEnergy(g, rhof_grad, rhob)
        return rotPE
    
    def rotational_kinetic_energy(self):
        rotKE = self._Particle.rotationalKineticEnergy()
        return rotKE
    #writeVector(par.relPoints, "r.txt");
