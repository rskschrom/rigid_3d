from .pybind11_lib._mesh_particle import MeshParticle as _MeshParticle
from .pybind11_lib._quat import quat
import numpy as np
import pyvista as pv

class Particle():
    '''
    Class for particles defined by a closed, triangulated, surface mesh.
    
    Parameters
    ----------
    vertices : ndarray (N,3) of float
        The coordinates of the N vertices of the mesh.
    faces : ndarray (M,3) of int
        The indices of the M mesh faces. Each index corresponds to a point within the `vertices` array.
    rho : float
        The density of the particle.
    '''
    def __init__(self, vertices, faces, rho):
        self._MeshParticle = _MeshParticle(vertices, faces, rho)
        
    @classmethod
    def from_file(cls, filepath, rho):
        '''
        Read the mesh from a `.stl` file and initialize the object.
        
        Parameters
        ----------
        filepath : str
            The path to the mesh `.stl` file.
        rho : float
            The density of the particle.
            
        Returns
        -------
        The `Particle` object.
        '''
        # read mesh file
        mesh = pv.read(filepath)
        vertices = mesh.points
        faces = mesh.regular_faces
        
        return cls(vertices, faces, rho)
    
    def get_mat_inerm(self):
        return self._MeshParticle.getMatInerm()

    def get_vertices(self):
        return self._MeshParticle.getVertices()

        
