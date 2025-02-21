from rigidpy.pybind11_lib._mesh_particle import MeshParticle
import numpy as np
import pyvista as pv

# create mesh particle
def create_mesh_particle():
    # read mesh file
    mesh = pv.read('data/crystal_remesh.stl')
    print(mesh)
    vertices = mesh.points
    faces = mesh.regular_faces
    density = 920.
    
    # create mesh particle object
    mp = MeshParticle(vertices, faces, density)

    return mp, mesh

# test area
def test_area():
    mp, mesh = create_mesh_particle()
    mp.calculateFaceAreas()
    areas = mp.faceAreas
    
    assert np.abs(np.sum(areas)-mesh.area)<1.e-16
    
