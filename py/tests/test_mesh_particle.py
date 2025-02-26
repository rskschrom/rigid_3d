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

# test area and normal vectors
def test_area_norms():
    mp, mesh = create_mesh_particle()
    mp.calculateFaceAreasNorms()
    areas = mp.faceAreas
    norms = mp.faceNorms
        
    # compute normals with pyvista
    mesh.compute_normals(inplace=True)
    norm_pv = mesh['Normals']
    
    # test dot product
    dotp = np.einsum('ki,ki->k', norms, norm_pv)
    
    assert np.abs(np.sum(dotp)-float(len(dotp)))<1.e-16
    assert np.abs(np.sum(areas)-mesh.area)<1.e-16

# test mass
def test_mass():
    mp, mesh = create_mesh_particle()
    mp.calculateFaceAreasNorms()
    mass = mp.totalMass()

    assert np.abs(mass-mesh.volume*920.)/(mesh.volume*920.)<1.e-3
        
# test triangle integral
def test_triangle_integral():
    mp, mesh = create_mesh_particle()
    
    # test with right triangle
    x1 = 0.5
    y1 = 0.3
    area = 0.5*x1*y1
    alp = np.array([0.,x1,0.])
    bet = np.array([0.,0.,y1])
    integral = mp.triAlp2BetIntegral(alp, bet, area)
    
    assert np.abs(x1**3.*y1**2./60.-integral)<1.e-7