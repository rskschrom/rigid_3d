from rigidpy.pybind11_lib._mesh_particle import MeshParticle
import numpy as np
import pyvista as pv

# create mesh particle
def create_mesh_particle():
    # read mesh file
    mesh = pv.read('data/crystal_remesh.stl')
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
    
# test inertia moment tensor
def test_inertia():
    #mp, mesh = create_mesh_particle()
    
    # create cuboid
    dx = 0.25
    dy = 0.06
    dz = 0.4
    
    mesh = pv.Cube(x_length=dx, y_length=dy, z_length=dz)
    mesh = mesh.triangulate()
    
    vertices = mesh.points
    faces = mesh.regular_faces
    density = 920.
        
    # create mesh particle object
    mp = MeshParticle(vertices, faces, density)
    mp.calculateFaceAreasNorms()
    mass = mp.totalMass()
    
    # evaluate diagonal tensor compared to known value
    itens = mp.inertiaMomentTensor()
    
    itens_true = np.array([mass/12.*(dy**2.+dz**2.),
                           mass/12.*(dx**2.+dz**2.),
                           mass/12.*(dx**2.+dy**2.)])
    print(itens_true)
    err = np.sqrt((itens[0,0]-itens_true[0])**2.+
                  (itens[1,1]-itens_true[1])**2.+
                  (itens[2,2]-itens_true[2])**2.)
    assert err<1.e-8

# test integral over individual mesh faces
def test_face_integral():
    # create cuboid
    dx = 0.25
    dy = 0.06
    dz = 0.4
    
    mesh = pv.Cube(x_length=dx, y_length=dy, z_length=dz)
    mesh = mesh.triangulate()
    
    vertices = mesh.points
    faces = mesh.regular_faces
    nface = faces.shape[0]
    density = 0.92
        
    # create mesh particle object
    mp = MeshParticle(vertices, faces, density)
    mp.calculateFaceAreasNorms()

    # pick a face integral with constant x
    fi = 2
    verts = vertices[faces[fi,:],:]
    area = mp.faceAreas[fi]
    
    x = verts[:,0]
    ixx = mp.triAlp2BetIntegral(x, x, area)
    assert (ixx-area*x[0]**3.)/(area*x[0]**3.)*100.<1.e-3