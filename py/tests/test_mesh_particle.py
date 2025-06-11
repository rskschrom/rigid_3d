from rigidpy.particle import Particle
import numpy as np
import pyvista as pv

# create mesh particle
def create_mesh_particle():
    # get mesh as well
    mesh = pv.read('data/crystal_remesh.stl')

    # create mesh particle object
    par = Particle.from_file('data/crystal_remesh.stl', 920.)

    return par, mesh

# test area and normal vectors
def test_area_norms():
    par, mesh = create_mesh_particle()
    par._MeshParticle.calculateFaceAreasNorms()
    areas = par._MeshParticle.faceAreas
    norms = par._MeshParticle.faceNorms
        
    # compute normals with pyvista
    mesh.compute_normals(inplace=True)
    norm_pv = mesh['Normals']
    
    # test dot product
    dotp = np.einsum('ki,ki->k', norms, norm_pv)
    
    assert np.abs(np.sum(dotp)-float(len(dotp)))<1.e-16
    assert np.abs(np.sum(areas)-mesh.area)<1.e-16

# test mass
def test_mass():
    par, mesh = create_mesh_particle()
    par._MeshParticle.calculateFaceAreasNorms()
    mass = par._MeshParticle.totalMass()

    assert np.abs(mass-mesh.volume*920.)/(mesh.volume*920.)<1.e-3
        
# test triangle integral
def test_triangle_integral():
    par, mesh = create_mesh_particle()
    
    # test with right triangle
    x1 = 0.5
    y1 = 0.3
    area = 0.5*x1*y1
    alp = np.array([0.,x1,0.])
    bet = np.array([0.,0.,y1])
    integral = par._MeshParticle.triAlp2BetIntegral(alp, bet, area)
    assert np.abs(x1**3.*y1**2./60.-integral)<1.e-7
    
# test inertia moment tensor
def test_inertia():
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
    par = Particle(vertices, faces, density)
    par._MeshParticle.calculateFaceAreasNorms()
    mass = par._MeshParticle.totalMass()
    
    # evaluate diagonal tensor compared to known value
    itens = par._MeshParticle.inertiaMomentTensor()
    
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
    par = Particle(vertices, faces, density)
    par._MeshParticle.calculateFaceAreasNorms()

    # pick a face integral with constant x
    fi = 2
    verts = vertices[faces[fi,:],:]
    area = par._MeshParticle.faceAreas[fi]
    
    x = verts[:,0]
    ixx = par._MeshParticle.triAlp2BetIntegral(x, x, area)
    assert (ixx-area*x[0]**3.)/(area*x[0]**3.)*100.<1.e-3