import pyvista as pv
import numpy as np

# quaternion multiplication
def qmult(q1, q2):
    r1 = q1[:,0]
    v1x = q1[:,1]
    v1y = q1[:,2]
    v1z = q1[:,3]
  
    r2 = q2[:,0]
    v2x = q2[:,1]
    v2y = q2[:,2]
    v2z = q2[:,3]

    r3 = r1*r2-v1x*v2x-v1y*v2y-v1z*v2z;
    v3x = r1*v2x+r2*v1x+v1y*v2z-v1z*v2y;
    v3y = r1*v2y+r2*v1y-v1x*v2z+v1z*v2x;
    v3z = r1*v2z+r2*v1z+v1x*v2y-v1y*v2x;

    q3 = np.vstack((r3,v3x,v3y,v3z)).T
    return q3
    
# quaternion rotation and translation
def transform(points, orient, com):    
    qpoints = np.hstack((np.zeros([points.shape[0],1]),points))
    orient.shape = (1,4)
    q1 = qmult(orient, qpoints)
    orient[:,1:] = -orient[:,1:]
    q2 = qmult(q1, orient) 
    
    points_tr = q2[:,1:]
    points_tr[:,0] = points_tr[:,0]+com[0]
    points_tr[:,1] = points_tr[:,1]+com[1]
    points_tr[:,2] = points_tr[:,2]+com[2]
    return points_tr

# read data
r = np.genfromtxt('r.txt')
orient = np.genfromtxt('orientation.txt')
pos = np.genfromtxt('position.txt')
orient.shape = (int(orient.shape[0]/4),4)
pos.shape = (int(pos.shape[0]/3),3)

x = r[::3]
y = r[1::3]
z = r[2::3]

print(np.min(x), np.max(x))
print(np.min(y), np.max(y))
print(np.min(z), np.max(z))

tf = 50
pos = pos[::tf,:]
orient = orient[::tf,:]
nt = pos.shape[0]
print(orient.shape, x.shape)

points = np.vstack((x,y,z)).T
pc = pv.PolyData(points)

print(points.shape)

# plot
dp = 4.
pl = pv.Plotter(notebook=False, off_screen=True)
pl.open_gif('test.gif', fps=20)
pl.window_size = [800,800]

pl.enable_eye_dome_lighting()
pl.add_mesh(pc, render_points_as_spheres=True, point_size=10, color='m')
pl.add_axes()

pl.set_background('w')

# camera settings
pl.camera.position = (0.01, 0., 0.)
pl.camera.focal_point = (0., 0., 0.)

for i in range(nt):
    # transform points
    #pointst = transform(points, orient[i,:], pos[i,:])
    pointst = transform(points, orient[i,:], [0.,0.,0.])
    pl.update_coordinates(pointst, render=False)
    pl.write_frame()

pl.close()
