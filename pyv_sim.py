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
r1 = np.genfromtxt('r1.txt')
orient1 = np.genfromtxt('orient1.txt')
pos1 = np.genfromtxt('pos1.txt')

r2 = np.genfromtxt('r2.txt')
orient2 = np.genfromtxt('orient2.txt')
pos2 = np.genfromtxt('pos2.txt')

x1 = r1[::3]
y1 = r1[1::3]
z1 = r1[2::3]

x2 = r2[::3]
y2 = r2[1::3]
z2 = r2[2::3]

ti = -1
pos1 = pos1[ti,:]
pos2 = pos2[ti,:]
orient1 = orient1[ti,:]
orient2 = orient2[ti,:]
nt = pos1.shape[0]
print(orient1.shape, x1.shape)

points1 = np.vstack((x1,y1,z1)).T
points2 = np.vstack((x2,y2,z2)).T

pointst1 = transform(points1, orient1, pos1)
pointst2 = transform(points2, orient2, pos2)
pointst = np.concatenate((pointst1,pointst2))

pc = pv.PolyData(pointst)

print(pointst.shape)

# plot
dp = 4.
pl = pv.Plotter()
pl.window_size = [800,800]

pl.enable_eye_dome_lighting()
pl.add_mesh(pc, render_points_as_spheres=True, point_size=5,)
pl.add_axes()

pl.set_background('k')

pl.show()
