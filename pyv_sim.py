import pyvista as pv
import numpy as np

# read data
sf = 1
x1 = np.genfromtxt('x1.txt')[::sf]
y1 = np.genfromtxt('y1.txt')[::sf]
x2 = np.genfromtxt('x2.txt')[::sf]
y2 = np.genfromtxt('y2.txt')[::sf]

data1 = np.genfromtxt('rigid1.txt')
data2 = np.genfromtxt('rigid2.txt')

tf = 5
xcom1 = data1[::tf,0]
ycom1 = data1[::tf,1]
omeg1 = data1[::tf,2]
xcom2 = data2[::tf,0]
ycom2 = data2[::tf,1]
omeg2 = data2[::tf,2]

xcom = np.vstack((xcom1,xcom2))
ycom = np.vstack((ycom1,ycom2))
omeg = np.vstack((omeg1,omeg2))

nt = xcom1.shape[0]

# reshape arrays for animation plot
np1 = len(x1)
np2 = len(x2)
bind = np.zeros([np1+np2])
bind[np1:] = 1.

x = np.concatenate((x1.T,x2.T)).T
y = np.concatenate((y1.T,y2.T)).T

# plot
dp = 5.
pl = pv.Plotter(notebook=False, off_screen=True)
pl.open_gif('test.gif', fps=30)

chart = pv.Chart2D()
chart.grid = False
chart.background_color = (1., 1., 1.)
chart.x_range = [-dp, dp]
chart.y_range = [-dp, dp]
chart.hide_axes()

plot = chart.scatter(x, y, color='m', size=5)
pl.add_chart(chart)

for i in range(nt):
    # transform points
    x[bind==0] = x1*np.cos(omeg[0,i])-y1*np.sin(omeg[0,i])+xcom[0,i]
    y[bind==0] = x1*np.sin(omeg[0,i])+y1*np.cos(omeg[0,i])+ycom[0,i]
    x[bind==1] = x2*np.cos(omeg[1,i])-y2*np.sin(omeg[1,i])+xcom[1,i]
    y[bind==1] = x2*np.sin(omeg[1,i])+y2*np.cos(omeg[1,i])+ycom[1,i]
 
    plot.update(x, y)
    pl.write_frame()
pl.close()
