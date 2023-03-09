import numpy as np
import matplotlib
matplotlib.use("TkCairo")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.animation import PillowWriter
from matplotlib import colors

# read data
x1 = np.genfromtxt('x1.txt')
y1 = np.genfromtxt('y1.txt')
x2 = np.genfromtxt('x2.txt')
y2 = np.genfromtxt('y2.txt')

data1 = np.genfromtxt('rigid1.txt')
data2 = np.genfromtxt('rigid2.txt')

tf = 2
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

#fx = np.concatenate((fx1.T,fx2.T)).T
#fy = np.concatenate((fy1.T,fy2.T)).T
#fmag = np.sqrt(fx**2.+fy**2.)
#print(np.min(fmag), np.max(fmag))

# plot
fig = plt.figure(figsize=(7,7))
plts = []
dp = 5.

mycmap = colors.ListedColormap(['c', 'm'])

ax = fig.add_subplot(1,1,1)
for i in range(nt):
    # transform points
    x[bind==0] = x1*np.cos(omeg[0,i])-y1*np.sin(omeg[0,i])+xcom[0,i]
    y[bind==0] = x1*np.sin(omeg[0,i])+y1*np.cos(omeg[0,i])+ycom[0,i]
    x[bind==1] = x2*np.cos(omeg[1,i])-y2*np.sin(omeg[1,i])+xcom[1,i]
    y[bind==1] = x2*np.sin(omeg[1,i])+y2*np.cos(omeg[1,i])+ycom[1,i]

    #pl1, = ax.plot([px1[i],px2[i]], [py1[i],py2[i]], 'k-')
    pl2 = ax.scatter(x, y, c=bind, s=10, cmap=mycmap)
    #pl2 = ax.scatter(x[i,:], y[i,:], c=fmag[i,:], s=10, cmap='Spectral_r', vmin=0., vmax=10.)
    #plt.scatter(x2, y2, c='c', s=10)
    ax.set_aspect(1.)
    ax.set_xlim([-dp,dp])
    ax.set_ylim([-dp,dp])

    #plts.append([pl1,pl2])
    plts.append([pl2])
    
# save animation
ani = animation.ArtistAnimation(fig, plts, interval=1, blit=False, repeat_delay=500)

writer = PillowWriter(fps=30)
ani.save('rigid.gif', writer=writer, dpi=90)
