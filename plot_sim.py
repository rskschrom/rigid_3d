import numpy as np
import matplotlib
matplotlib.use("TkCairo")
import matplotlib.pyplot as plt

# read data
x1 = np.genfromtxt('x1.txt')
y1 = np.genfromtxt('y1.txt')
x2 = np.genfromtxt('x2.txt')
y2 = np.genfromtxt('y2.txt')

hx1_1d = np.genfromtxt('hx1.txt')
hy1_1d = np.genfromtxt('hy1.txt')
hx2_1d = np.genfromtxt('hx2.txt')
hy2_1d = np.genfromtxt('hy2.txt')

print(hx1_1d.shape)
hx1, hy1 = np.meshgrid(hx1_1d, hy1_1d, indexing='ij')
hx2, hy2 = np.meshgrid(hx2_1d, hy2_1d, indexing='ij')

h1 = np.genfromtxt('hash1.txt')
h2 = np.genfromtxt('hash2.txt')

data1 = np.genfromtxt('rigid1.txt')
data2 = np.genfromtxt('rigid2.txt')

ti = -1
xcom1 = data1[ti,0]
ycom1 = data1[ti,1]
omeg1 = data1[ti,2]
xcom2 = data2[ti,0]
ycom2 = data2[ti,1]
omeg2 = data2[ti,2]

# transform points
x1t = x1*np.cos(omeg1)-y1*np.sin(omeg1)+xcom1
y1t = x1*np.sin(omeg1)+y1*np.cos(omeg1)+ycom1
x2t = x2*np.cos(omeg2)-y2*np.sin(omeg2)+xcom2
y2t = x2*np.sin(omeg2)+y2*np.cos(omeg2)+ycom2

hx1t = hx1*np.cos(omeg1)-hy1*np.sin(omeg1)+xcom1
hy1t = hx1*np.sin(omeg1)+hy1*np.cos(omeg1)+ycom1
hx2t = hx2*np.cos(omeg2)-hy2*np.sin(omeg2)+xcom2
hy2t = hx2*np.sin(omeg2)+hy2*np.cos(omeg2)+ycom2

x = np.concatenate((x1t.T,x2t.T)).T
y = np.concatenate((y1t.T,y2t.T)).T
h = np.concatenate((h1.T,h2.T)).T
hx = np.concatenate((hx1t.T,hx2t.T)).T
hy = np.concatenate((hy1t.T,hy2t.T)).T

#fx = np.concatenate((fx1.T,fx2.T)).T
#fy = np.concatenate((fy1.T,fy2.T)).T
#fmag = np.sqrt(fx**2.+fy**2.)
#print(np.min(fmag), np.max(fmag))

# plot
fig = plt.figure(figsize=(7,7))
dp = 5.

ax = fig.add_subplot(1,1,1)

plt.scatter(hx, hy, c='k', s=20.)
plt.scatter(x, y, c=h, s=1, cmap='Spectral_r')
plt.colorbar()
ax.set_aspect(1.)
ax.set_xlim([-dp,dp])
ax.set_ylim([-dp,dp])
plt.savefig('sim.png', dpi=150)
