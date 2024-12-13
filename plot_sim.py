import numpy as np
import matplotlib
matplotlib.use("TkCairo")
import matplotlib.pyplot as plt

# read data
rt1 = np.genfromtxt('rt1.txt')
xt1 = rt1[::3]
yt1 = rt1[1::3]
zt1 = rt1[2::3]

# plot
fig = plt.figure(figsize=(7,7))
dp = 3.

ax = fig.add_subplot(1,1,1)

plt.scatter(-yt1, zt1, c='m', s=1, alpha=0.2)
ax.set_aspect(1.)
ax.set_xlim([-dp,dp])
ax.set_ylim([-dp,dp])
plt.savefig('sim.png', dpi=150)
