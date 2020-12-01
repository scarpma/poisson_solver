import numpy as np
import matplotlib.pyplot as plt
import sys

N = int(sys.argv[1])
nx, ny = N, N
D = 0.5
L = 0.3
boxC = 0.5
# load flattened data with fortran order (y1x1 y1x2 y1x3 ... y2x1 y2x2 y2x3 ...)
# (x is contiguous in memory)
data = np.loadtxt('sol.dat')
# reshape and transpose to convert to "standard" c order
# (y is contiguous in memory)
data = data.reshape((ny,nx)).T
x = np.linspace(0, 1, nx)
y = np.linspace(0, 1, ny)
X, Y = np.meshgrid(x, y, indexing='ij')
Ey, Ex = np.gradient(-data, y, x)

# create figure and axis
fig, ax = plt.subplots()
ax.set_xlabel('y')
ax.set_ylabel('x')

# plot 2d solution with contour lines
im = ax.imshow(data, extent=[x[0],x[-1],y[0],y[-1]], origin='lower', interpolation='nearest')
fig.colorbar(im, ax=ax)

# plot streamlines of the gradient field (electric field)
Em = np.sqrt(Ey**2. + Ex**2.)
lw = 8. * Em / Em.max() # linewidth depending on magnitude
ax.streamplot(x, y, Ex, Ey, color='white', linewidth=lw, arrowsize=0.7, density=1.2)

cntr = ax.contour(data, [-0.7, -0.25, -0.05, 0.05, 0.25, 0.7], colors='red', extent=[0,1,0,1])
ax.clabel(cntr, cntr.levels, fontsize=10, colors='red') # plot contour values

# plot bars inside domain
ax.vlines(x=boxC-D/2., ymin=boxC-L/2., ymax=boxC+L/2., linewidth=2, color='b')
ax.vlines(x=boxC+D/2., ymin=boxC-L/2., ymax=boxC+L/2., linewidth=2, color='b')

plt.show()
