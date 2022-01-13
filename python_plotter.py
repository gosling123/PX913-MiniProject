import netCDF4 as NC
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
# matplotlib.rcParams['text.usetex'] = True
from matplotlib import cm
import matplotlib.patches as mpatches
import os


## Plotting Parameters 
from pylab import rcParams
rcParams['figure.figsize'] = 16, 7
rcParams['lines.linewidth'] = 2
rcParams['font.size'] = 14
rcParams['axes.titlesize'] = 20
rcParams['axes.labelsize'] = 16
rcParams['xtick.labelsize'] = 14
rcParams['ytick.labelsize'] =  14


### CHECK FILE EXISTS
if not os.path.isfile("particle_pusher_data.nc"):
    print("Python_plotter: No Data File Found")
    quit()



### READ NETCDF FILE
dat = NC.Dataset("particle_pusher_data.nc", "r", format="NETCDF4")

### ADD VARIABLES TO ARRAYS
grid_rho = dat.variables['rho']
grid_phi = dat.variables['phi']
grid_E = dat.variables['E']


grid_r = dat.variables['r']
grid_v = dat.variables['v']
grid_a = dat.variables['a']


### SIMULATION SETUP
ver_iter = dat.Verlet_Iteration
x_spacing = dat.Delta_x
y_spacing = dat.Delta_y
problem = dat.Prob


# GRIDS FOR COLOURMAP
x_grid = np.arange(-1, 1, x_spacing) 
y_grid = np.arange(-1, 1, y_spacing)


rho = grid_rho
phi = grid_phi
E = grid_E

print(grid_r.shape)
print(ver_iter)

### DEFINED UP TO NUMBER OF VERLET ITERATIONS AS IT STOPS ONCE OUT OF THE GRID
r = grid_r[:, :ver_iter]

v = grid_v[:, :ver_iter]

a = grid_a[:, :ver_iter]


rho = np.array(rho).T
phi = np.array(phi).T
E = np.array(E)
E = np.swapaxes(E, 1, 2)
r = np.array(r)
v = np.array(v)
a = np.array(a)



### PLOTS

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

fig.suptitle(r"Values for rho, phi, and electric field for problem '" + problem + "'")

y, x = np.meshgrid(y_grid, x_grid)

ax1.set_xlabel(r"$X$")
ax1.set_ylabel(r"$Y$")
pos1 = ax1.pcolor(x, y, rho)
fig.colorbar(pos1, ax=ax1, label=r'$\rho$')
ax1.plot(r[0,:], r[1, :], "r--", label=r"particle trajectory")
ax1.set_xlim(-1,1)
ax1.set_ylim(-1,1)



ax2.set_xlabel(r"$X$")
ax2.set_ylabel(r"$Y$")
pos2 = ax2.pcolor(x, y, phi)
fig.colorbar(pos2, ax=ax2, label=r'$\phi$')
ax2.set_xlim(-1,1)
ax2.set_ylim(-1,1)

ax3.set_xlabel(r"$X$")
ax3.set_ylabel(r"$Y$")
pos3 = ax3.pcolor(x, y, E[0,:,:])
fig.colorbar(pos3, ax=ax3, label=r'$E_x$')
ax3.set_xlim(-1,1)
ax3.set_ylim(-1,1)

ax4.set_xlabel(r"$X$")
ax4.set_ylabel(r"$Y$")
pos4 = ax4.pcolor(x, y, E[1,:,:])
fig.colorbar(pos4, ax=ax4, label=r'$E_y$')
ax4.set_xlim(-1,1)
ax4.set_ylim(-1,1)

ax1.legend()

plt.show()



plt.figure()
plt.scatter(r[0, :], r[1, :])
plt.xlabel(r'X')
plt.ylabel(r'Y')
plt.title(r'Scatter Plot For X vs Y')


plt.show()