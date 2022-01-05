import netCDF4 as NC
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import numpy as np
matplotlib.rcParams['text.usetex'] = True
from matplotlib import cm
import matplotlib.patches as mpatches


## Plotting Parameters 
from pylab import rcParams
rcParams['figure.figsize'] = 18, 9
rcParams['lines.linewidth'] = 3
rcParams['font.size'] = 16
rcParams['axes.titlesize'] = 24
rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 16
rcParams['ytick.labelsize'] =  16



### READ NETCDF FILE

dat = NC.Dataset("particle_pusher_data.nc", "r", format="NETCDF4")

# print(dat.__dict__)

### ADD VARIABLES TO ARRAYS
grid_rho = dat.variables['rho']
grid_phi = dat.variables['phi']
grid_E_x = dat.variables['E_x']
grid_E_y = dat.variables['E_y']
grid_x = dat.variables['x']
grid_y = dat.variables['y']
grid_v_x = dat.variables['v_x']
grid_v_y = dat.variables['v_y']
grid_a_x = dat.variables['a_x']
grid_a_y = dat.variables['a_y']


### SIMULATION SETUP
ver_iter = dat.Verlet_Iteration
timestep = dat.Time_step
x_spacing = dat.Delta_x
y_spacing = dat.Delta_y
problem = dat.Prob
n_x = dat.n_x
n_y = dat.n_y

### AXIS FOR PLOTS
time = np.arange(0, ver_iter)*timestep

### GRIDS FOR COLORMAP
x_grid_rp = np.arange(-1-x_spacing, 1+x_spacing, x_spacing)
y_grid_rp = np.arange(-1-y_spacing, 1+y_spacing , y_spacing)

x_grid_E = np.arange(-1, 1, x_spacing)
y_grid_E = np.arange(-1, 1, y_spacing)


rho = grid_rho[:]
phi = grid_phi[:]
E_x = grid_E_x[:]
E_y = grid_E_y[:]


### DEFINED UP TO NUMBER OF VERLET ITERATIONS AS IT STOPS ONCE OUT OF THE GRID
x = grid_x[:ver_iter]
y = grid_y[:ver_iter]

v_x = grid_v_x[:ver_iter]
v_y = grid_v_y[:ver_iter]

a_x = grid_a_x[:ver_iter]
a_y = grid_a_y[:ver_iter]


rho = np.array(rho)
phi = np.array(phi)
E_x = np.array(E_x)
E_y = np.array(E_y)
x = np.array(x)
y = np.array(y)
v_x = np.array(v_x)
v_y = np.array(v_y)
a_x = np.array(a_x)
a_y = np.array(a_y)



#### PLOTS

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

y_rp, x_rp = np.meshgrid(y_grid_rp, x_grid_rp)
y_E, x_E = np.meshgrid(y_grid_E, x_grid_E)

ax1.set_xlabel(r"$X$")
ax1.set_ylabel(r"$Y$")
pos1 = ax1.pcolor(x_rp, y_rp, rho)
fig.colorbar(pos1, ax=ax1, label=r'$\rho$')

ax2.set_xlabel(r"$X$")
ax2.set_ylabel(r"$Y$")
pos2 = ax2.pcolor(x_rp, y_rp, phi)
fig.colorbar(pos2, ax=ax2, label=r'$\phi$')

ax3.set_xlabel(r"$X$")
ax3.set_ylabel(r"$Y$")
pos3 = ax3.pcolor(x_E, y_E, E_x)
fig.colorbar(pos3, ax=ax3, label=r'$E_x$')

ax4.set_xlabel(r"$X$")
ax4.set_ylabel(r"$Y$")
pos4 = ax4.pcolor(x_E, y_E, E_y)
fig.colorbar(pos4, ax=ax4, label=r'$E_y$')

plt.show()

fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(3, 2)

pos1 = ax1.plot(time, x)
pos2 = ax2.plot(time, y)
pos3 = ax3.plot(time, v_x)
pos4 = ax4.plot(time, v_y)
pos5 = ax5.plot(time, a_x)
pos6 = ax6.plot(time, a_y)

ax5.set_xlabel(r'Time')
ax6.set_xlabel(r'Time')


ax1.set_ylabel(r'$X$')
ax2.set_ylabel(r'$Y$')
ax3.set_ylabel(r'$V_x$')
ax4.set_ylabel(r'$V_y$')
ax5.set_ylabel(r'$a_x$')
ax6.set_ylabel(r'$a_y$')


plt.show()