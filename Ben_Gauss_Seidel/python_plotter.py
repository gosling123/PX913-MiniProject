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



### loop to read all the files and put all useful data into the lits above

dat_x = NC.Dataset("E_x.nc", "r", format="NETCDF4")
dat_y = NC.Dataset("E_y.nc", "r", format="NETCDF4")
grid_Ex = dat_x.variables['grid_data']
grid_Ey = dat_y.variables['grid_data']
Ex = grid_Ex[:]
Ey = grid_Ey[:]

dat_phi = NC.Dataset("phi.nc", "r", format="NETCDF4")
dat_rho = NC.Dataset("rho.nc", "r", format="NETCDF4")
grid_phi = dat_phi.variables['grid_data']
grid_rho = dat_rho.variables['grid_data']
rho = grid_rho[:]
phi = grid_phi[:]

phi = np.array(phi)


norm_phi = np.linalg.norm(phi)
print(norm_phi)

##create colours to map onto -1 and +1 spins
# cmap = cm.copper
# values = [-1, +1] ### spins 
# i=0 ### to set first frame of animation

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)


pos1 = ax1.imshow(Ex)
ax1.set_title("Elec X")
fig.colorbar(pos1, ax=ax1)

ax2.set_title("Elec Y")
pos2 = ax2.imshow(Ey)
fig.colorbar(pos2, ax=ax2)


pos3 = ax3.imshow(rho)
ax3.set_title("rho")
fig.colorbar(pos3, ax=ax3)

ax4.set_title("phi")
# ax4.imshow(np.log(np.abs(phi)))
pos4 = ax4.imshow(phi)
fig.colorbar(pos4, ax=ax4)



# ###create legend
# colors = [ im.cmap(im.norm(value)) for value in values]
# patches = [mpatches.Patch(color=colors[k], label="Spin {l}".format(l=values[k]) ) for k in range(len(values)) ]
# # put those patched as legend-handles into the legend
# plt.legend(handles=patches, bbox_to_anchor=(0.2, 1.1), ncol=2 , borderaxespad=0,  fontsize=16)




plt.show()
