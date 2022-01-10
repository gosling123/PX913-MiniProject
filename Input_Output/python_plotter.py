import netCDF4 as NC
import matplotlib.pyplot as plt

#Read the netcdf file into variable dat
dat = NC.Dataset("particle_pusher_data.nc","r",format = "NETCDF4")

#Get the 2D-grid and time history for plotting
Ex = dat.variables['Ex'][:]
x = dat.variables['x'][:]
y = dat.variables['y'][:]

#Subplot with 1 row 2 columns
#Multi panel plot as asked in extension
fig, axes = plt.subplots(1,2)

#Using imshow for displaying grid
axes[0].pcolor(Ex)
axes[0].set_title('Ex')
axes[0].set_xlabel('X')
axes[0].set_ylabel('Y')
#Creating custom legend

#Using plot for plotting time history
axes[1].scatter(x,y)
axes[1].set_title('Scatter Plot')
axes[1].set_xlabel('xlabel')
axes[1].set_ylabel('ylabel')

#Display the plot
plt.show()                                                                                                                                                                                    
