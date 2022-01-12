import netCDF4 as NC
import matplotlib.pyplot as plt

#Read the netcdf file into variable dat
dat = NC.Dataset("particle_pusher_data.nc","r",format = "NETCDF4")

#Get x-component of electric field
#and x and y particle positions
Ex = dat.variables['Ex'][:]
x = dat.variables['x'][:]
y = dat.variables['y'][:]

#Subplots with 1 row 2 columns
fig, axes = plt.subplots(1,2)

#pseudocolour plot of Ex
pcolor=axes[0].pcolor(Ex)
axes[0].set_title('Pseudocolor plot of Ex')
axes[0].set_xlabel('X')
axes[0].set_ylabel('Y')
fig.colorbar(pcolor,ax=axes[0],label='Ex')


#Scatter plot of particle position x vs particle position y
axes[1].scatter(x,y)
axes[1].set_title('Scatter Plot of x vs y')
axes[1].set_xlabel('X')
axes[1].set_ylabel('Y')

#Display the plot
plt.show()                                                                                                                                                                                    
