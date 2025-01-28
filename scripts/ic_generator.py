import f90nml
import numpy as np
import xarray as xr

config = f90nml.read('config.nml')

#=======================================================
# Domain information

nx = 200
ny = 200
Lx = 1500e3
Ly = 1500e3

dx = Lx / nx
dy = Ly / ny

mass_xs, mass_ys = np.meshgrid(
    np.linspace((dx-Lx)/2, (Lx-dx)/2, nx),
    np.linspace((dy-Ly)/2, (Ly-dy)/2, ny))
face_xs = mass_xs - dx/2
face_ys = mass_ys - dy/2

#=======================================================
# Initialise data arrays

h = np.zeros((nx,ny))
u = np.zeros((nx,ny))
v = np.zeros((nx,ny))

#=======================================================
# Set data arrays

# This creates a balanced flow around a low pressure region at (x_center,y_center)
# embedded in a background flow in the positive x-direction.

x_center = 0
y_center = 0
radius = lambda x, y: np.sqrt((x - x_center)**2 + (y - y_center)**2)

height = 20
decay = 8e-11
h[:,:] = 10e3 - height * np.exp(-decay * radius(mass_xs, mass_ys)**2)

f = config['physics_parameters']['f']
g = config['physics_parameters']['g']

r = radius(mass_xs, face_ys)
u[:,:] = -(face_ys - y_center) * np.exp(-decay * r**2)
u[:,:] /= np.sqrt(u**2 + ((mass_xs - x_center) * np.exp(-decay * r**2))**2)
u[:,:] *= -f*r/2 + np.sqrt((f*r/2)**2 + (g * 2*height*decay * r**2 * np.exp(-decay*r**2.)))

u[:,:] += 10

r = radius(face_xs, mass_ys)
v[:,:] = +(face_xs - x_center) * np.exp(-decay * r**2)
v[:,:] /= np.sqrt(v**2 + ((mass_ys - y_center) * np.exp(-decay * r**2))**2)
v[:,:] *= -f*r/2 + np.sqrt((f*r/2)**2 + (g * 2*height*decay * r**2 * np.exp(-decay*r**2.)))

#=======================================================
# Create output file

xr.Dataset(
    data_vars=dict(
        h=(['xc', 'yc', 't'], h[..., np.newaxis]),
        u=(['xc', 'yf', 't'], u[..., np.newaxis]),
        v=(['xf', 'yc', 't'], v[..., np.newaxis]),
    ),
    coords=dict(
        t=('t', [0]),
        xc=('xc', mass_xs[0,:]),
        yc=('yc', mass_ys[:,0]),
        xf=('xf', face_xs[0,:]),
        yf=('yf', face_ys[:,0]),
    ),
).to_netcdf('initial.nc', unlimited_dims=['t'])