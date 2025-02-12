import f90nml
import numpy as np
import xarray as xr

config = f90nml.read('config.nml')

#=======================================================
# Domain information

nx = 200
ny = 200
Lx = 15000e3
Ly = 15000e3

nlev = 2

dx = Lx / nx
dy = Ly / ny

mass_xs, mass_ys = np.meshgrid(
    np.linspace((dx-Lx)/2, (Lx-dx)/2, nx),
    np.linspace((dy-Ly)/2, (Ly-dy)/2, ny))
face_xs = mass_xs - dx/2
face_ys = mass_ys - dy/2

#=======================================================
# Initialise data arrays

dp = np.zeros((ny,nx,nlev)) # pressure thickness
pt = np.zeros((ny,nx,nlev)) # potential temperature
u = np.zeros((ny,nx,nlev)) # eastward wind
v = np.zeros((ny,nx,nlev)) # northward wind
gzs = np.zeros((ny,nx)) # surface geopotential height

#=======================================================
# Set data arrays

# We initialise the atmosphere to be in hydrostatic balance over some
# specified terrain with a constant surface temperature and lapse rate

f = config['physics_parameters']['coriolis_parameter']
g = config['physics_parameters']['gravity']
top_pressure = config['physics_parameters']['top_pressure']
reference_pressure = config['physics_parameters']['reference_pressure']
dry_gas_constant = config['physics_parameters']['dry_gas_constant']
dry_heat_capacity = config['physics_parameters']['dry_heat_capacity']
kappa = dry_gas_constant / dry_heat_capacity
hydrostatic_constant = dry_heat_capacity / reference_pressure**kappa

surface_temperature = 288 # K = 15 C
lapse_rate = 20/3e3 # K / m
top_height = 10e3

repeat_along_z = lambda a, n: np.repeat(a[:,:,np.newaxis], n, axis=2)

surface_height = (
    1e3 * np.exp(-(mass_xs**2 + mass_ys**2) / (Lx/8)**2)
    * (1 + np.cos(2*np.pi * np.sqrt(mass_xs**2 + mass_ys**2) / (Lx/10))))

pressure = np.zeros((ny,nx,nlev+1))
pressure_kappa = np.zeros((ny,nx,nlev+1))

pressure[:,:,0] = top_pressure * (
    (surface_temperature - lapse_rate * surface_height)
    / (surface_temperature - lapse_rate * top_height)
)**(g/lapse_rate/dry_gas_constant)

# we split the pressure difference between top and bottom into `nlev` equal layers
dp[:,:,:] = repeat_along_z((pressure[:,:,0] - top_pressure) / nlev, nlev)

# calculate the pressure on the layer interfaces
pressure[:,:,-1] = top_pressure
pressure[:,:,:-1] = top_pressure + np.cumsum(dp[:,:,::-1], axis=2)[:,:,::-1]

# calculate the geopotential on interfaces for our constant lapse rate atmosphere
geopotential = np.zeros((ny,nx,nlev+1))
geopotential[:,:,0] = g * surface_height[:,:]
geopotential[:,:,1:] = (g/lapse_rate) * (
    surface_temperature
    - (surface_temperature - lapse_rate * repeat_along_z(surface_height, nlev))
    * (pressure[:,:,1:] / repeat_along_z(pressure[:,:,0], nlev))
    ** (lapse_rate*dry_gas_constant/g)
)

# calculate the average potential temperature for each layer given the calculated
# geopotential and pressure thicknesses
pt[:,:,:] = abs(
    np.diff(geopotential)
    / np.diff(pressure**kappa, axis=2)
    / hydrostatic_constant)

u[:,:,:] = 0
v[:,:,:] = 0

#=======================================================
# Create output file

xr.Dataset(
    data_vars=dict(
        dp=(['yc', 'xc', 'lev', 't'], dp[..., np.newaxis]),
        pt=(['yc', 'xc', 'lev', 't'], pt[..., np.newaxis]),
        u=(['yc', 'xf', 'lev', 't'], u[..., np.newaxis]),
        v=(['yf', 'xc', 'lev', 't'], v[..., np.newaxis]),
        gzs=(['yc', 'xc', 't'], geopotential[:,:,0, np.newaxis]),
    ),
    coords=dict(
        t=('t', [0]),
        xc=('xc', mass_xs[0,:]),
        yc=('yc', mass_ys[:,0]),
        xf=('xf', face_xs[0,:]),
        yf=('yf', face_ys[:,0]),
        lev=('lev', list(range(1,nlev+1)))
    ),
).to_netcdf('initial.nc', unlimited_dims=['t'])