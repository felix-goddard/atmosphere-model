import f90nml
import numpy as np
import xarray as xr

config = f90nml.read('config.nml')

# ==============================================================================
# Domain and setup information

# Horizontal grid --------------------------------------------------------------

nx = 200
ny = 200
Lx = 15000e3
Ly = 15000e3

dx = Lx / nx
dy = Ly / ny

mass_xs, mass_ys = np.meshgrid(
    np.linspace((dx-Lx)/2, (Lx-dx)/2, nx),
    np.linspace((dy-Ly)/2, (Ly-dy)/2, ny))
face_xs = mass_xs - dx/2
face_ys = mass_ys - dy/2

# Vertical grid ----------------------------------------------------------------

nlay = 20

top_pressure = 26250  # Pa
top_height = 10e3  # m

# the vertical coordinate we remap to is a hybrid σ-p coordinate, defined as
#   p_k = A_k + B_k * (p_surface - p_top)
# we compute the coefficients A_k and B_k according to some target set of
# pressure thicknesses and the algorithm of Eckermann (2009)

kp = 5  # number of isobaric layers at model top
ks = 0  # number of sigma layers at model bottom

rp = 2.2  # coefficients of the hybrid σ-p coordinate, see Eckermann (2009)
rs = 1.2
S = 5

ref_ps = 100_000  # this should be close to surface pressures encountered
                  # in the model run, Pa

splay = .5  # the thickess of pressure layers at the top and bottom relative to
            # the middle of the model domain
dp_func = lambda p: 1 - 4 * (1 - splay) * (p - .5)**2

dp = dp_func(np.linspace(0, 1, nlay))
dp = dp * (ref_ps - top_pressure) / dp.sum()

plev = np.ones(nlay+1) * top_pressure
plev[:-1] += np.cumsum(dp[::-1])[::-1]

eta = (plev - top_pressure) / (ref_ps - top_pressure)

Ak = np.zeros(nlay+1)
Bk = np.zeros(nlay+1)

bk = (eta - eta[nlay-kp]) / (1 - eta[nlay-kp])
rk = rp + (rs - rp) * np.arctan(S * bk) / np.arctan(S)

if kp > 0:
    for k in range(nlay, nlay-kp, -1):
        Ak[k] = (top_pressure + (eta[k]/eta[nlay-kp])
                 * (plev[nlay-kp] - top_pressure))
        Bk[k] = 0

for k in range(nlay-kp, ks-1, -1):
    Bk[k] = bk[k]**rk[k]
    Ak[k] = (top_pressure + (eta[k] - Bk[k])
             * (ref_ps - top_pressure))

if ks > 0:
    for k in range(ks-1, -1, -1):
        Bk[k] = bk[k]
        Ak[k] = (top_pressure + (eta[k] - Bk[k])
                 * (ref_ps - top_pressure))

# Physics constants ------------------------------------------------------------

R_dry = 287.052874  # specific gas constant of dry air, J/kg/K
cp_dry = 1005  # constant pressure heat capacity of dry air, J/kg/K

f = config['physics_parameters']['coriolis_parameter']
g = config['physics_parameters']['gravity']
kappa = R_dry / cp_dry

lapse_rate = 20/3e3  # initial environmental lapse rate, K/m

# ==============================================================================
# Initialise data arrays

dp = np.zeros((ny, nx, nlay))  # pressure thickness
pt = np.zeros((ny, nx, nlay))  # potential temperature
u = np.zeros((ny, nx, nlay))  # eastward wind
v = np.zeros((ny, nx, nlay))  # northward wind
gzs = np.zeros((ny, nx))  # surface geopotential height
ts = np.zeros((ny, nx))  # surface temperature

# ==============================================================================
# Set data arrays

def repeat_along_z(a, n): return np.repeat(a[:, :, np.newaxis], n, axis=2)

# We initialise the atmosphere to be in hydrostatic balance over some
# specified terrain with a constant surface temperature and lapse rate

ts[:, :] = 288  # K = 15 C

surface_height = np.zeros((ny, nx))

pressure = np.zeros((ny, nx, nlay+1))
pressure_kappa = np.zeros((ny, nx, nlay+1))

pressure[:, :, 0] = top_pressure * (
    (ts - lapse_rate * surface_height) / (ts - lapse_rate * top_height)
)**(g/lapse_rate/R_dry)

pressure[:, :, :] = (
    Ak + np.einsum('ij,k->ijk', pressure[:, :, 0] - top_pressure, Bk))

dp = -np.diff(pressure, axis=2)

# calculate the geopotential on interfaces for our constant lapse rate atmosphere
geopotential = np.zeros((ny, nx, nlay+1))
geopotential[:, :, 0] = g * surface_height[:, :]
geopotential[:, :, 1:] = (g/lapse_rate) * (
    repeat_along_z(ts, nlay)
    - repeat_along_z(ts - lapse_rate * surface_height, nlay)
    * (pressure[:, :, 1:] / repeat_along_z(pressure[:, :, 0], nlay))
    ** (lapse_rate*R_dry/g)
)

# calculate the average potential temperature for each layer given the calculated
# geopotential and pressure thicknesses
pt[:, :, :] = abs(
    np.diff(geopotential)
    / np.diff(pressure**kappa, axis=2)
    / cp_dry)

u[:, :, :] = 0
v[:, :, :] = 0

# ==============================================================================
# Create file

xr.Dataset(
    data_vars=dict(
        dp=(['yc', 'xc', 'lay', 't'], dp[..., np.newaxis]),
        pt=(['yc', 'xc', 'lay', 't'], pt[..., np.newaxis]),
        u=(['yc', 'xf', 'lay', 't'], u[..., np.newaxis]),
        v=(['yf', 'xc', 'lay', 't'], v[..., np.newaxis]),
        gzs=(['yc', 'xc', 't'], geopotential[:, :, 0, np.newaxis]),
        ts=(['yc', 'xc', 't'], (ts[:, :])[..., np.newaxis]),
        ak=(['lev'], Ak),
        bk=(['lev'], Bk),
    ),
    coords=dict(
        t=('t', [0]),
        xc=('xc', mass_xs[0, :]),
        yc=('yc', mass_ys[:, 0]),
        xf=('xf', face_xs[0, :]),
        yf=('yf', face_ys[:, 0]),
        lay=('lay', list(range(1, nlay+1))),
        lev=('lev', list(range(1, nlay+2))),
    ),
).to_netcdf('initial.nc', unlimited_dims=['t'])
