import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = xr.open_dataset('output/output.nc')

gravity = 9.807 # make sure this matches the config

hmean = data.h.mean().values
hrange = (data.h.max() - data.h.min()).values / 2

norm = mpl.colors.CenteredNorm(vcenter=hmean, halfrange=hrange)
cmap = 'bwr'

fig, ax = plt.subplots()
ax.set_aspect('equal')
mesh = ax.pcolormesh(data.x, data.y, data.h[0,:,:], norm=norm, cmap=cmap)
plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)

def init():
    ax.set_title(f'time$=0$')
    ax.set_xlim(data.x.min(), data.x.max())
    ax.set_ylim(data.y.min(), data.y.max())
    return mesh,

def update(t):
    ax.set_title(f'time$={t/3600:.2f}$ h')
    mesh.set_array(data.h.sel(t=t).values.ravel())
    return mesh,

ani = FuncAnimation(fig, update, frames=data.t.values, init_func=init, interval=20)
plt.show()

mass = data.h.sum(['x', 'y'])
x_momentum = (data.h * data.u).sum(['x', 'y'])
y_momentum = (data.h * data.v).sum(['x', 'y'])
energy = .5 * (gravity * data.h**2 + data.h * (data.u**2 + data.v**2)).sum(['x', 'y'])

time = data.t / 3600
fig, axs = plt.subplots(3, 1, sharex=True)
axs[0].plot(time, mass - mass.isel(t=0))
axs[1].plot(time, x_momentum - x_momentum.isel(t=0))
axs[1].plot(time, y_momentum - y_momentum.isel(t=0))
axs[2].plot(time, energy - energy.isel(t=0))
axs[2].set_xlim(time.min(), time.max())
axs[2].set_xlabel(r'$t$ (h)')
plt.show()