import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = xr.open_dataset('output/output.nc')

hmean = 10e3 # data.h.mean().values
hrange = 5 # (data.h.max() - data.h.min()).values / 2

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
