import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

data = xr.open_dataset('output/output.nc')

# make sure these values match the config
Lx = Ly = 1500e3
gravity = 9.807
coriolis = 1e-5

hmean = data.h.mean().values
hrange = (data.h.max() - data.h.min()).values / 2

xs = data.x.values / (Lx / 2)
ys = data.y.values / (Ly / 2)

xgrid, ygrid = np.meshgrid(xs, ys)

norm = mpl.colors.CenteredNorm(vcenter=hmean, halfrange=hrange)
cmap = 'bwr'

fig, ax = plt.subplots(layout='constrained')
ax.set_aspect('equal')
mesh = ax.pcolormesh(xs, ys, data.h[0,:,:], norm=norm, cmap=cmap)
# plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=ax)

quiverskip = 5
quiver = ax.quiver(xgrid[::quiverskip,::quiverskip],
                   ygrid[::quiverskip,::quiverskip],
                   data.u.values[0,::quiverskip,::quiverskip],
                   data.v.values[0,::quiverskip,::quiverskip], scale=1)

def init():
    ax.set_title(f'time$=0$')
    ax.set_xlim(xs.min(), xs.max())
    ax.set_ylim(ys.min(), ys.max())
    return mesh,

def update(i):
    t = data.t.values[i]
    ax.set_title(f'time$={t/3600:.2f}$ h')
    mesh.set_array(data.h.sel(t=t).values.ravel())
    quiver.set_UVC(data.u.values[i,::quiverskip,::quiverskip],
                   data.v.values[i,::quiverskip,::quiverskip])
    return mesh, quiver,

ani = FuncAnimation(fig, update, frames=range(len(data.t)), init_func=init, interval=1)
# ani.save('animation.mp4')
plt.show()

mass = data.h.sum(['x', 'y'])
x_momentum = (data.h * data.u).sum(['x', 'y'])
y_momentum = (data.h * data.v).sum(['x', 'y'])
energy = .5 * (gravity * data.h**2 + data.h * (data.u**2 + data.v**2)).sum(['x', 'y'])
potential_vorticity = ((coriolis - data.u.differentiate('y') + data.v.differentiate('x')) / data.h).sum(['x','y'])

roll = lambda x: x.rolling(t=180, center=True).mean()
time = roll(data.t) / 3600
fig, axs = plt.subplots(4, 1, layout='constrained', sharex=True)
axs[0].plot(time, roll(mass) / mass.isel(t=0) - 1)
axs[1].plot(time, roll(x_momentum) - x_momentum.isel(t=0))
axs[1].plot(time, roll(y_momentum) - y_momentum.isel(t=0))
axs[2].plot(time, roll(energy) / energy.isel(t=0) - 1)
axs[3].plot(time, roll(potential_vorticity) / potential_vorticity.isel(t=0) - 1)
axs[3].set_xlim(time.min(), time.max())
axs[3].set_xlabel(r'$t$ (h)')
plt.show()