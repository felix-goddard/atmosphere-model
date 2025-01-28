import warnings
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

xs = data.x.values / (Lx / 2)
ys = data.y.values / (Ly / 2)
xgrid, ygrid = np.meshgrid(xs, ys)

mass = data.h.sum(['x', 'y'])
x_momentum = (data.h * data.u).sum(['x', 'y'])
y_momentum = (data.h * data.v).sum(['x', 'y'])
potential_energy = .5 * gravity * (data.h**2).sum(['x', 'y'])
kinetic_energy = .5 * (data.h * (data.u**2 + data.v**2)).sum(['x', 'y'])
total_energy = potential_energy + kinetic_energy
potential_vorticity = (
    (coriolis - data.u.differentiate('y') + data.v.differentiate('x')) / data.h).sum(['x','y'])

hmean = data.h.mean().values
hrange = (data.h.max() - data.h.min()).values / 2

fig, axs = plt.subplots(5,1, figsize=(5,7), layout='constrained', height_ratios=[10,1,1,1,1])

axs[0].set_aspect('equal')
norm = mpl.colors.CenteredNorm(vcenter=hmean, halfrange=hrange)
cmap = 'bwr'
mesh = axs[0].pcolormesh(xs, ys, data.h[0,:,:], norm=norm, cmap=cmap)
plt.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axs[0], fraction=.05, pad=.02)

quiverskip = 5
quiver = axs[0].quiver(xgrid[::quiverskip,::quiverskip],
                       ygrid[::quiverskip,::quiverskip],
                       data.u.values[0,::quiverskip,::quiverskip],
                       data.v.values[0,::quiverskip,::quiverskip], scale=400)

mass_line, = axs[1].plot(data.t[0], mass[0], 'k-')
momx_line, = axs[2].plot(data.t[0], x_momentum[0], 'r-')
momy_line, = axs[2].plot(data.t[0], y_momentum[0], 'g-')
pe_line, = axs[3].plot(data.t[0], potential_energy[0], 'b-')
ke_line, = axs[3].plot(data.t[0], kinetic_energy[0], 'r-')
te_line, = axs[3].plot(data.t[0], total_energy[0], 'k-', lw=2)
pv_line, = axs[4].plot(data.t[0], potential_vorticity[0], 'k-')

def init():
    axs[0].set_title(f'time$=0$')
    axs[0].set_xlim(xs.min(), xs.max())
    axs[0].set_ylim(ys.min(), ys.max())
    axs[1].set_title('Mass anomaly', fontsize='x-small')
    axs[2].set_title('Momentum anomaly', fontsize='x-small')
    axs[3].set_title('Energy anomaly', fontsize='x-small')
    axs[4].set_title('Potential vorticity anomaly', fontsize='x-small')
    axs[4].set_xlabel(r'$t$ (h)')
    return mesh,

def update(i):
    t = data.t.values[i]
    axs[0].set_title(f'time$={t/3600:.2f}$ h')
    mesh.set_array(data.h.sel(t=t).values.ravel())
    quiver.set_UVC(data.u.values[i,::quiverskip,::quiverskip],
                   data.v.values[i,::quiverskip,::quiverskip])
    time = data.t[:i] / 3600
    if i > 0:
        for line, d, ax in [
            (mass_line, mass, axs[1]),
            (momx_line, x_momentum, axs[2]),
            (momy_line, y_momentum, axs[2]),
            (pe_line, potential_energy, axs[3]),
            (ke_line, kinetic_energy, axs[3]),
            (te_line, total_energy, axs[3]),
            (pv_line, potential_vorticity, axs[4]),
        ]:
            line.set_data(time, d[:i] - d[0])
            low = (d[:i] - d[0]).min()
            high = (d[:i] - d[0]).max()
            widen = .1 * abs(low-high) if low != high else 1
            if i > 1:
                current_low, current_high = ax.get_ylim()
                ax.set_ylim(bottom=min(current_low, low - widen), top=max(current_high, high + widen))
            else:
                ax.set_ylim(bottom=low - widen, top=high + widen)
        if time.max() > time.min():
            for i in range(1,5):
                axs[i].set_xlim(time.min(), time.max())
                if i != 4: axs[i].set_xticklabels([])

ani = FuncAnimation(fig, update, frames=range(len(data.t)), init_func=init, interval=50)
# ani.save('animation.mp4', dpi=150)
plt.show()