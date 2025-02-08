import warnings
import numpy as np
import xarray as xr
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.axes_grid1 import Divider, Size

data = xr.open_dataset('output/output.nc')

#=======================================================
# Grid definitions

# make sure these values match the config
Lx = (data['x'].max() - data['x'].min() + data['x'].diff('x')[0]).values
Ly = (data['y'].max() - data['y'].min() + data['y'].diff('y')[0]).values
gravity = 9.807
coriolis = 1e-5

xs = data.x.values / (Lx / 2)
ys = data.y.values / (Ly / 2)
xgrid, ygrid = np.meshgrid(xs, ys)

#=======================================================
# Definitions of quantities to plot

mass = data.h.sum(['x', 'y'])
x_momentum = (data.h * data.u).sum(['x', 'y'])
y_momentum = (data.h * data.v).sum(['x', 'y'])
potential_energy = .5 * gravity * (data.h**2).sum(['x', 'y'])
kinetic_energy = .5 * (data.h * (data.u**2 + data.v**2)).sum(['x', 'y'])
total_energy = potential_energy + kinetic_energy
absolute_vorticity = (
    coriolis - data.u.differentiate('y') + data.v.differentiate('x')).sum(['x','y'])

#=======================================================
# Plotting setup -- axes, normalisation, colorbars, etc.

fig_height = 6 # figure height, inches
nrows = 5 # 1 map + `nrows-1` time series
map_ratio = 10 # the map is `map_ratio` times taller than each of the time series
x_span = (.12, .88) # start and end x coordinates in figure fraction
y_span = (.08, .95) # start and end y coordinates in figure fraction
balance = .15

y_spacing = .05
cbar_fraction = .025

h = (y_span[1] - y_span[0] - (nrows-1)*y_spacing) / (map_ratio + nrows - 1)
fig_width = fig_height * (map_ratio * h)/(x_span[1] - x_span[0])
fig = plt.figure(figsize=(fig_width, fig_height))

rows = [Size.Scaled(y_span[0]), Size.Scaled(h),
        *[Size.Scaled(y_spacing*(1-balance)), Size.Scaled(h)] * (nrows-2),
        Size.Scaled(y_spacing*(1+balance*(nrows-2))), Size.Scaled(map_ratio*h),
        Size.Scaled(1-y_span[1])]
cols = [Size.Scaled(x_span[0]), Size.Scaled(x_span[1]-x_span[0]-1.5*cbar_fraction),
        Size.Scaled(cbar_fraction/2),
        Size.Scaled(cbar_fraction), Size.Scaled(1-x_span[1])]

divider = Divider(fig, (0, 0, 1, 1), cols, rows, aspect=False)

map_ax = fig.add_axes(
    divider.get_position(),
    axes_locator=divider.new_locator(nx=1, ny=2*nrows-1))

map_ax.tick_params(labelsize='small')
map_ax.set_xticks([])
map_ax.set_yticks([])
map_ax.set_xlabel(r'$x\rightarrow$', fontsize='small', loc='left')
map_ax.set_ylabel(r'$y\rightarrow$', fontsize='small', loc='bottom')
map_ax.set_xlim(-1, +1)
map_ax.set_ylim(-1, +1)
map_ax.set_title(f'time$=0.00$ hours', fontsize='small', loc='left', pad=0)

cbar_ax = fig.add_axes(
    divider.get_position(),
    axes_locator=divider.new_locator(nx=3, ny=2*nrows-1))

cbar_ax.tick_params(labelsize='small')

line_axs = []
for i in range(nrows-1):
    ax = fig.add_axes(
        divider.get_position(),
        axes_locator=divider.new_locator(nx=1, nx1=4, ny=2*i+1))
    ax.tick_params(labelsize='small')
    ax.ticklabel_format(axis='y', scilimits=(0,0))
    ax.yaxis.get_offset_text().set_fontsize('small')
    if i != 0:
        ax.set_xticklabels([])
    line_axs.append(ax)

line_axs[0].set_xlabel('$t$ (hours)', fontsize='small')

hrange = ((data.h - data.h.mean())**2).mean() ** .5
norm = mpl.colors.CenteredNorm(vcenter=0, halfrange=hrange)
cmap = 'bwr'
plt.colorbar(
    mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax, fraction=.05, pad=.02,
    format=lambda x, _: rf'${x:.0f}$' if x==0 else rf'${x:+.1f}$' if int(x)!=x else rf'${x:+.0f}$')

quiverskip = 5
windspeed = np.sqrt(data.u**2 + data.v**2)

mesh = map_ax.pcolormesh(xs, ys, np.zeros_like(data.h[0,:,:]), norm=norm, cmap=cmap)
quiver = map_ax.quiver(
    xgrid[::quiverskip,::quiverskip], ygrid[::quiverskip,::quiverskip],
    data.u.values[0,::quiverskip,::quiverskip], data.v.values[0,::quiverskip,::quiverskip], 
    scale=10*windspeed[::quiverskip,::quiverskip].max().values)
map_ax.quiverkey(quiver, X=.8, Y=1.02, U=50, label='$50$ m s$^{-1}$',
                 labelpos='E', fontproperties=dict(size='small'))

momentum_magnitude = np.sqrt(x_momentum**2 + y_momentum**2)
timeseries = [
    (3, 'Fractional mass anomaly', mass/mass[0]-1, dict(c='k', ls='-')),
    (2, None, (x_momentum - x_momentum[0]) / momentum_magnitude[0], dict(c='r', ls='-')),
    (2, 'Fractional momentum anomaly', (y_momentum - y_momentum[0]) / momentum_magnitude[0], dict(c='g', ls='-')),
    # (1, None, potential_energy/total_energy[0], dict(c='b', ls='-')),
    # (1, None, kinetic_energy/total_energy[0], dict(c='r', ls='-')),
    (1, 'Fractional energy anomaly', total_energy/total_energy[0] - 1, dict(c='k', ls='-', lw=2)),
    (0, 'Fractional absolute vorticity anomaly', absolute_vorticity/absolute_vorticity[0] - 1, dict(c='k', ls='-')),
]

timeseries_lines = []

for axn, name, d, args in timeseries:
    if name:
        line_axs[axn].set_title(name, fontsize='x-small', loc='right', pad=0)
    line, = line_axs[axn].plot(d.t[0], d[0], **args)
    timeseries_lines.append(line)

#=======================================================
# Create the animation

window = np.inf # hours

def update(i):
    ti = data.t.values[i]
    map_ax.set_title(f'time$={ti/3600:.2f}$ hours', fontsize='small', loc='left', pad=0)
    mesh.set_array((data.h.sel(t=ti) - data.h[0,:,:]).values.ravel())
    quiver.set_UVC(data.u.values[i,::quiverskip,::quiverskip],
                   data.v.values[i,::quiverskip,::quiverskip])
    time = data.t.sel(t=slice(max(data.t.min(), ti-3600*window), ti))
    if i > 0:
        for (axn, _, d, _), line in zip(timeseries, timeseries_lines):
            ax = line_axs[axn]
            line.set_data(time / 3600, d.sel(t=time))
            low = d.sel(t=time).min()
            high = d.sel(t=time).max()
            widen = .1 * abs(low-high) if low != high else .01 * low if low != 0 else 1e-9
            if i > 1:
                current_low, current_high = ax.get_ylim()
                ax.set_ylim(bottom=min(current_low, low - widen), top=max(current_high, high + widen))
            else:
                ax.set_ylim(bottom=low - widen, top=high + widen)
            if time.max() > time.min():
                ax.set_xlim(time.min() / 3600, time.max() / 3600)

ani = FuncAnimation(fig, update, frames=range(len(data.t)), interval=50)
# ani.save('animation.mp4', dpi=150)
plt.show()