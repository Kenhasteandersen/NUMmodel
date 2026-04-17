#!/usr/bin/env python3
"""One panel per output state variable (depth × time)."""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta

ds = nc.Dataset('NUMmodel_watercolumn.nc')

t_dates = [datetime(2000, 1, 1) + timedelta(seconds=float(s))
           for s in ds.variables['time'][:]]
t_num   = np.array(mdates.date2num(t_dates))
z       = np.array(ds.variables['z'][0, :, 0, 0])

variables = [
    ('temp',              'Temperature',       '°C',           'RdYlBu_r'),
    ('salt',              'Salinity',           'PSU',          'viridis'),
    ('num_model_N',       'DIN (N)',            'µg N L⁻¹',    'viridis'),
    ('num_model_DOC',     'DOC',                'µg C L⁻¹',    'YlOrBr'),
    ('num_model_Si',      'Silicate (Si)',       'µg Si L⁻¹',   'PuBuGn'),
    ('num_model_ProdGross','Gross production',   'mg C m⁻³ d⁻¹','YlGn'),
    ('num_model_ProdNet', 'Net production',      'mg C m⁻³ d⁻¹','YlGn'),
    ('num_model_ProdHTL', 'HTL production',      'mg C m⁻³ d⁻¹','Oranges'),
    ('num_model_Bpico',   'Pico biomass',        'mg C m⁻³',    'Blues'),
    ('num_model_Bnano',   'Nano biomass',        'mg C m⁻³',    'Greens'),
    ('num_model_Bmicro',  'Micro biomass',       'mg C m⁻³',    'RdPu'),
]

data = {v[0]: np.array(ds.variables[v[0]][:, :, 0, 0]) for v in variables}
ds.close()

ncols = 3
nrows = -(-len(variables) // ncols)   # ceiling division

fig, axes = plt.subplots(nrows, ncols, figsize=(14, 3.2 * nrows),
                         sharex=True, sharey=True)

for idx, (varname, title, units, cmap) in enumerate(variables):
    ax  = axes.flat[idx]
    arr = data[varname]

    vmin = np.nanpercentile(arr, 2)
    vmax = np.nanpercentile(arr, 98)
    if vmax <= vmin:
        vmax = vmin + 1e-9

    pcm = ax.pcolormesh(t_num, z, arr.T, shading='nearest',
                        cmap=cmap, vmin=vmin, vmax=vmax)
    cbar = fig.colorbar(pcm, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label(units, fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    ax.set_title(title, fontsize=9, fontweight='bold')
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.tick_params(labelsize=7)
    if idx % ncols == 0:
        ax.set_ylabel('Depth (m)', fontsize=8)

# hide unused panels
for idx in range(len(variables), nrows * ncols):
    axes.flat[idx].set_visible(False)

fig.suptitle('NUMmodel / GOTM — 59°N, year 2000', fontsize=12, y=1.01)
plt.tight_layout()
plt.savefig('all_variables.png', dpi=150, bbox_inches='tight')
print('Saved all_variables.png')
plt.show()
