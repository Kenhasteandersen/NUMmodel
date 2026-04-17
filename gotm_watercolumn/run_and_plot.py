#!/usr/bin/env python3
"""
Run the NUMmodel/GOTM water column example and plot all output variables.

Run from the gotm_watercolumn/ directory with the GOTM executable already built:
    python3 run_and_plot.py [--gotm PATH]

The script:
  1. Generates synthetic meteorological forcing (meteo.dat)
  2. Runs the GOTM simulation
  3. Plots all output state variables as depth-time panels
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--gotm', default='build/gotm',
                    help='path to GOTM executable (default: build/gotm)')
args = parser.parse_args()

here = Path(__file__).parent.resolve()
os.chdir(here)

gotm_exe = Path(args.gotm)
if not gotm_exe.exists():
    sys.exit(f'ERROR: GOTM executable not found at {gotm_exe}\n'
             f'Build it first — see README.md.')

# ---------------------------------------------------------------------------
# 1. Generate forcing
# ---------------------------------------------------------------------------
print('==> Generating meteorological forcing ...')
subprocess.run([sys.executable, 'generate_forcing.py'], check=True)

# ---------------------------------------------------------------------------
# 2. Run GOTM
# ---------------------------------------------------------------------------
print('==> Running GOTM ...')
subprocess.run([str(gotm_exe)], check=True)

# ---------------------------------------------------------------------------
# 3. Plot
# ---------------------------------------------------------------------------
print('==> Plotting output ...')

ds = nc.Dataset('NUMmodel_watercolumn.nc')

t_dates = [datetime(2000, 1, 1) + timedelta(seconds=float(s))
           for s in ds.variables['time'][:]]
t_num = np.array(mdates.date2num(t_dates))
z     = np.array(ds.variables['z'][0, :, 0, 0])

variables = [
    ('temp',               'Temperature',              '°C',           'RdYlBu_r'),
    ('salt',               'Salinity',                 'PSU',          'viridis'),
    ('num_model_N',        'DIN (N)',                  'µg N L⁻¹',    'viridis'),
    ('num_model_DOC',      'DOC',                      'µg C L⁻¹',    'YlOrBr'),
    ('num_model_Si',       'Silicate (Si)',             'µg Si L⁻¹',   'PuBuGn'),
    ('num_model_ProdGross', 'Gross production',         'mg C m⁻³ d⁻¹','YlGn'),
    ('num_model_ProdNet',  'Net production',            'mg C m⁻³ d⁻¹','YlGn'),
    ('num_model_ProdHTL',  'HTL production',            'mg C m⁻³ d⁻¹','Oranges'),
    ('num_model_Bpico',    'Pico biomass',              'µg C L⁻¹',    'Blues'),
    ('num_model_Bnano',    'Nano biomass',              'µg C L⁻¹',    'Greens'),
    ('num_model_Bmicro',   'Micro biomass',             'µg C L⁻¹',    'RdPu'),
    ('num_model_BGen1',    'Generalists (total)',       'µg C L⁻¹',    'YlGn'),
    ('num_model_BDia1',    'Diatoms (total)',           'µg C L⁻¹',    'PuBuGn'),
    ('num_model_BPCop1',   'Passive copepods grp 1',   'µg C L⁻¹',    'BuPu'),
    ('num_model_BPCop2',   'Passive copepods grp 2',   'µg C L⁻¹',    'BuPu'),
    ('num_model_BACop1',   'Active copepods grp 1',    'µg C L⁻¹',    'OrRd'),
    ('num_model_BACop2',   'Active copepods grp 2',    'µg C L⁻¹',    'OrRd'),
    ('num_model_BACop3',   'Active copepods grp 3',    'µg C L⁻¹',    'OrRd'),
    ('num_model_BPOM1',    'POM',                      'µg C L⁻¹',    'YlOrBr'),
]

data = {v[0]: np.array(ds.variables[v[0]][:, :, 0, 0]) for v in variables}
ds.close()

ncols = 3
nrows = -(-len(variables) // ncols)

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

for idx in range(len(variables), nrows * ncols):
    axes.flat[idx].set_visible(False)

fig.suptitle('NUMmodel / GOTM — 59°N, 2000–2001', fontsize=12, y=1.01)
plt.tight_layout()

outfile = 'all_variables.png'
plt.savefig(outfile, dpi=150, bbox_inches='tight')
print(f'==> Saved {outfile}')

# ---------------------------------------------------------------------------
# Plot 2: depth-integrated time series — all state variables in one panel
# ---------------------------------------------------------------------------
# z is negative (surface ~ 0, bottom ~ -100 m); trapz with negative z gives
# a negative integral, so negate to get positive column totals.
_skip_integrated = {
    'num_model_ProdGross', 'num_model_ProdNet', 'num_model_ProdHTL',
    'num_model_Bpico', 'num_model_Bnano', 'num_model_Bmicro',
}
integrated_vars = [(v, t, u, c) for v, t, u, c in variables if v not in _skip_integrated]

integrated = {}
for varname, *_ in integrated_vars:
    arr = data[varname]                          # (time, z)
    integrated[varname] = -np.trapezoid(arr, z, axis=1)

colors = plt.get_cmap('tab20').colors

fig2, ax2 = plt.subplots(figsize=(14, 5))

for idx, (varname, title, units, cmap) in enumerate(integrated_vars):
    ax2.plot(t_num, integrated[varname], label=title,
             color=colors[idx % len(colors)], linewidth=1.4)

ax2.set_yscale('symlog', linthresh=1.0)
ax2.xaxis_date()
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b %Y'))
ax2.xaxis.set_major_locator(mdates.MonthLocator(bymonth=[1, 4, 7, 10]))
ax2.tick_params(labelsize=8)
ax2.set_ylabel('Depth-integrated value (unit × m)', fontsize=9)
ax2.set_title('NUMmodel / GOTM — depth-integrated state variables', fontsize=10)
ax2.legend(fontsize=7, ncol=2, loc='upper right', framealpha=0.7)
plt.tight_layout()

outfile2 = 'integrated_variables.png'
fig2.savefig(outfile2, dpi=150, bbox_inches='tight')
print(f'==> Saved {outfile2}')

plt.show()
