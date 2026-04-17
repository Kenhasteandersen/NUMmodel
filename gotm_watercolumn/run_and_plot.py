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

fig.suptitle('NUMmodel / GOTM — 59°N, year 2000', fontsize=12, y=1.01)
plt.tight_layout()

outfile = 'all_variables.png'
plt.savefig(outfile, dpi=150, bbox_inches='tight')
print(f'==> Saved {outfile}')
plt.show()
