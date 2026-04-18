#!/usr/bin/env python3
"""
Run the NUMmodel/GOTM OWS Papa example and plot output.

Run from the gotm_ows_papa/ directory with the GOTM executable already built:
    python3 run_and_plot.py [--gotm PATH]
"""

import argparse
import subprocess
import sys
import os
from pathlib import Path

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.dates as mdates
from datetime import datetime, timedelta

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
parser = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('--gotm', default='../gotm_build/gotm',
                    help='path to GOTM executable (default: ../gotm_build/gotm)')
args = parser.parse_args()

here = Path(__file__).parent.resolve()
os.chdir(here)

gotm_exe = Path(args.gotm)
if not gotm_exe.exists():
    sys.exit(f'ERROR: GOTM executable not found at {gotm_exe}\n'
             f'Build it first — see README.md.')

# ---------------------------------------------------------------------------
# Run GOTM
# ---------------------------------------------------------------------------
print('==> Running GOTM ...')
subprocess.run([str(gotm_exe)], check=True)

# ---------------------------------------------------------------------------
# Load output
# ---------------------------------------------------------------------------
print('==> Plotting output ...')

ds = nc.Dataset('ows_papa_NUMmodel.nc')

t_dates = [datetime(2011, 3, 21) + timedelta(seconds=float(s))
           for s in ds.variables['time'][:]]
t_num = np.array(mdates.date2num(t_dates))
z     = np.array(ds.variables['z'][0, :, 0, 0])          # layer centres (150,)
nuh   = np.array(ds.variables['nuh'][:, :, 0, 0])        # interfaces  (time, 151)

# Interface depths for nuh
z_i = np.concatenate([[0.0], 0.5*(z[:-1] + z[1:]), [z[-1]]])

# Column 1: physics and nutrients  (individual colour scales)
col1 = [
    ('temp',          'Temperature',   '°C'),
    ('nuh',           'Diffusivity',   'm² s⁻¹'),
    ('salt',          'Salinity',      'PSU'),
    ('num_model_N',   'DIN',           'µg N L⁻¹'),
    ('num_model_Si',  'Silicate',      'µg Si L⁻¹'),
    ('num_model_DOC', 'DOC',           'µg C L⁻¹'),
]

# Column 2: biology (shared colour scale)
col2 = [
    ('num_model_BGen1',  'Generalists',        'µg C L⁻¹'),
    ('num_model_BDia1',  'Diatoms',            'µg C L⁻¹'),
    ('num_model_BPCop1', 'Passive copepods 1', 'µg C L⁻¹'),
    ('num_model_BPCop2', 'Passive copepods 2', 'µg C L⁻¹'),
    ('num_model_BACop1', 'Active copepods 1',  'µg C L⁻¹'),
    ('num_model_BACop2', 'Active copepods 2',  'µg C L⁻¹'),
    ('num_model_BACop3', 'Active copepods 3',  'µg C L⁻¹'),
]

# Load all needed variables
bio_names  = [v for v, *_ in col2]
phys_names = [v for v, *_ in col1 if v != 'nuh']
data = {v: np.array(ds.variables[v][:, :, 0, 0])
        for v in bio_names + phys_names}
ds.close()

# Shared log colour scale for all biology panels
bio_all = np.concatenate([data[v].ravel() for v in bio_names])
bio_pos = bio_all[bio_all > 0]
bio_norm = mcolors.LogNorm(vmin=np.nanpercentile(bio_pos, 2),
                           vmax=np.nanpercentile(bio_pos, 98))

CMAP = 'viridis'

# ---------------------------------------------------------------------------
# Plot 1: two-column space-time panels
# ---------------------------------------------------------------------------
nrows = max(len(col1), len(col2))   # 7

fig, axes = plt.subplots(nrows, 2, figsize=(11, 2.8 * nrows),
                         sharex=True, sharey=True)

def fmt_axes(ax, row, col):
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.tick_params(labelsize=7)
    if col == 0:
        ax.set_ylabel('Depth (m)', fontsize=8)

# --- Column 1 ---
for row, (varname, title, units) in enumerate(col1):
    ax = axes[row, 0]

    if varname == 'nuh':
        arr, zc = nuh, z_i
        norm = mcolors.LogNorm(vmin=max(nuh[nuh > 0].min(), 1e-5),
                               vmax=nuh.max())
    else:
        arr, zc = data[varname], z
        vmin = np.nanpercentile(arr, 2)
        vmax = np.nanpercentile(arr, 98)
        if vmax <= vmin:
            vmax = vmin + 1e-9
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    pcm = ax.pcolormesh(t_num, zc, arr.T, shading='nearest',
                        cmap=CMAP, norm=norm)
    cbar = fig.colorbar(pcm, ax=ax, pad=0.02, fraction=0.046)
    cbar.set_label(units, fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    ax.set_title(title, fontsize=9, fontweight='bold')
    fmt_axes(ax, row, 0)

# hide unused rows in column 1
for row in range(len(col1), nrows):
    axes[row, 0].set_visible(False)

# --- Column 2: shared colour scale ---
for row, (varname, title, units) in enumerate(col2):
    ax = axes[row, 1]
    arr = data[varname]
    pcm = ax.pcolormesh(t_num, z, arr.T, shading='nearest',
                        cmap=CMAP, norm=bio_norm)
    ax.set_title(title, fontsize=9, fontweight='bold')
    fmt_axes(ax, row, 1)

# Single shared colorbar for column 2, spanning all its rows
cbar2 = fig.colorbar(pcm, ax=axes[:len(col2), 1], pad=0.02, fraction=0.046)
cbar2.set_label('µg C L⁻¹', fontsize=7)
cbar2.ax.tick_params(labelsize=6)

fig.suptitle('NUMmodel / GOTM — OWS Papa (50°N, 145°W)', fontsize=11)
plt.tight_layout()

outfile = 'spacetime.png'
plt.savefig(outfile, dpi=150, bbox_inches='tight')
print(f'==> Saved {outfile}')

# ---------------------------------------------------------------------------
# Plot 2: depth-integrated biology time series
# ---------------------------------------------------------------------------
bio_labels = {
    'num_model_BGen1':  'Generalists',
    'num_model_BDia1':  'Diatoms',
    'num_model_BPCop1': 'Passive copepods 1',
    'num_model_BPCop2': 'Passive copepods 2',
    'num_model_BACop1': 'Active copepods 1',
    'num_model_BACop2': 'Active copepods 2',
    'num_model_BACop3': 'Active copepods 3',
}

colors = plt.get_cmap('tab10').colors

fig2, ax2 = plt.subplots(figsize=(12, 4))

for idx, (varname, label) in enumerate(bio_labels.items()):
    integrated = -np.trapezoid(data[varname], z, axis=1)
    ax2.plot(t_num, integrated, label=label,
             color=colors[idx % len(colors)], linewidth=1.4)

ax2.set_yscale('symlog', linthresh=1.0)
ax2.xaxis_date()
ax2.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax2.xaxis.set_major_locator(mdates.YearLocator())
ax2.tick_params(labelsize=8)
ax2.set_ylabel('Depth-integrated biomass (µg C L⁻¹ × m)', fontsize=9)
ax2.set_title('NUMmodel / GOTM OWS Papa — depth-integrated biomass', fontsize=10)
ax2.legend(fontsize=8, ncol=2, loc='upper left', framealpha=0.7)
plt.tight_layout()

outfile2 = 'integrated_biomass.png'
fig2.savefig(outfile2, dpi=150, bbox_inches='tight')
print(f'==> Saved {outfile2}')

plt.show()
