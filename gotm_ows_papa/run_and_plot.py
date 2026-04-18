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

# Column 1: physics and nutrients (individual colour scales)
col1 = [
    ('temp',          'Temperature',   '°C'),
    ('nuh',           'Diffusivity',   'm² s⁻¹'),
    ('salt',          'Salinity',      'PSU'),
    ('num_model_N',    'DIN',           'µg N L⁻¹'),
    ('num_model_Si',   'Silicate',      'µg Si L⁻¹'),
    ('num_model_DOC',  'DOC',           'µg C L⁻¹'),
    ('num_model_BPOM1','POM',           'µg C L⁻¹'),
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

bio_names  = [v for v, *_ in col2] + ['num_model_BPOM1']
phys_names = [v for v, *_ in col1 if v not in ('nuh', 'num_model_BPOM1')]
prod_names = ['num_model_ProdGross', 'num_model_ProdNet', 'num_model_ProdHTL']
data = {v: np.array(ds.variables[v][:, :, 0, 0])
        for v in bio_names + phys_names + prod_names}
ds.close()

# Shared log colour scale for all biology panels
bio_all = np.concatenate([data[v].ravel() for v in bio_names])
bio_pos = bio_all[bio_all > 0]
bio_norm = mcolors.LogNorm(vmin=1e-5,
                           vmax=np.nanpercentile(bio_pos, 98))

CMAP  = 'viridis'
TITLE_FS = 14
LABEL_FS = 13
TICK_FS  = 11

# ---------------------------------------------------------------------------
# Plot 1: two-column space-time panels
# ---------------------------------------------------------------------------
nrows = max(len(col1), len(col2))   # 7

fig, axes = plt.subplots(nrows, 2, figsize=(12, 3.0 * nrows),
                         sharex=True, sharey=True)

# Reserve space on the right for the shared col-2 colorbar
fig.subplots_adjust(left=0.07, right=0.86, top=0.96, bottom=0.04,
                    hspace=0.32, wspace=0.08)

def add_xlabels(ax):
    ax.xaxis_date()
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
    ax.xaxis.set_major_locator(mdates.YearLocator())
    ax.tick_params(axis='x', labelsize=TICK_FS, labelbottom=True)

# --- Column 1: individual colour scales ---
for row, (varname, title, units) in enumerate(col1):
    ax = axes[row, 0]

    if varname == 'nuh':
        arr, zc = nuh, z_i
        norm = mcolors.LogNorm(vmin=1e-5, vmax=nuh.max())
    elif varname == 'num_model_BPOM1':
        arr, zc = data[varname], z
        norm = bio_norm   # shared scale with biology column
    else:
        arr, zc = data[varname], z
        vmin = np.nanpercentile(arr, 2)
        vmax = np.nanpercentile(arr, 98)
        if vmax <= vmin:
            vmax = vmin + 1e-9
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    pcm = ax.pcolormesh(t_num, zc, arr.T, shading='nearest',
                        cmap=CMAP, norm=norm)
    if varname != 'num_model_BPOM1':
        cbar = fig.colorbar(pcm, ax=ax, pad=0.03, fraction=0.05)
        cbar.set_label(units, fontsize=LABEL_FS)
        cbar.ax.tick_params(labelsize=TICK_FS)
    ax.set_title(title, fontsize=TITLE_FS, fontweight='bold')
    ax.set_ylabel('Depth (m)', fontsize=LABEL_FS)
    ax.tick_params(axis='y', labelsize=TICK_FS)

# x-axis labels only on last visible row of col 1 (DOC) and bottom of col 2
add_xlabels(axes[len(col1) - 1, 0])   # POM panel (last in col 1)

# Hide unused rows in column 1
for row in range(len(col1), nrows):
    axes[row, 0].set_visible(False)

# --- Column 2: shared colour scale ---
for row, (varname, title, _) in enumerate(col2):
    ax = axes[row, 1]
    pcm2 = ax.pcolormesh(t_num, z, data[varname].T, shading='nearest',
                         cmap=CMAP, norm=bio_norm)
    ax.set_title(title, fontsize=TITLE_FS, fontweight='bold')
    ax.tick_params(axis='y', labelsize=TICK_FS)

add_xlabels(axes[len(col2) - 1, 1])   # Active copepods 3 panel

# Single shared colorbar for column 2, placed to the right of the figure
fig.canvas.draw()                       # fix axes positions before reading them
pos_top = axes[0, 1].get_position()
pos_bot = axes[len(col2) - 1, 1].get_position()
cbar_ax = fig.add_axes([0.88, pos_bot.y0, 0.02,
                         pos_top.y1 - pos_bot.y0])
sm = plt.cm.ScalarMappable(cmap=CMAP, norm=bio_norm)
cbar2 = fig.colorbar(sm, cax=cbar_ax)
cbar2.set_label('µg C L⁻¹', fontsize=LABEL_FS)
cbar2.ax.tick_params(labelsize=TICK_FS)

# Y-axis: surface at top (0 m), bottom at depth
axes[0, 0].set_ylim(z.min(), 0)

fig.suptitle('NUMmodel / GOTM — OWS Papa (50°N, 145°W)',
             fontsize=TITLE_FS + 2, y=0.99)

outfile = 'spacetime.png'
fig.savefig(outfile, dpi=150, bbox_inches='tight')
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

prod_labels = {
    'num_model_ProdGross': 'Gross production',
    'num_model_ProdNet':   'Net production',
    'num_model_ProdHTL':   'HTL production',
}

colors = plt.get_cmap('tab10').colors

fig2, (ax2, ax3) = plt.subplots(2, 1, figsize=(12, 8), sharex=True)

for idx, (varname, label) in enumerate(bio_labels.items()):
    integrated = -np.trapezoid(data[varname], z, axis=1)
    ax2.plot(t_num, integrated, label=label,
             color=colors[idx % len(colors)], linewidth=1.6)

ax2.set_yscale('symlog', linthresh=1.0)
ax2.tick_params(labelsize=TICK_FS + 1)
ax2.set_ylabel('Depth-integrated biomass\n(µg C L⁻¹ × m)', fontsize=LABEL_FS + 1)
ax2.set_title('NUMmodel / GOTM OWS Papa — depth-integrated biomass',
              fontsize=TITLE_FS + 1)
ax2.legend(fontsize=TICK_FS + 1, ncol=2, loc='upper left', framealpha=0.7)

for idx, (varname, label) in enumerate(prod_labels.items()):
    # Production is in mg C m⁻³ d⁻¹; integrate over depth → mg C m⁻² d⁻¹
    integrated = -np.trapezoid(data[varname], z, axis=1)
    ax3.plot(t_num, integrated, label=label,
             color=colors[idx % len(colors)], linewidth=1.6)

ax3.xaxis_date()
ax3.xaxis.set_major_formatter(mdates.DateFormatter('%Y'))
ax3.xaxis.set_major_locator(mdates.YearLocator())
ax3.tick_params(labelsize=TICK_FS + 1)
ax3.set_ylabel('Depth-integrated production\n(mg C m⁻² d⁻¹)', fontsize=LABEL_FS + 1)
ax3.set_title('Depth-integrated production rates', fontsize=TITLE_FS + 1)
ax3.legend(fontsize=TICK_FS + 1, loc='upper left', framealpha=0.7)

plt.tight_layout()

outfile2 = 'integrated_biomass.png'
fig2.savefig(outfile2, dpi=150, bbox_inches='tight')
print(f'==> Saved {outfile2}')

plt.show()
