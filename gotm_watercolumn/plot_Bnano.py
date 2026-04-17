#!/usr/bin/env python3
"""Plot num_model_Bnanoerature as a function of depth and time."""

import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from datetime import datetime, timedelta

ds = nc.Dataset('NUMmodel_watercolumn.nc')

t_secs = ds.variables['time'][:]
t_dates = [datetime(2000, 1, 1) + timedelta(seconds=float(s)) for s in t_secs]

z = np.array(ds.variables['z'][0, :, 0, 0])
T = np.array(ds.variables['num_model_Bnano'][:, :, 0, 0])

ds.close()

t_num = mdates.date2num(t_dates)

fig, ax = plt.subplots(figsize=(11, 5))

pcm = ax.pcolormesh(t_num, z, T.T, shading='nearest', cmap='YlGn')

cbar = fig.colorbar(pcm, ax=ax, pad=0.02)
cbar.set_label('Bnano biomass  (mg C m⁻³)', fontsize=11)

ax.set_ylabel('Depth  (m)', fontsize=11)
ax.set_title('Bnano biomass — NUMmodel / GOTM  (59°N, year 2000)', fontsize=12)
ax.xaxis_date()
ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
ax.xaxis.set_major_locator(mdates.MonthLocator())
fig.autofmt_xdate(rotation=0, ha='center')

plt.tight_layout()
plt.savefig('num_model_Bnano_depth_time.png', dpi=150)
print('Saved num_model_Bnano_depth_time.png')
plt.show()
