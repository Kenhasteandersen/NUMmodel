#!/usr/bin/env python3
"""
Generate synthetic meteorological forcing for the GOTM+NUMmodel North Sea example.

Produces meteo.dat with daily values for year 2000.
Uses only the Python standard library (no NumPy required).

Output columns (space-separated, GOTM format):
  YYYY-MM-DD HH:MM:SS  u10[m/s]  v10[m/s]  airt[°C]  hum[%]  cloud[-]  swr[W/m²]

Run from the gotm_watercolumn/ directory:
    python3 generate_forcing.py
"""

import math
import random
from datetime import datetime, timedelta

# --------------------------------------------------------------------------
# Site and simulation parameters
# --------------------------------------------------------------------------
LAT = 59.0           # degrees N
YEAR = 2000
CLOUD_MEAN = 0.65    # mean cloud fraction (typical North Sea)
SOLAR_CONST = 1368.0 # W/m²

lat_rad = math.radians(LAT)

ndays = 366 if (YEAR % 4 == 0 and (YEAR % 100 != 0 or YEAR % 400 == 0)) else 365
start = datetime(YEAR, 1, 1, 12, 0, 0)

rng = random.Random(42)   # fixed seed for reproducibility

# --------------------------------------------------------------------------
# Helper: clamp a value to [lo, hi]
# --------------------------------------------------------------------------
def clamp(x, lo, hi):
    return max(lo, min(hi, x))

# --------------------------------------------------------------------------
# Main loop: compute forcing for each day
# --------------------------------------------------------------------------
out_file = "meteo.dat"
swr_vals = []
airt_vals = []

with open(out_file, "w") as fh:
    fh.write("! Synthetic meteorological forcing — North Sea 59°N, year 2000\n")
    fh.write("! Columns: datetime  u10[m/s]  v10[m/s]  airt[°C]  "
             "hum[%]  cloud[-]  swr[W/m²]\n")

    for i in range(ndays):
        doy = i + 1   # 1-based day of year
        t   = start + timedelta(days=i)
        phi = 2.0 * math.pi

        # ------------------------------------------------------------------
        # Solar declination (radians)
        # ------------------------------------------------------------------
        decl = math.radians(23.45 * math.sin(phi * (doy - 81.0) / 365.0))

        # ------------------------------------------------------------------
        # Daily mean extra-terrestrial radiation (W/m²)
        # ------------------------------------------------------------------
        cos_hs = clamp(-math.tan(lat_rad) * math.tan(decl), -1.0, 1.0)
        hs = math.acos(cos_hs)   # sunset hour angle (radians)

        I0 = (SOLAR_CONST / math.pi) * (
            hs * math.sin(lat_rad) * math.sin(decl)
            + math.cos(lat_rad) * math.cos(decl) * math.sin(hs)
        )
        I0 = max(I0, 0.0)

        # ------------------------------------------------------------------
        # Cloud cover: larger in winter (North Sea)
        # ------------------------------------------------------------------
        cloud = CLOUD_MEAN + 0.15 * math.sin(phi * (doy - 180.0) / 365.0)
        cloud = clamp(cloud, 0.0, 1.0)

        # Cloud transmission (Atwater & Ball 1981)
        transmission = 1.0 - 0.75 * cloud**3.4
        swr = I0 * transmission

        # ------------------------------------------------------------------
        # Wind: westerly, slightly stronger in winter
        # ------------------------------------------------------------------
        u10 = 5.0 + 3.0 * math.cos(phi * (doy - 10.0) / 365.0)
        u10 += rng.gauss(0, 0.5)

        v10 = 0.5 * math.sin(phi * doy / 80.0)
        v10 += rng.gauss(0, 0.3)

        # ------------------------------------------------------------------
        # Air temperature: sinusoidal seasonal cycle (~4–17 °C)
        # ------------------------------------------------------------------
        airt = 10.5 + 6.5 * math.sin(phi * (doy - 55.0) / 365.0)

        # ------------------------------------------------------------------
        # Relative humidity (~76–88 %, higher in winter)
        # ------------------------------------------------------------------
        hum = 82.0 - 6.0 * math.sin(phi * (doy - 100.0) / 365.0)

        fh.write(
            f"{t.strftime('%Y-%m-%d %H:%M:%S')}  "
            f"{u10:7.3f}  {v10:7.3f}  "
            f"{airt:6.2f}  {hum:5.1f}  "
            f"{cloud:.3f}  {swr:7.2f}\n"
        )

        swr_vals.append(swr)
        airt_vals.append(airt)

print(f"Written {ndays} daily records to '{out_file}'")
print(f"  SWR range : {min(swr_vals):.1f} – {max(swr_vals):.1f} W/m²")
print(f"  Airt range: {min(airt_vals):.1f} – {max(airt_vals):.1f} °C")
