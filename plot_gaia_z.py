import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import seaborn as sns

#load Gaia table
df = pd.read_csv("22jan_full_cepheids.csv")   # need ra, dec, parallax, parallax_error (parallax in mas)

#mask negative parallaxes
mask = df['parallax'] > 0.1
df = df.loc[mask].copy()

#distance
df['distance_pc'] = 1000.0 / df['parallax']

#convertign to to galactic lat b and calculating Z
coords = SkyCoord(ra=df['ra'].values*u.deg, dec=df['dec'].values*u.deg, distance=df['distance_pc'].values*u.pc, frame='icrs')
gal = coords.galactic
df['b_deg'] = gal.b.deg
df['Z_pc'] = df['distance_pc'] * np.sin(np.deg2rad(df['b_deg']))

#subtract Sun offset
sun_z = 20.8  #pc
df['Z_pc_nosun'] = df['Z_pc'] - sun_z

#scatter coloured by above/below
df['above'] = df['Z_pc_nosun'] > 0

#plt.figure(figsize=(8,6))
# scatter: projected distance on sky (longitude as x)
#plt.scatter(gal.l.wrap_at(180*u.deg).deg, df['Z_pc_nosun'], s=8, c=df['Z_pc_nosun'], cmap='viridis', alpha=0.7)
#plt.axhline(0, color='k', lw=0.8)
#plt.xlabel("Galactic longitude (deg)")
#plt.ylabel("Z (pc) — height relative to midplane")
#plt.title("Gaia objects: Z vs Galactic longitude")
#plt.colorbar(label='Z (pc)')
#plt.grid()
#plt.gca().invert_xaxis()
#plt.tight_layout()

#histogram of Z (above vs below)
plt.figure(figsize=(6,4))
plt.hist(df['Z_pc_nosun'], bins='fd')
plt.xlabel('Z (pc)')
plt.ylabel('Number of sources')
plt.title('Distribution of Z')
plt.grid()
plt.tight_layout()

plt.figure(figsize=(8,6))
# scatter: projected distance on sky (longitude as x)
plt.scatter(gal.l.wrap_at(180*u.deg).deg, gal.b.deg, s=8, c=df['Z_pc_nosun'], cmap='viridis', alpha=0.7)
plt.axhline(0, color='k', lw=0.8)
plt.xlabel(r"Galactic longitude, $l$ (deg)")
plt.ylabel(r"Galactic latitude, $b$ (deg)")
#plt.title("Gaia objects: b vs l, latitude v. longitude")
plt.colorbar(label='Distance to midplane (pc)')
plt.gca().invert_xaxis()
plt.tight_layout()

plt.savefig('plot_gaia_z.png')
