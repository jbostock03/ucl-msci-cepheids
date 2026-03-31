import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
import seaborn as sns

sns.set_theme(style="whitegrid")

# load Gaia table
df = pd.read_csv("22jan_full_cepheids.csv")

# mask negative parallaxes
mask = df['parallax'] > 0.1
df = df.loc[mask].copy()

# distance
df['distance_pc'] = 1000.0 / df['parallax']

# convert to galactic coordinates
coords = SkyCoord(
    ra=df['ra'].values * u.deg,
    dec=df['dec'].values * u.deg,
    distance=df['distance_pc'].values * u.pc,
    frame='icrs'
)

gal = coords.galactic

df['l_deg'] = gal.l.wrap_at(180 * u.deg).deg
df['b_deg'] = gal.b.deg

# Z height
df['Z_pc'] = df['distance_pc'] * np.sin(np.deg2rad(df['b_deg']))

# subtract Sun offset
sun_z = 20.8
df['Z_pc_nosun'] = df['Z_pc'] - sun_z

df['above'] = df['Z_pc_nosun'] > 0

# l vs b kde coloured by Z
plt.figure(figsize=(8,6))

sns.scatterplot(
    data=df,
    x="l_deg",
    y="b_deg",
    s=7,
    alpha=0.95,
    color="gray"
)

sns.kdeplot(
    data=df,
    x="l_deg",
    y="b_deg",
    levels=20,
    color="black",
    bw_adjust=0.7,
    alpha=0.2
)

plt.axhline(0, color='k', lw=0.8)
plt.xlim(-190,190)
plt.xlabel(r"Galactic longitude, $l$ (deg)")
plt.ylabel(r"Galactic latitude, $b$ (deg)")
plt.gca().invert_xaxis()

plt.tight_layout()
plt.savefig("kde_lb.pdf")
