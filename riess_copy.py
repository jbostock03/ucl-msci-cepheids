import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit

#load Riess et al. table
df = pd.read_csv("riess2018gaiadr2.csv", delimiter='  ', engine='python')

#data
#masking 4 w/ dodgy stats: T MON too close to saturation threshold; others too large uncertainties
#["RY VEL", "RW CAM", "SV PER", "T MON"]
mask = df["FinalAnalysisFlag"] == 0
df46 = df[mask]

logP = df46["LogPeriod"].values
mWH = df46["mWH"].values
parallax = df46["parallax"].values       # mas
parallax_err = df46["parallax_error"].values

#band magnitudes
m_f160w = df46["F160W"].values
m_f555w = df46["F555W"].values
m_f814w = df46["F814W"].values

#band magnitude errors
m_f160w_err = df46["F160W_error"].values
m_f555w_err = df46["F555W_error"].values
m_f814w_err = df46["F814W_error"].values

# count-rate non-linearity correction
crnl = 0.052
crnl_err = 0#.014    # excluded in Riess et al. ***

# apparent Wesenheit magnitude
mWH_calc = m_f160w - 0.386*(m_f555w - m_f814w) + crnl    

mWH_err = np.sqrt(m_f160w_err**2 + (0.386*m_f555w_err)**2 + (0.386*m_f814w_err)**2 + crnl_err**2)
    
#Riess et al. calibration
# M_W = a (logP − 1) + b

def PL(logP, m, c):
    return m * (logP - 1.0) + c

#corrected parallaxes
zero_point = -0.046 #mas
zp_err = 0.006      #mas

parallax_corr = parallax - zero_point
parallax_corr_err = np.sqrt(parallax_err**2 + zp_err**2)

#absolute Wesenheit magnitude
mu = 5*np.log10(1e3/parallax_corr) - 5
MWH = mWH_calc - mu

#error in abs M_WH
MWH_err = np.sqrt( (5/(np.log(10)) * (parallax_corr_err / parallax_corr))**2 + mWH_err**2 )

#least squares
params, cov = curve_fit(
    PL,
    logP,
    MWH,
    sigma=MWH_err,
)

m, c = params
m_err, c_err = np.sqrt(np.diag(cov))

print(f"M^W_H = ({m:.3f} +- {m_err:.3f}) * (logP - 1) + ({c:.3f} +- {c_err:.3f})")

#plot
plt.errorbar(logP, MWH, yerr=MWH_err, fmt='.', alpha=0.7, capsize=3)
x = np.linspace(min(logP), max(logP), 200)
plt.plot(x, PL(x, m, c))

plt.gca().invert_yaxis()
plt.xlabel(r"log $P$ (days)")
plt.ylabel(r"$M_W$ (Absolute Wesenheit magnitude)")
plt.title("Milky Way Cepheid PL Relation (Gaia DR2)")
plt.show()