import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy.optimize import curve_fit
from zero_point import zpt

def PL_relation_ls_fit(dr=2, mask_binaries=True, color=None, label=None):
    """
    Plots PL relation from least squares using either Gaia DR2 or DR3 data

    input:
    dr=2: Data Release 2
    dr=3: Data Release 3
    color & label: plot characteristics
    """

    if dr == 2 or dr == 3:
        df = pd.read_csv("xmatch_riess_col_corr_gaia_source-result.csv")

    else:
        raise ValueError('dr argument must equal 2 or 3!')

    #data
    #masking 4 w/ dodgy stats: T MON too close to saturation threshold; others too large uncertainties
    #["RY VEL", "RW CAM", "SV PER", "T MON"]
    mask = df["FinalAnalysisFlag"] == 0
    
    #NEW BINARIES: ['FO CAR', 'MY PUP', 'R CRU', 'VY PER', 'VZ PUP']
    #OLD BINARIES: ['DL CAS', 'FF AQL', 'S MUS', 'S SGE', 'SY NOR', 'TX MON', 'U VUL', 'V1334 CYG', 'YZ CAR']
    #EXTRA 2: ['X PUP', 'XX SAG']
    binaries = ['FO CAR', 'MY PUP', 'R CRU', 'VY PER', 'VZ PUP', 'DL CAS', 'FF AQL', 'S MUS', 'S SGE', 'SY NOR', 'TX MON', 'U VUL', 'V1334 CYG', 'YZ CAR', 'X PUP', 'XX SAG']
    
    df46 = df[mask].reset_index(drop=True)
    if mask_binaries == True:
        df46 = df46[~df46['Name'].isin(binaries)]
    logP = df46["LogPeriod"].values
    mWH = df46["mWH"].values

    if dr == 2:
        parallax = df46["parallax"].values       # mas
        parallax_err = df46["parallax_error"].values

        #zero point offset
        zero_point = -0.046 #mas
        zp_err = 0.006      #mas

    #for using dr3
    if dr == 3:
        parallax = df46["parallax.1"].values       # mas
        parallax_err = df46["parallax_error.1"].values

        #zp using gaiadr3_zeropoint, Lindegren+20. in mas
        zpt.load_tables()
        zero_point_nan = zpt.get_zpt(
            df46['phot_g_mean_mag.1'], 
            df46['nu_eff_used_in_astrometry'],
            df46['pseudocolour'],
            df46['ecl_lat'],
            df46['astrometric_params_solved']
        )
        #removing nan values in zp
        zero_point = np.nan_to_num(zero_point_nan, nan=np.nanmean(zero_point_nan))
        zp_err = 0.0 #mas

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
    crnl_err = 0.014

    # apparent Wesenheit magnitude
    mWH_calc = m_f160w - 0.386*(m_f555w - m_f814w) + crnl    

    mWH_err = np.sqrt(m_f160w_err**2 + (0.386*m_f555w_err)**2 + (0.386*m_f814w_err)**2 + crnl_err**2)
    
    #Riess et al. calibration
    # M_W = m(logP − 1) + c

    def PL(logP, m, c):
        return m * (logP - 1.0) + c

    #corrected parallaxes
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
    x = np.linspace(min(logP), max(logP), 200)
    if color == None:
        plt.errorbar(logP, MWH, yerr=MWH_err, fmt='.', alpha=0.7, capsize=3)
        plt.plot(x, PL(x, m, c), label=label)
    else:
        plt.errorbar(logP, MWH, yerr=MWH_err, fmt='.', alpha=0.7, capsize=3, color=color)
        plt.plot(x, PL(x, m, c), label=label)

    plt.gca().invert_yaxis()
    plt.xlabel("log P (days)")
    plt.ylabel("M_W (Absolute Wesenheit magnitude)")
    plt.title(f"Milky Way Cepheid PL Relation (Gaia DR{dr})")
    plt.grid()
    if not label == None:
        plt.legend()
    plt.show()

    return m, c, logP, MWH, MWH_err
        