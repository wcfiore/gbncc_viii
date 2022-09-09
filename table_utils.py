import astropy.constants as const
# See pygedm (https://github.com/FRBs/pygedm)
# > conda install -c conda-forge f2c
# > pip install pygedm
import pygedm
import uncertainties as unc
import astropy.units as u
import numpy as np


def uf(param,digits=1):
    """Simple uncertainty formatting
    
    Input e.g. model param: model['ELONG']
    Output e.g. string: 65.90834(5)
    """
    val = param.value
    err = param.uncertainty_value
    x = unc.ufloat(val, err)
    
    fstr = "{:.%duSL}" % digits
    return fstr.format(x)

def ufve(val,err,digits=1):
    x = unc.ufloat(val, err)
    
    fstr = "{:.%duSL}" % digits
    return fstr.format(x)

def get_dmdists(gal_coord,dm):
    """ Get YMW17/NE2001 DM distances (kpc)
    
    Input: astropy galactic coord object; dm quantity
    Output: DM distance list NE/YWM (kpc)
    """
    dmdist_ne, tau_sc = pygedm.dm_to_dist(gal_coord.l,gal_coord.b,dm.value,method='ne2001')
    dmdist_ymw, tau_sc = pygedm.dm_to_dist(gal_coord.l,gal_coord.b,dm.value,method='ymw16')
    return dmdist_ne.to(u.kpc).value, dmdist_ymw.to(u.kpc).value

def format_ra(eqcoord_obj):
    """ LaTex-format RA
    E.g.: $20^{\rm h}\, 22^{\rm m}\, 33\, \fs256$
    Input: astropy equatorial coord object
    """
    hms = eqcoord_obj.ra.hms
    s = f"{hms.s:.2f}"
    # Really janky, but deals with s1 0 padding
    s1 = s.split('.')[0]
    s2 = s.split('.')[1]
    hms_form = f"${int(hms.h):02}^{{\\rm h}}\\, {int(hms.m):02}^{{\\rm m}}\\, {int(s1):02}\\, \\fs{s2}$"
    return(hms_form)

def format_dec(eqcoord_obj):
    """ LaTex-format Dec
    E.g.: $+25\arcdeg\, 34\arcmin\, 42\, \farcs5$
    Input: astropy equatorial coord object
    """
    dms = eqcoord_obj.dec.dms
    # {:+} adds a sign regardless of -/+
    # abs() to remove '-' from Dec m, s
    s = f"{abs(dms.s):.1f}"
    # Really janky, but deals with s1 0 padding
    s1 = s.split('.')[0]
    s2 = s.split('.')[1]
    dms_form = f"${int(dms.d):+03}\\arcdeg\\, {int(abs(dms.m)):02}\\arcmin\\, {int(s1):02}\\, \\farcs{s2}$"
    return(dms_form)

def PMtot_err(PMx,PMy):
    """Calculate total proper motion and uncertainty
    
    Input: PMx, PMy (PINT model quantities)
    """
    PMx_err, PMy_err = (PMx.uncertainty_value,PMy.uncertainty_value)
    pmtot = np.sqrt(PMx.value**2+PMy.value**2)
    pmtot_err = np.sqrt(PMx.value**2*PMx_err**2+PMy.value**2*PMy_err**2)/pmtot
    return (pmtot, pmtot_err) * (u.mas / u.yr)

def Vtrans_err(PMtot,PMtot_err,D,D_err):
    """Calculate transverse velocity and uncertainty
    
    Input: PMtot, PMtot error, distance, and D_err (all quantities)
    Output: VT, VTerr (quantities; km/s)
    """
    with u.set_enabled_equivalencies(u.dimensionless_angles()): # Fancy trick from pint/derived_quantities.py
        vtrans = (PMtot*D).to(u.km/u.s)
        vtrans_err = np.sqrt((D*PMtot_err)**2+(PMtot*D_err)**2).to(u.km/u.s)
    return vtrans,vtrans_err

def pd_gal(P0,gal_coords,D):
    """Calculate P1 contribution due to Galactic potential
    
    Based on Guo et al (2021); Section 7
    https://ui.adsabs.harvard.edu/abs/2021A%26A...654A..16G/abstract
    
    Inputs: spin period quantity, astropy galactic coord object, distance quantity
    """
    R0 = 8.275*u.kpc # Sun to GC distance 8.275(34) kpc; Gravity Collaboration (2021)
    phi0 = 240.5*u.km/u.s # Galactic circular velocity at Sun's location 240.5(41) km/s; GC (2021)
    gl_rad,gb_rad = (gal_coords.l.rad,gal_coords.b.rad)
    beta = (D/R0)*np.cos(gb_rad)-np.cos(gl_rad)
    z = np.abs(D*np.sin(gb_rad)).to(u.kpc)
    Kz = (2.27*z.value+3.68*(1-np.exp(-4.3*z.value)))*1e-9*u.cm/u.s**2 # Only for z <= 1.5 kpc??
    gal_vert = -1.0*Kz*np.abs(np.sin(gb_rad))/const.c
    gal_plane = -1.0*phi0**2/(const.c*R0)*(np.cos(gl_rad)+beta/(np.sin(gl_rad)**2+beta**2))*np.cos(gb_rad)

    pd_gal = (gal_vert+gal_plane)*P0.decompose()
    return pd_gal
