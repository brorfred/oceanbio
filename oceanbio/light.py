
import numpy as np

def depthmean(irr0, k490=None, eup=None, zlev=None):
    """Calculate MEAN irradiance down to zlev.

    Calculate the MEAN light field in a water body. Depth-resolved
    light attenuations is estimated from two datapoints - surface light (irr0)
    and the depth where on 1% of the light remains. The lower datapoints
    can be given either as k490/kPAR or the 'euphotic depth' (one
    percent light level).

    irr0: Surface Irradiance (mol quanta m-2 d-1)
    k490: Diffuse attenuation coefficient at 490 nm (m-1)
     eup: One percen light level (euphotic depth) (m)
    zlev: Depth to avereage irr to (m)

    return
     irr: Mean light down to zlev (mol quanta m-2 d-1)


    Hood et al: Remote estimation of nitrogen fixation by Trichodesmium.
        Deep Sea Research Part II - Topical Studies in Oceanography (2002) 
        vol. 49 (1-3) pp. 123-147
    Cloern et al. An empirical model of the phytoplankton chlorophyll:carbon
        ratio - The conversion factor between productivity and growth rate.
        Limnology and Oceanography (1995) vol. 40 (7) pp. 1313-1321     
    """
    zmax = np.nanmax(zlev)
    if (k490 is not None) & (eup is None):
        nanmask = ~np.isnan(zlev + k490 + irr0)
        kpar = 0.0304 + 0.893 * k490    #Eq 10, Hood 2002
    elif (k490 is None) & (eup is not None):
        nanmask = ~np.isnan(zlev + eup + irr0)
        kpar = (4.6/eup)
    elif (k490 is None) & (eup is None):
        raise ValueError, "Either k490 or eup has to be set"
    else:
        raise ValueError, "Both k490 or eup can't be set"

    try:
        len(irr0)
    except TypeError:
        irr0 = np.array([irr0,])
        kpar = np.array([kpar,])
        zlev = np.array([zlev,])
    if not (len(irr0) == len(kpar) == len(zlev)):
        raise IndexError, "All inport values must have the same length"
    
            
    zmat = np.zeros((len(irr0), zmax+1)) * np.nan
    for z in np.arange(zmax): zmat[z<=zlev, z] = z

    irrmat = irr0[:,None] * np.exp(-kpar[:,None] * zmat ) #Eq 11, Cloern 1995
    return np.nanmean(irrmat, axis=1)


def depthmedian(irr0, k490=None, eup=None, zlev=2.):
    """Calculate MEDIAN irradiance down to zlev.

    Calculate the MEDIAN light field in a water body. Depth-resolved
    light attenuations is estimated from two datapoints - surface light (irr0)
    and the depth where on 1% of the light remains. The lower datapoints
    can be given either as k490/kPAR or the 'euphotic depth' (one
    percent light level).

    irr0: Surface Irradiance (mol quanta m-2 d-1)
    k490: Diffuse attenuation coefficient at 490 nm (m-1)
     eup: One percent light level (euphotic depth) (m)
    zlev: Depth to avereage irr to (m)

    return
     irr: Median light down to zlev (mol quanta m-2 d-1)
     
    Hood et al: Remote estimation of nitrogen fixation by Trichodesmium.
        Deep Sea Research Part II - Topical Studies in Oceanography (2002) 
        vol. 49 (1-3) pp. 123-147
    Cloern et al. An empirical model of the phytoplankton
        chlorophyll:carbon ratio - The conversion factor between producti-
        vity and growth rate. Limnology and Oceanography (1995)
        vol. 40 (7) pp. 1313-1321     
    """
    if (eup is None) & (k490 is not None):
        kpar = 0.0304 + 0.893 * k490    #Eq 10, Hood 2002
        eup = 4.6/k490
    elif (eup is not None) & (k490 is None):
        kpar = (4.6/eup)
    elif (eup is None) & (k490 is None):
        raise KeyError, "Neither k490 or eup is provided"
    else:
        kpar = 0.0304 + 0.893 * k490    #Eq 10, Hood 2002
    return irr0 * np.exp(-kpar * (eup/zlev) )    #Eq 11, Cloern 1995

def eup_morel(chl):
    """ Calculate euphotic depth from chl with Morel's Case I model"""
    chl_tot = 40.2 * chl**0.507
    mask = chl < 1
    chl_tot[mask] = 38.0 * chl[mask]**0.425
    eup = 200.0 * chl_tot**-0.293
    mask = eup <= 102.0
    eup[mask] = 568.2 * chl_tot[mask]**-0.746
    return eup

