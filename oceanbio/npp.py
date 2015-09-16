"""Functions to calculate growth rates or  NPP from properties that can
   be derived from climatologies or satellite fields.

   Most function will return calculated growth rates, that can be converted
   to NPP with biomass (chl or POC): NPP = mu * biomass
"""

import numpy as np

def arrigo(temp, irr):
    """Calculate growthrates using the Arrigo algorithm

    Growthrates estimated from temperature and PAR using an algorithm
    optimized for low-temorature conditions in the Southern Ocean.
    Combine the restulting values with estimates of biomass (e.g. Chl)
    to get Net Primary Production. 

    temp: Temperature (deg C)
     irr: Representative irridiance (mol quanta m-2 d-1) 

     Arrigo, K.R., Worthen, D., Schnell, A. & Lizotte, M., 1998.
         Primary production in Southern Ocean waters. JGR-Oceans
         DOI:10.1029/98JC00930
    """
    mu_max  = 0.59 * np.exp(temp * 0.0633)
    E_k_max = 80.
    B = np.exp(10.089 - 2.12 * np.log(E_k_max))
    E_k = E_k_max / (1 + 10 * np.exp(-B * irr * 1e6 / (12*60*60) ) )
    return mu_max * (1 - np.exp( -(irr * 1e6 / (12*60*60) ) / E_k)) 

