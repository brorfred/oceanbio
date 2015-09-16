
import pylab as pl

from njord import winds

class persist:
    pass

def seawinds(date, lonvec, latvec):

    jd = pl.datestr2num(date) if type(date) is str else date
    if not hasattr(persist, "swn"):
        persist.swn = winds.Seawinds()
        persist.swn.add_kd()
    ivec,jvec = persist.swn.ll2ij(lonvec, latvec)
    persist.swn.load(jd)
    return persist.swn.nwnd[jvec, ivec]
