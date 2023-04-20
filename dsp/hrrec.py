from .dsp import SignalCutSegment
import numpy as np
from scipy.optimize import curve_fit
from scipy.interpolate import splrep, BSpline, splev, spalde, sproot

###################################################################
##################Logistic curve###################################
###################################################################

def InitialGuessForLogisticFunParams(x, y):
    nknock = len(x)*15
    xs = np.arange(x[0], x[-1], 0.1)
    s = splrep(x, y, s=nknock)
    ys = BSpline(*s)(xs)
    sd = spalde(xs, s)
    sd1 = [d[1] for d in sd] #first deriv
    ymax = splev(x[0], s)
    ymin = splev(x[-1], s)
    #ymdl = (ymax+ymin)/2
    #xmdl = sproot()
    xmdl=xs[ np.argmin( sd1 ) ]
    return(ymax,ymin,xmdl)

def LogisticFun(x, ymin, ymax, xmdl, tau):
    return ymin + (ymax - ymin) / ( 1 + np.exp( tau *( x - xmdl ) ) ) 

def FitLogisticFun(tHR, HR, t1=None, t2=None):
    tHR, HR = SignalCutSegment(tHR, HR, t1, t2, time_from_zero=True)
    (HRmax0,HRmin0,tHRmdl0) = InitialGuessForLogisticFunParams(tHR, HR)
    tau0 = 0.1
    p0 = [HRmin0, HRmax0, tHRmdl0, tau0]
    bounds=(0., [300, 400, 300, 10])
    params_fit, covr = curve_fit(LogisticFun, tHR, HR, p0, bounds=bounds)
    HR_fit = LogisticFun(tHR, *params_fit)
    #HR_fit = infodict['fvec']
    return(params_fit, HR_fit, covr, tHR, HR)