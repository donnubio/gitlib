import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import FloatVector
import os
#import plotly.express as px
#import plotly.graph_objects as go


def get_hrvprms(rr, fd_interp_rri=4, td_win_size=60, fd_stft_win_size=60, fd_stft_win_shift=5):

    file_path = os.path.dirname(os.path.realpath(__file__))
    ro.r['source'](file_path + '/hrv.R')
    fun_get_hrvprms_r = ro.globalenv['get_hrvprms'] #rr, fd_interp_rri=4, td_win_size=60, fd_stft_win_size=60, fd_stft_win_shift=5
    rr_r = FloatVector( rr )
    hrv_r = fun_get_hrvprms_r(rr_r, fd_interp_rri, td_win_size, fd_stft_win_size, fd_stft_win_shift)

    hrvp = { key : np.array( hrv_r.rx(key)[0] ) for key in hrv_r.names }
    for n in hrvp.keys():
        if hrvp[n].size==1:
            hrvp[n]=hrvp[n][0]

    return(hrvp)


    return(D,A, ba_pr, ba_pr_L, ba_pr_R)


 