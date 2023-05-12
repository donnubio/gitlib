import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import FloatVector
import os
#import plotly.express as px
#import plotly.graph_objects as go


def GetHRVIndices(rr, 
                  fd_interp_rri=4, 
                  td_win_size=60, 
                  fd_stft_win_size=60, 
                  fd_stft_win_shift=5, 
                  fd_stft_win_extend=None,
                  fd_time_mean=True,
                  ULF_L = 0, ULF_R = 0.03,
                  VLF_L = 0.03, VLF_R = 0.05,
                  LF_L = 0.05, LF_R = 0.15,
                  HF_L = 0.15, HF_R = 0.4,
                  wavelet = "d4",
                  bandtolerance = 0.01):
    
    file_path = os.path.dirname(os.path.realpath(__file__))
    ro.r['source'](file_path + '/hrv.R')
    fun_get_hrvprms_r = ro.globalenv['get_hrvprms'] #rr, fd_interp_rri=4, td_win_size=60, fd_stft_win_size=60, fd_stft_win_shift=5
    rr_r = FloatVector( rr )
    if fd_stft_win_extend is None:
        fd_stft_win_extend = ro.r("NULL")
    hrv_r = fun_get_hrvprms_r(rr_r, 
                              fd_interp_rri, td_win_size, 
                              fd_stft_win_size, fd_stft_win_shift, fd_stft_win_extend, fd_time_mean, 
                              ULF_L, ULF_R,
                              VLF_L, VLF_R,
                              LF_L, LF_R,
                              HF_L, HF_R,
                              wavelet,
                              bandtolerance)

    hrvp = { key : np.array( hrv_r.rx(key)[0] ) for key in hrv_r.names }
    for n in hrvp.keys():
        if hrvp[n].size==1:
            hrvp[n]=hrvp[n][0]

    return(hrvp)

def SlidingHRV(t_NN, NN, 
               win_segm, 
               win_ovlp,  
               fd_interp_rri=4, 
               td_win_size=60, 
               fd_stft_win_size=60, 
               fd_stft_win_shift=5, 
               fd_time_mean=True,
               ULF_L = 0, ULF_R = 0.03,
               VLF_L = 0.03, VLF_R = 0.05,
               LF_L = 0.05, LF_R = 0.15,
               HF_L = 0.15, HF_R = 0.4,
               wavelet = "d4",
               bandtolerance = 0.01):
    file_path = os.path.dirname(os.path.realpath(__file__))
    ro.r['source'](file_path + '/hrv.R')
    SlidingHRV_r = ro.globalenv['SlidingHRV'] #rr, fd_interp_rri=4, td_win_size=60, fd_stft_win_size=60, fd_stft_win_shift=5
    t_NN_r = FloatVector( t_NN )
    NN_r = FloatVector( NN )
    hrv_r = SlidingHRV_r( t_NN_r, NN_r, 
                            win_segm, win_ovlp,  
                            fd_interp_rri, 
                            td_win_size, 
                            fd_stft_win_size, 
                            fd_stft_win_shift,
                            ULF_L, ULF_R,
                            VLF_L, VLF_R,
                            LF_L, LF_R,
                            HF_L, HF_R,
                            wavelet,
                            bandtolerance)   
    hrvp = { key : np.array( hrv_r.rx(key)[0] ) for key in hrv_r.names }
    for n in hrvp.keys():
        if hrvp[n].size==1:
            hrvp[n]=hrvp[n][0]

    return(hrvp)
  