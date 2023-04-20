import numpy as np

import scipy
import scipy.io as sio
from scipy import signal
from scipy.optimize import curve_fit
from scipy.interpolate import splrep, BSpline, splev, spalde, sproot

import bokeh.plotting as bh
from bokeh.io import output_notebook, save
from bokeh.palettes import Dark2_5 as palette
import itertools
output_notebook()

import plotly.express as px
import plotly.graph_objects as go
import plotly.offline as pyo
#from plotly.offline import init_notebook_mode #to plot in jupyter notebook
#init_notebook_mode() # init plotly in jupyter notebook

import matplotlib.pyplot as plt

##interactive matplotlib in jupyter:
#%matplotlib widget
##in case of GoogleCollab:
#!pip install ipympl -q
#from google.colab import output
#output.enable_custom_widget_manager()

import heartpy as hp
from heartpy.datautils import rolling_mean

import neurokit2 as nk

import hrvanalysis as hrva
#pip install hrv-analysis -q

import sys

####################################################################  

def LoadData(filename):
    m = sio.loadmat(filename, struct_as_record=False)
    t_events = None; y_events = None
    fs_ppg = None; t_ppg = None; y_ppg = None; y_ppg_r = None
    fs_ecg = None; t_ecg = None; y_ecg = None

    if 't_ppg' in m:
        fs_ppg = m['fs_ppg'][0][0]
        t_ppg = m['t_ppg'][0]
        y_ppg = m['y_ppg'][0]
        y_ppg_r = m['y_ppg_r'][0]
    if 'events' in m:
        t_events = m['events'][:,1]
        y_events = m['events'][:,0] #[:,1:]
    if 't_ecg' in m:
        fs_ecg = m['fs_ecg'][0][0]
        t_ecg = m['t_ecg'][0]
        y_ecg = m['y_ecg'][0]
    
    return(t_ppg, y_ppg, y_ppg_r, fs_ppg, t_ecg, y_ecg, fs_ecg, t_events, y_events)


def Filtration(y, fs, fc, ord=2, ftype='lowpass'):

    """
    Фільтрація сигналу

    Parameters
    ----------
    y - синал, що має чатоту дискретизації fs (в Гц)
    fc - частота зрізу фільтра в Гц (corner freaquency)
    ord - порядок фільтру, ціле число
    ftype='lowpass' - Фільтр низьких частот (ФНЧ), пропускає низькі частоти, та послаблює частоти, розташовані вище частоти зрізу фільтру (fc). 
                     Наприклад фільтр з fc=10 пропустить сигнали з частотами <10Гц
    ftype='highpass' - Фільтр верхніх частот (ФВЧ), пропускає високочастотні сигнали, але послаблює (зменшує амплітуду) сигналів з частотами нижче частоти зрізу. 
                      Наприклад фільтр з fc=10 пропустить сигнали з частотами >10Гц 
    ftype='bandpass' - Смуговий фільтр, пропускає сигнали в певному діапазоні (смузі) частот, і послаблює (вирізає) сигнали частот за межами цієї смуги. 
                      Наприклад, смуговий фільтр з fc=[10,20] пропускає тільки сигнали, частота яких лежить в інтервалі 10-20 Гц. 

    Results
    ----------
    yf

    Examples
    --------
    """

    if isinstance(fc, list): fc = np.array(fc)
    fnq = fs / 2
    fc_nrm = fc / fnq
    b, a = signal.butter(ord, fc_nrm, btype=ftype)
    yf = signal.filtfilt(b, a, y)
    return yf  

def GetHeartRate(t, yf, fs=500, ppg_systolic_peaks=True, signal_type='ppg'):
  # # завантажуємо сигнал фотоплетизмограми (ФПГ) в масиви t (час в сек), y ("об'єм крові в пальці", що міняється в часі)
  # # частота дискретизаціх сигналу 500 Гц
  # # в експерименті піддослідний почав пити чай на 100 секунді, завершив пити на 500 (?) секунді
  # d=pmat.read_mat(filename)
  # t = d['t_ppg']
  # y = d['y_ppg']
  # # фільтрація сигналу ФПГ (сигнал після фільтрації - yf)
  # yf = filt_butter(y, fs=500, fc=[0.5, 10], ord=2, ftype='bandpass')
  # #шукаємо систолічні ( nk.ppg_findpeaks( yf, ...) )
  # #або діастолічні піки ( nk.ppg_findpeaks( -yf, ...) ) сигналу ФПГ
  if signal_type=='ppg':
    if ppg_systolic_peaks:
        pks = nk.ppg_findpeaks( yf , sampling_rate=fs) #method="bishop" ["elgendi"]
    else:
        pks = nk.ppg_findpeaks( -yf , sampling_rate=fs) #method="bishop" ["elgendi"]
    ip=pks['PPG_Peaks']
  if signal_type=='ecg':
    pks = nk.ecg_findpeaks( yf , sampling_rate=fs) 
    ip=pks['ECG_R_Peaks']  
  tP = t[ip]
  P = yf[ip]
  PP = np.diff(tP)
  tPP = tP[:-1]
  tHR = tPP
  HR = 60 / PP
  return(tP,P, tPP,PP, tHR,HR)


def OtliersCorrectionNNI(tNN,NN, p=95, dbg=False):
    dNN = np.abs(np.diff(NN))
    thr = np.percentile(dNN, p)
    io = np.where( dNN>thr )[0]
    iodel = np.where( np.diff(io)>1 )[0] + 1
    iodel = np.insert(iodel,0, 0)
    io = np.delete(io,iodel)
    if dbg:
        PPlot([tNN,tNN[io]],[dNN,dNN[io]],mode=['lines+markers','markers'])
        PPlot([tNN,tNN[io]],[NN,NN[io]],mode=['lines+markers','markers'])
    NN2=NN.copy()
    NN2[io]=np.nan
    NN2 = hrva.interpolate_nan_values(NN2, interpolation_method='linear')
    if dbg:
        PPlot([tNN,tNN],[NN,NN2],mode=['lines','lines'])
    return( np.array(NN2) )

def SignalCutSegment(t,y,t1,t2=None,time_from_zero=0):
    if (t1 is not None) and (t2 is None):
        i = np.where( t>t1 )
    elif (t1 is None) and (t2 is not None):
        i = np.where( t<=t2 )
    elif (t1 is not None) and (t2 is not None):
        i = np.where( (t>t1) & (t<=t2) )
    else:
        i =  range(len(y))
    segm_t=t[i]
    segm_y=y[i]
    if time_from_zero:
        segm_t -= segm_t[0]
    return(segm_t, segm_y)     

def SimpleFFT(t,y, t1=None,t2=None):
  fs = 1/np.mean(np.diff(t))
  if not t1:
    t1=t[0]
  if not t2:
    t2=t[-1]
  y = y[(t>=t1)&(t<=t2)]
  t = t[(t>=t1)&(t<=t2)]
  A = scipy.fft.rfft(y, norm='forward')
  A = np.abs( A ) * 2
  Fr = scipy.fft.rfftfreq( len(y), 1 / fs)
  return(Fr,A)  

def SpectrogramSTFT(y, sf, frband=None, segmleng=1000, plottype='plotly'):
    #The spectrogram is a two-dimensional representation of the squared magnitude of the STFT:

    freqs, bins, Zxx = signal.stft(y, sf, nperseg=segmleng)
    Zxx=np.abs(Zxx)**2

    if frband:
        ifr = np.where( (freqs>frband[0]) & (freqs<frband[1])  )[0]
        freqs = freqs[ifr]
        Zxx = Zxx[ifr,:]

    if plottype.casefold()=='plotly':
        trace = [go.Heatmap(
            x= bins,
            y= freqs,
            z= Zxx,#10*np.log10(Zxx),
            colorscale='Jet',
            )]
        layout = go.Layout(
            #title = 'Spectrogram with plotly',
            yaxis = dict(title = 'Frequency [Hz]'), # x-axis label
            xaxis = dict(title = 'Time [sec]'), # y-axis label
            )
        fig = go.Figure(data=trace, layout=layout)
        pyo.iplot(fig, filename='Spectrogram')

    if plottype.casefold()=='matplotlib':
        #plt.pcolormesh(bins, freqs, np.abs(Zxx), vmin=0, vmax=amp, shading='gouraud')
        plt.pcolormesh(bins, freqs, np.abs(Zxx), shading='gouraud')
        #plt.title('STFT Magnitude')
        plt.ylabel('Frequency [Hz]')
        plt.xlabel('Time [sec]')
        #plt.ylim(0.01, 0.3)
        plt.show() 

def InterpNNI(x, y, new_sample_rate):
    xnew = np.arange(x[0], x[-1], 1/new_sample_rate)
    ynew = np.interp(xnew, x, y)
    return(xnew, ynew)

####################################################################  

def BPlot(x,y,label=None,h=150,w=1000,marker=None,**kargs):
  multlin=0
  if ( type(y[0]) is np.ndarray) | ( type(y[0]) is list):
      multlin=1
      
  if multlin & (not label):
      label = [str(l) for l in np.arange(1,len(y)+1)]
  colors = itertools.cycle(palette)
  fig=bh.figure(height=h, tools=['box_zoom','hover','reset','undo','redo','pan','save','crosshair'],sizing_mode="scale_width")
  fig.toolbar.active_drag = None
  if multlin:
      for i in np.arange(len(y)):
          col = next(colors)
          fig.line(x[i],y[i], legend_label=label[i], color=col, **kargs)
          if marker:
              fig.circle(x[i],y[i], color=col, legend_label=label[i])
          fig.legend.click_policy = "hide" #"mute"
  else:
      fig.line(x,y, **kargs)
      if marker:
          fig.circle(x,y)
          
  bh.show(fig)

################################################################################

def PPlot(x=None, y=None, 
          err_up=None, err_dw=None, abserr=1, label=None, 
          h=300, w=None, mode=["lines+markers"], 
          data=None, data_err=None,
          xlabel=None, ylabel=None, title=None, fontsize=12, showlegend=True):

    if data:
        x=[]; y=[]
        for k in range(len(data)):
            x.append(data[k][0])
            y.append(data[k][1])
    if data_err:
        err_up=[]; err_dw=[]
        for k in range(len(data_err)):
            err_dw.append(data_err[k][0])
            err_up.append(data_err[k][1])
    
    multlin=0
    if ( type(y[0]) is np.ndarray) | ( type(y[0]) is list):
        multlin=1
    
    if multlin & (not label):
        label = [str(l) for l in np.arange(1,len(y)+1)]
    elif not multlin:
        label = [None]
        x=[x]
        y=[y]
        if err_up: err_up = [err_up]
        if err_dw: err_dw = [err_dw]

    if (err_up is not None) & (err_dw is None):
        err_dw = err_up

    if np.isscalar(mode):
        mode = [mode]
    if multlin & (len(mode)==1):
        mode *= len(y)

    layout = go.Layout(autosize=True, width=w, height=h,  margin=dict(l=10, r=10, t=20, b=10))
    fig = go.Figure(layout=layout)

    for i in np.arange(len(y)):
        error_y = None
        if abserr & (err_up is not None):
            error_y=dict(type='data', 
                         array=np.array(err_up[i])-np.array(y[i]), 
                         arrayminus=np.array(y[i])-np.array(err_dw[i]), 
                         visible=True)
        elif (not abserr) & (err_up is not None):
            error_y=dict(type='data', 
                         array=err_up[i], 
                         arrayminus=err_dw[i], 
                         visible=True)
        fig.add_scatter(x=x[i], y=y[i], name=label[i], error_y=error_y, mode=mode[i])#**kargs)

        fig.update_layout(xaxis_title=xlabel,
                          yaxis_title=ylabel,
                          #legend_title="Legend Title",
                          title=title,
                          font=dict(
                          #    family="Courier New, monospace",
                              size=fontsize,
                          #    color="RebeccaPurple"
                          ),
                          showlegend=showlegend 
                         )
        
    fig.show()
