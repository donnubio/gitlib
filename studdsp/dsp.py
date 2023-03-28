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

import heartpy as hp
from heartpy.datautils import rolling_mean

import neurokit2 as nk

import sys

####################################################################  

def ЗавантажитиРезультатиРеєстраціїФПГ(імя_файла):
  m = sio.loadmat(імя_файла, struct_as_record=False)
  fs_ppg = m['fs_ppg'][0][0]
  t_ppg = m['t_ppg'][0]
  y_ppg = m['y_ppg'][0]
  y_ppg_r = m['y_ppg_r'][0]
  if 'events' in m:
  	t_events = m['events'][:,0]
  	y_events = m['events'][:,1] #[:,1:]
  else:
  	t_events = None
  	y_events = None

  return(t_ppg, y_ppg, y_ppg_r, fs_ppg, t_events, y_events)


def ФільтрБаттерворта(y, fs, fc, ord=2, ftype='lowpass'):
  #y - синал, що має чатоту дискретизації fs (в Гц)
  #fc - частота зрізу фільтра в Гц (corner freaquency)
  #ord - порядок фільтру, ціле число
  #ftype='lowpass' - Фільтр низьких частот (ФНЧ), пропускає низькі частоти, та послаблює частоти, розташовані вище частоти зрізу фільтру (fc). 
  #                  Наприклад фільтр з fc=10 пропустить сигнали з частотами <10Гц
  #ftype='highpass' - Фільтр верхніх частот (ФВЧ), пропускає високочастотні сигнали, але послаблює (зменшує амплітуду) сигналів з частотами нижче частоти зрізу. 
  #                  Наприклад фільтр з fc=10 пропустить сигнали з частотами >10Гц 
  #ftype='bandpass' - Смуговий фільтр, пропускає сигнали в певному діапазоні (смузі) частот, і послаблює (вирізає) сигнали частот за межами цієї смуги. 
  #                  Наприклад, смуговий фільтр з fc=[10,20] пропускає тільки сигнали, частота яких лежить в інтервалі 10-20 Гц. 
  if isinstance(fc, list): fc = np.array(fc)
  fnq = fs / 2
  fc_nrm = fc / fnq
  b, a = signal.butter(ord, fc_nrm, btype=ftype)
  yf = signal.filtfilt(b, a, y)
  return yf  

def getHeartRate(yf, частота_дискретизації=500, peaks_max=1):
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
  if peaks_max:
    pks = nk.ppg_findpeaks( yf , sampling_rate=частота_дискретизації) #method="bishop" ["elgendi"]
  else:
    pks = nk.ppg_findpeaks( -yf , sampling_rate=частота_дискретизації) #method="bishop" ["elgendi"]
  ip=pks['PPG_Peaks']
  tP = t[ip]
  P = yf[ip]
  PP = np.diff(tP)
  tPP = tP[:-1]
  tHR = tPP
  HR = 60 / PP
  return(t,yf, tP,P, tPP,PP, tHR,HR)


def signal_segment(t,y,t1,t2,time_from_zero=0):
  i = np.where( (t>t1) & (t<=t2) )
  segm_t=t[i]
  segm_y=y[i]
  if time_from_zero:
      segm_t -= segm_t[0]
  return(segm_t, segm_y)     


def ПеретворенняФурє(t,y, t1=None,t2=None):
  fs = 1/np.mean(np.diff(t))
  if not t1:
    t1=t[0]
  if not t2:
    t2=t[-1]
  y = y[(t>=t1)&(t<=t2)]
  t = t[(t>=t1)&(t<=t2)]
  амплітуди = scipy.fft.rfft(y, norm='forward')
  амплітуди = np.abs( амплітуди ) * 2
  частоти = scipy.fft.rfftfreq( len(y), 1 / fs)
  return(частоти,амплітуди)  

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

def PPlot(data=None, x=None,y=None, data_err=None, err_up=None, err_dw=None, abserr=1, label=None, h=300, w=None, **kargs):

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
        fig.add_scatter(x=x[i],y=y[i],name=label[i],error_y=error_y, **kargs)
        
    fig.show()
