__version__ = 'dev'

import numpy as np
#import scipy.io as sio
from scipy import signal


def filt_butter(y, fs, fc, ord=2, ftype='lowpass'):
    #функція filt_butter здійснює цифрову фільтрацію сигнала
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
    y = signal.filtfilt(b, a, y)
    return y


def signal_segment( t,y, t1, t2 ):
    i = np.where( (t>=t1) & (t<=t2) )
    segm_t=t[i]
    segm_y=y[i]
    return(segm_t, segm_y)
