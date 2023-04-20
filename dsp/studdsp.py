import dsp
import dsp.hrrec as hrrec

####################################################################  

def ЗавантажитиДані( ІмяФайлу ):
   (t_ppg, y_ppg, y_ppg_r, fs_ppg, t_ecg, y_ecg, fs_ecg, t_events, y_events) = dsp.LoadData(ІмяФайлу)
   return(t_ppg, y_ppg, y_ppg_r, fs_ppg, t_ecg, y_ecg, fs_ecg, t_events, y_events)

def ФільтраціяСигналу(сигнал, частота_дискретизац, частота_зрізу, порядок_фільтра=2, тип_фільтра='lowpass'):
   відфільтрований_сигнал = dsp.Filtration(сигнал, частота_дискретизац, частота_зрізу, ord=порядок_фільтра, ftype=тип_фільтра)
   return відфільтрований_сигнал
  #y - синал, що має чатоту дискретизації fs (в Гц)
  #fc - частота зрізу фільтра в Гц (corner freaquency)
  #ord - порядок фільтру, ціле число
  #ftype='lowpass' - Фільтр низьких частот (ФНЧ), пропускає низькі частоти, та послаблює частоти, розташовані вище частоти зрізу фільтру (fc). 
  #                  Наприклад фільтр з fc=10 пропустить сигнали з частотами <10Гц
  #ftype='highpass' - Фільтр верхніх частот (ФВЧ), пропускає високочастотні сигнали, але послаблює (зменшує амплітуду) сигналів з частотами нижче частоти зрізу. 
  #                  Наприклад фільтр з fc=10 пропустить сигнали з частотами >10Гц 
  #ftype='bandpass' - Смуговий фільтр, пропускає сигнали в певному діапазоні (смузі) частот, і послаблює (вирізає) сигнали частот за межами цієї смуги. 
  #                  Наприклад, смуговий фільтр з fc=[10,20] пропускає тільки сигнали, частота яких лежить в інтервалі 10-20 Гц. 

def ОбчислЧастотуСерцевихСкорочень(час, сигнал, частота_дискретизац, шукати_систолічні_піки_ФПГ=True, тип_сигналу='ppg'):
   (tP,P, tPP,PP, tHR,HR) = dsp.GetHeartRate(час, сигнал, частота_дискретизац, шукати_систолічні_піки_ФПГ, тип_сигналу)
   return(tP,P, tPP,PP, tHR,HR)

def ВирізатиФрагментСигналу(час, сигнал, час_початок_фрагменту, час_кінець_фрагменту, час_фрагменту_з_нуля=False):
   (segm_t, segm_y) = dsp.SignalCutSegment(час, сигнал, час_початок_фрагменту, час_кінець_фрагменту, час_фрагменту_з_нуля)
   return(segm_t, segm_y)

def КорекціяВикідівNNI(час, сигналNNI, імовірність=95, будувати_графіки=False):
   сигналNNI = dsp.OtliersCorrectionNNI(час,сигналNNI, p=імовірність, dbg=будувати_графіки)
   return(сигналNNI)

def ПеретворенняФурє(час, сигнал, час_початок_фрагменту=None, час_кінець_фрагменту=None):
   частоти, амплітуди = dsp.SimpleFFT(час, сигнал, час_початок_фрагменту, час_кінець_фрагменту)
   return(частоти, амплітуди)  

def Спектрограма(сигнал, частота_дискретизац, смуга_частот=None, довжина_сегмента=1000, plottype='plotly'):
   dsp.SpectrogramSTFT(y=сигнал, sf=частота_дискретизац, frband=смуга_частот, segmleng=довжина_сегмента, plottype='plotly')

def ІнтерполяціяNNI(час_NNI, сигнал_NNI, нова_частота_дискретизації):
   час_NNI_і, сигнал_NNI_і = dsp.InterpNNI(час_NNI, сигнал_NNI, нова_частота_дискретизації)
   return(час_NNI_і, сигнал_NNI_і)
   
def ЛогістичнаФункція(x, ymin, ymax, xmdl, tau):
   return hrrec.LogisticFun(x, ymin, ymax, xmdl, tau)

def ПрипасуватиЛогістичнуФункціюДоВідновленняЧСС(час, ЧСС, початок_фрагмента=None, кінець_фрагмента=None):
   (params_fit, HR_fit, covr, tHR, HR) = hrrec.FitLogisticFun(tHR=час, HR=ЧСС, t1=початок_фрагмента, t2=кінець_фрагмента)
   return(params_fit, HR_fit, covr, tHR, HR)

####################################################################  
from dsp import BPlot
from dsp import PPlot
