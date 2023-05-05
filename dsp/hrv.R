#library(R.matlab) # nolint
library(pacman)
#p_load(rmatio, gdata, boot, ggplot2, RHRV, clipr, Metrics, tictoc, nortest, quantreg)
p_load(RHRV)

get_hrvprms <- function(rr, fd_interp_rri=4, td_win_size=60, fd_stft_win_size=60, fd_stft_win_shift=5, fd_stft_mean = T)
{
    ######################
    #Frequency Analysis
    ######################
    FR_INTERP_RRI <- fd_interp_rri #4 #Hz 
    ULF_L <- 0;    ULF_R <- 0.04
    VLF_L <- 0;    VLF_R <- 0.04
    LF_L <- 0.04;  LF_R <- 0.15 
    HF_L <- 0.15;  HF_R <- 0.4
    ######################
    #STFT
    ######################
    STFT_WIN_SEC <- fd_stft_win_size #default 5min
    STFT_WIN_SHIFT_SEC <- fd_stft_win_shift #default 30
    ######################
    #periodogram using a FFT
    ######################
    RRI_DETRENDING <- T
    WIN_DANIELL_SMOOTH <- 5

    
    hd = CreateHRVData()
    hd = SetVerbose(hd,FALSE)
    hd$Beat = data.frame( Time=cumsum( c(0,rr) )*1e-3 )
    hd = BuildNIHR(hd)
    hd$Beat = hd$Beat[-1,]


    #########################################
    ####Time Analysis
    #########################################
    hd = CreateTimeAnalysis(hd,size=td_win_size,interval = 7.8125)
    
    #########################################
    ####Frequency Analysis
    #########################################
    #hd = FilterNIHR(hd, long=50, last=10, minbpm=25, maxbpm=180)
    if(FR_INTERP_RRI>0) hd = InterpolateNIHR (hd, freqhr = FR_INTERP_RRI)


    ###################################################
    ####Frequency Analysis of Nonstationary Signals
    ####Short-Time Fourier Transform (STFT)
    ###################################################
    hd = CreateFreqAnalysis(hd)
    hd = CalculatePowerBand(hd, indexFreqAnalysis= 1, size = STFT_WIN_SEC, shift = STFT_WIN_SHIFT_SEC,
                            ULFmin = ULF_L, ULFmax = ULF_R, VLFmin = VLF_L, VLFmax = VLF_R,
                            LFmin = LF_L, LFmax = LF_R, HFmin = HF_L, HFmax = HF_R)
    TF <- hd$FreqAnalysis[1][[1]]$HRV
    ULF <- hd$FreqAnalysis[1][[1]]$ULF
    VLF <- hd$FreqAnalysis[1][[1]]$VLF
    LF <- hd$FreqAnalysis[1][[1]]$LF
    HF <- hd$FreqAnalysis[1][[1]]$HF
    LFHF <- hd$FreqAnalysis[1][[1]]$LFHF
    HFn <- HF/(LF+HF)*100 #HFn <- HF/(TF-VLF)*100
    LFn <- LF/(LF+HF)*100 #LFn <- LF/(TF-VLF)*100

    ###################################################
    ####Frequency Analysis of Stationary Signals
    ####periodogram using a FFT (estimates the PSD using the FFT and optionally smooths the estimate with Daniell smoothers)
    ###################################################

    hd = CreateFreqAnalysis(hd)
    hd=CalculatePSD(hd,indexFreqAnalysis=2, method="pgram", doPlot = F,
                    detrend=RRI_DETRENDING, #taper=PR_TAP, 
                    spans = WIN_DANIELL_SMOOTH,
                    ULFmin = ULF_L, ULFmax = HF_R, VLFmin = VLF_L, VLFmax = VLF_R,
                    LFmin = LF_L, LFmax = LF_R, HFmin = HF_L, HFmax = HF_R)
    E_PG <- CalculateEnergyInPSDBands(hd, indexFreqAnalysis=2,
                    ULFmin = ULF_L, ULFmax = HF_R, VLFmin = VLF_L, VLFmax = VLF_R,
                    LFmin = LF_L, LFmax = LF_R, HFmin = HF_L, HFmax = HF_R)


    ###################################################
    ####Frequency Analysis of Stationary Signals
    ####Lomb-Scargle periodogram (estimates the frequencies present on the signal by applying a least squares fit of sinusoids to the data samples)
    ###################################################
    hd2 = CreateHRVData()
    hd2 = SetVerbose(hd,FALSE)
    hd2$Beat = data.frame( Time=cumsum( c(0,rr) )*1e-3 )
    hd2 = BuildNIHR(hd2)
    hd2$Beat = hd2$Beat[-1,]

    hd2 = CreateFreqAnalysis(hd2)
    hd2=CalculatePSD(hd2,indexFreqAnalysis=3,method="lomb", doPlot = F)
    E_LB <- CalculateEnergyInPSDBands(hd2, indexFreqAnalysis=3,
                    ULFmin = ULF_L, ULFmax = HF_R, VLFmin = VLF_L, VLFmax = VLF_R,
                    LFmin = LF_L, LFmax = LF_R, HFmin = HF_L, HFmax = HF_R)

    

    #########RES2###############
    if(fd_stft_mean){
      TF_ST <- mean(TF)
      ULF_ST <- mean(ULF)
      VLF_ST <- mean(VLF)
      LF_ST <- mean(LF)
      HF_ST <- mean(HF)
      LFHF_ST <- mean(LFHF)
      HFn_ST <- mean(HFn) 
      LFn_ST <- mean(LFn)   
    }else{
      TF_ST <- TF
      ULF_ST <- ULF
      VLF_ST <- VLF
      LF_ST <- LF
      HF_ST <- HF
      LFHF_ST <- LFHF
      HFn_ST <- HFn 
      LFn_ST <- LFn   
    }
    ########################### 

    #########RES3###############
    ULF_PG <- E_PG[1]
    VLF_PG <- E_PG[2]
    LF_PG <- E_PG[3]
    HF_PG <- E_PG[4]
    LFn_PG <- LF_PG/(LF_PG+HF_PG)*100
    HFn_PG <- HF_PG/(LF_PG+HF_PG)*100
    LFHF_PG <- LF_PG/HF_PG
    ###########################

    #########RES4###############
    ULF_LB <- E_LB[1]
    VLF_LB <- E_LB[2]
    LF_LB <- E_LB[3]
    HF_LB <- E_LB[4]
    LFn_LB <- LF_LB/(LF_LB+HF_LB)*100
    HFn_LB <- HF_LB/(LF_LB+HF_LB)*100
    LFHF_LB <- LF_LB/HF_LB
    ###########################

    r <- list( 
      SDNN=hd$TimeAnalysis[[1]]$SDNN, #1
      SDANN=hd$TimeAnalysis[[1]]$SDANN, #2
      SDNNIDX=hd$TimeAnalysis[[1]]$SDNNIDX, #3
      pNN50=hd$TimeAnalysis[[1]]$pNN50, #4
      SDSD=hd$TimeAnalysis[[1]]$SDSD, #5
      rMSSD=hd$TimeAnalysis[[1]]$rMSSD, #6
      IRRR=hd$TimeAnalysis[[1]]$IRRR, #7
      MADRR=hd$TimeAnalysis[[1]]$MADRR, #8
      TINN=hd$TimeAnalysis[[1]]$TINN, #9
      HRVi=hd$TimeAnalysis[[1]]$HRVi, #10
      mNN=mean(rr), #11
      mHR=mean(60000/rr), #12
      medianNN=median(rr), #13
      #####################################
      TF_ST=TF_ST, #14
      ULF_ST=ULF_ST, #15
      VLF_ST=VLF_ST, #16
      LF_ST=LF_ST, #17
      HF_ST=HF_ST, #18
      LFn_ST=LFn_ST, #19
      HFn_ST=HFn_ST, #20
      LFHF_ST=LFHF_ST, #21
      #####################################
      ULF_PG=ULF_PG, #22
      VLF_PG=VLF_PG, #23
      LF_PG=LF_PG, #24
      HF_PG=HF_PG, #25
      LFn_PG=LFn_PG, #26
      HFn_PG=HFn_PG, #27
      LFHF_PG=LFHF_PG, #28
      #####################################
      ULF_LB=ULF_LB, #29
      VLF_LB=VLF_LB, #30
      LF_LB=LF_LB, #31
      HF_LB=HF_LB, #32
      LFn_LB=LFn_LB, #33
      HFn_LB=HFn_LB, #34
      LFHF_LB=LFHF_LB, #35
      #####################################
      #PSD_ST_FR=hd$FreqAnalysis[[1]]$periodogram$freq, 
      #PSD_ST_PW=hd$FreqAnalysis[[1]]$periodogram$spec,
      PSD_PG_FR=hd$FreqAnalysis[[2]]$periodogram$freq, 
      PSD_PG_PW=hd$FreqAnalysis[[2]]$periodogram$spec,
      PSD_LB_FR=hd2$FreqAnalysis[[3]]$periodogram$freq, 
      PSD_LB_PW=hd2$FreqAnalysis[[3]]$periodogram$spec
    )   

    #names(r)  <- NMS_HRVPRMS0

    return(r)

}


get_td_hrvprms <- function(rr, td_win_size=60)
{
    hd = CreateHRVData()
    hd = SetVerbose(hd,FALSE)
    hd$Beat = data.frame( Time=cumsum( c(0,rr) )*1e-3 )
    hd = BuildNIHR(hd)
    hd$Beat = hd$Beat[-1,]

    #########################################
    ####Time Analysis
    #########################################
    hd = CreateTimeAnalysis(hd,size=td_win_size,interval = 7.8125)
    
    r <- c( hd$TimeAnalysis[[1]]$SDNN, #1
      hd$TimeAnalysis[[1]]$SDANN, #2
      hd$TimeAnalysis[[1]]$SDNNIDX, #3
      hd$TimeAnalysis[[1]]$pNN50, #4
      hd$TimeAnalysis[[1]]$SDSD, #5
      hd$TimeAnalysis[[1]]$rMSSD, #6
      hd$TimeAnalysis[[1]]$IRRR, #7
      hd$TimeAnalysis[[1]]$MADRR, #8
      hd$TimeAnalysis[[1]]$TINN, #9
      hd$TimeAnalysis[[1]]$HRVi, #10
      mean(rr), #11
      mean(60000/rr), #12
      median(rr) #13
    )   

    names(r)  <- NMS_HRVPRMS0[1:13]

    return(r)

}



calc_nonstat_prms <- function(rr, fd_interp_rri=4, fd_stft_win_size=60, fd_stft_win_shift=5)
{
    ######################
    #Frequency Analysis
    ######################
    FR_INTERP_RRI <- fd_interp_rri #4 #Hz 
    ULF_L <- 0;    ULF_R <- 0.04
    VLF_L <- 0;    VLF_R <- 0.04
    LF_L <- 0.04;  LF_R <- 0.15 
    HF_L <- 0.15;  HF_R <- 0.4
    ######################
    #STFT
    ######################
    STFT_WIN_SEC <- fd_stft_win_size #default 5min
    STFT_WIN_SHIFT_SEC <- fd_stft_win_shift #default 30
    ######################

    
    hd = CreateHRVData()
    hd = SetVerbose(hd,FALSE)
    hd$Beat = data.frame( Time=cumsum( c(0,rr) )*1e-3 )
    hd = BuildNIHR(hd)
    hd$Beat = hd$Beat[-1,]

    
    #########################################
    ####Frequency Analysis
    #########################################
    #hd = FilterNIHR(hd, long=50, last=10, minbpm=25, maxbpm=180)
    hd = InterpolateNIHR (hd, freqhr = FR_INTERP_RRI)

    ###################################################
    ####Frequency Analysis of Nonstationary Signals
    ####Short-Time Fourier Transform (STFT)
    ###################################################
    hd = CreateFreqAnalysis(hd)
    hd = CalculatePowerBand(hd, indexFreqAnalysis= 1, size = STFT_WIN_SEC, shift = STFT_WIN_SHIFT_SEC,
                            ULFmin = ULF_L, ULFmax = ULF_R, VLFmin = VLF_L, VLFmax = VLF_R,
                            LFmin = LF_L, LFmax = LF_R, HFmin = HF_L, HFmax = HF_R)
    TF <- hd$FreqAnalysis[1][[1]]$HRV
    #ULF <- hd$FreqAnalysis[1][[1]]$ULF
    VLF <- hd$FreqAnalysis[1][[1]]$VLF
    LF <- hd$FreqAnalysis[1][[1]]$LF
    HF <- hd$FreqAnalysis[1][[1]]$HF
    #LFHF <- hd$FreqAnalysis[1][[1]]$LFHF
    #HFn <- HF/(LF+HF)*100 #HFn <- HF/(TF-VLF)*100
    #LFn <- LF/(LF+HF)*100 #LFn <- LF/(TF-VLF)*100

    TF_S <- sd(TF)
    #ULF_S <- sd(ULF)
    VLF_S <- sd(VLF)
    LF_S <- sd(LF)
    HF_S <- sd(HF)
    #LFHF_S <- sd(LFHF)
    #HFn_S <- sd(HFn) 
    #LFn_S <- sd(LFn)   

    TF_IQR <- IQR(TF)
    VLF_IQR <- IQR(VLF)
    LF_IQR <- IQR(LF)
    HF_IQR <- IQR(HF)
  
    r <- c( 
      TF_S, #1
      VLF_S, #2
      LF_S, #3
      HF_S, #4
      TF_IQR, #5
      VLF_IQR, #6
      LF_IQR, #7
      HF_IQR #8
    )   

    names(r)  <- NMS_STATPRMS0

    return(r)
}


calc_td_nonstat_prms <- function(trr, rr, 
                                 td_win_size=60,   #
                                 segm_size=60,     #
                                 segm_shift=5)     #
{
    t1 <- trr[1]
    t2 <- trr[length(trr)]-segm_size
    if(t1>t2){
      t2 <- t1
      cat('warn calc_td_nonstat_prms, segm_size=',segm_size)
    }

    T1_SEGM <- seq(t1, t2, segm_shift)
    
    hrv_segm <- c()
    for( t1_segm in T1_SEGM ){
      t2_segm <- t1_segm + segm_size
      i_segm <- which( trr>=t1_segm & trr<=t2_segm )
      segm <- rr[ i_segm ]
      hrv_segm_ <- get_td_hrvprms(segm, td_win_size)
      hrv_segm <- rbind(hrv_segm, hrv_segm_)
    }

    r1 <- apply(hrv_segm,2,sd)
    r2 <- apply(hrv_segm,2,IQR)
    names(r1) <- paste0( names(r1), '_s')
    names(r2) <- paste0( names(r2), '_iqr')
    return( c(r1,r2) )
}


                 #   1      2       3      4      5        6         7        8      
NMS_STATPRMS0 <- c("TF_s","VLF_s","LF_s","HF_s","TF_iqr","VLF_iqr","LF_iqr","HF_iqr")
           
               #   1      2        3        4      5       6       7      8      
NMS_HRVPRMS0 <- c("SDNN","SDANN","SDNNIDX","pNN50","SDSD","rMSSD","IRRR","MADRR",
               #   9      10     11    12     13        
                 "TINN","HRVi","mNN","mHR","medianNN",  
               #     14        15     16      17     18        19       20       21
                 "TF_1","ULF_1","VLF_1","LF_1","HF_1", "LFn_1","HFn_1","LF/HF_1",
               #    22       23      24       25       26       27        28         
                 "ULF_2","VLF_2","LF_2","HF_2", "LFn_2","HFn_2","LF/HF_2",
               #    29       30      31       32       33       34        35         
                 "ULF_3","VLF_3","LF_3","HF_3", "LFn_3","HFn_3","LF/HF_3" )

           
               #   1      2        3        4      5       6       7      8      
NMS_HRVPRMS <- c("SDNN","SDANN","SDNNIDX","pNN50","SDSD","rMSSD","IRRR","MADRR",
               #   9      10     11    12     13        
                 "TINN","HRVi","mNN","mHR","medianNN",  
               #     14        15     16      17     18        19       20       21
                 "TF ST","ULF ST","VLF ST","LF ST","HF ST", "LFn ST","HFn ST","LF/HF ST",
               #    22       23      24       25       26       27        28         
                 "ULF PG","VLF PG","LF PG","HF PG", "LFn PG","HFn PG","LF/HF PG",
               #    29       30      31       32       33       34        35         
                 "ULF LB","VLF LB","LF LB","HF LB", "LFn LB","HFn LB","LF/HF LB" )
                #      1            2             3           4             5           6    
NMS_HRVPRMS2 <- c("SDNN (ms)","SDANN (ms)","SDNNIDX (ms)","pNN50 (%)","SDSD (ms)","rMSSD (ms)",
               #       7          8            9        10       11         12             13        
                 "IRRR (ms)","MADRR (ms)","TINN (ms)","HRVi","mNN (ms)","mHR (bmp)","medianNN (ms)",  
               #       14          15              16            17            18              19           20             21
                 "TF ST (ms2)","ULF ST (ms2)","VLF ST (ms2)","LF ST (ms2)","HF ST (ms2)", "LFn ST (%)","HFn ST (%)","LF/HF ST",
               #      22             23            24            25              26           27          28         
                 "ULF PG (ms2)","VLF PG (ms2)","LF PG (ms2)","HF PG (ms2)", "LFn PG (%)","HFn PG (%)","LF/HF PG",
               #      29             30            31            32              33           34          35         
                 "ULF LB (ms2)","VLF LB (ms2)","LF LB (ms2)","HF LB (ms2)", "LFn LB (%)","HFn LB (%)","LF/HF LB")

#######################################################################################################
#######################################################################################################
#######################################################################################################