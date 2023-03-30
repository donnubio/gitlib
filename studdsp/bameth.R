library(boot) 
library(quantreg)

BA_sfun <- function(data, idx, ba_ci_conf=0.95) 
{
# v1 
    D <- data$D[idx]
    A <- data$A[idx]

    al <- (1-ba_ci_conf)/2
  
    ###### 1 ######
    #D=const, параметрично
    Bias_1 <- mean(D)
    LoA_L_1 <- Bias_1 - 1.96*sd(D)
    LoA_U_1 <- Bias_1 + 1.96*sd(D)
    Min_1 <- min(D)
    Max_1 <- max(D)
    SD_1 <- sd(D)

    ###### 2 ######
    #D=const, непараметрично
    r <- GetQuantiles(D, c(al,0.5,1-al))
    Bias_2 <- r[2] #q0.5
    LoA_L_2 <- r[1] #q0.025
    LoA_U_2 <- r[3] #q0.975
    DLoA_2 <- LoA_U_2 - LoA_L_2
    #Min_2 <- 
    #Max_2 <- 

    ###### 3 ######
    #D=a*A+b, непараметрично
    r <- summary( rq(D ~ A, c(al, 0.5, 1-al)) )
    Bias_b <- r[[2]]$coefficients[1,1] #Intercept,  q0.5
    Bias_a <- r[[2]]$coefficients[2,1] #Slope,      q0.5
    #Bias_b_L <- r[[2]]$coefficients[1,2]
    #Bias_b_U <- r[[2]]$coefficients[1,3]
    #Bias_a_L <- r[[2]]$coefficients[2,2]
    #Bias_a_U <- r[[2]]$coefficients[2,3]
    LoA_L_b <- r[[1]]$coefficients[1,1] #Intercept, q0.025 
    LoA_L_a <- r[[1]]$coefficients[2,1] #Slope,     q0.025
    LoA_U_b <- r[[3]]$coefficients[1,1] #Intercept, q0.975
    LoA_U_a <- r[[3]]$coefficients[2,1] #Slope,     q0.975

    ###### 4 ######
    #|D|=const, непараметрично
    absD <- abs(D)
    r <- GetQuantiles(absD, c(1-al))
    TDI <- r #q0.975
    Max_4 <- max(absD)
    r <- summary( rq(absD ~ A, c(0.5, 0.95, 1-al)) )
    TDI_Q_05_b <- r[[1]]$coefficients[1,1] #Intercept,  q0.5
    TDI_Q_05_a <- r[[1]]$coefficients[2,1] #Slope,      q0.5
    TDI_Q_095_b <- r[[2]]$coefficients[1,1] #Intercept,  q0.95
    TDI_Q_095_a <- r[[2]]$coefficients[2,1] #Slope,      q0.95
    TDI_Q_0975_b <- r[[3]]$coefficients[1,1] #Intercept,  q0.975
    TDI_Q_0975_a <- r[[3]]$coefficients[2,1] #Slope,      q0.975

    ###### 5 ######
    #min max values for quantiles regression lines in range [Amin Amax]
    Amin <- min(data$A)
    Amax <- max(data$A)
    Bias_min <- Bias_a*Amin + Bias_b
    Bias_max <- Bias_a*Amax + Bias_b
    LoA_L_min <- LoA_L_a*Amin + LoA_L_b
    LoA_L_max <- LoA_L_a*Amax + LoA_L_b
    LoA_U_min <- LoA_U_a*Amin + LoA_U_b
    LoA_U_max <- LoA_U_a*Amax + LoA_U_b
    DLoA_min <- LoA_U_min - LoA_L_min
    DLoA_max <- LoA_U_max - LoA_L_max
    TDI_Q_05_min <- TDI_Q_05_a*Amin + TDI_Q_05_b
    TDI_Q_05_max <- TDI_Q_05_a*Amax + TDI_Q_05_b
    TDI_Q_095_min <- TDI_Q_095_a*Amin + TDI_Q_095_b
    TDI_Q_095_max <- TDI_Q_095_a*Amax + TDI_Q_095_b
    TDI_Q_0975_min <- TDI_Q_0975_a*Amin + TDI_Q_0975_b
    TDI_Q_0975_max <- TDI_Q_0975_a*Amax + TDI_Q_0975_b

    ###### 6 ######
    D0 <- D - (Bias_a*A + Bias_b)
    res <-cor.test(A, abs(D0),  method = "spearman")
    BA_ro <- res$estimate
    #res <-cor.test(A, abs(D),  method = "spearman")
    #TDI_ro <- res$estimate

    o <- c(
        #1
        Bias_1=unname(Bias_1),
        LoA_L_1=unname(LoA_L_1),
        LoA_U_1=unname(LoA_U_1),
        Min_1=unname(Min_1),
        Max_1=unname(Max_1),
        SD_1=unname(SD_1),
        #2
        Bias_2=unname(Bias_2),
        LoA_L_2=unname(LoA_L_2),
        LoA_U_2=unname(LoA_U_2),
        DLoA_2=unname(DLoA_2),
        #3
        Bias_b=unname(Bias_b),
        Bias_a=unname(Bias_a),
        LoA_L_b=unname(LoA_L_b),
        LoA_L_a=unname(LoA_L_a),
        LoA_U_b=unname(LoA_U_b),
        LoA_U_a=unname(LoA_U_a),
        #4
        TDI=unname(TDI),
        Max_4=unname(Max_4),
        TDI_Q_05_b=unname(TDI_Q_05_b),
        TDI_Q_05_a=unname(TDI_Q_05_a),
        TDI_Q_095_b=unname(TDI_Q_095_b),
        TDI_Q_095_a=unname(TDI_Q_095_a),
        TDI_Q_0975_b=unname(TDI_Q_0975_b),
        TDI_Q_0975_a=unname(TDI_Q_0975_a),
        #5
        Bias_min=unname(Bias_min),
        Bias_max=unname(Bias_max),
        LoA_L_min=unname(LoA_L_min),
        LoA_L_max=unname(LoA_L_max),
        LoA_U_min=unname(LoA_U_min),
        LoA_U_max=unname(LoA_U_max),
        DLoA_min=unname(DLoA_min),
        DLoA_max=unname(DLoA_max),
        TDI_Q_05_min=unname(TDI_Q_05_min),
        TDI_Q_05_max=unname(TDI_Q_05_max),
        TDI_Q_095_min=unname(TDI_Q_095_min),
        TDI_Q_095_max=unname(TDI_Q_095_max),
        TDI_Q_0975_min=unname(TDI_Q_0975_min),
        TDI_Q_0975_max=unname(TDI_Q_0975_max),
        #6
        BA_ro=unname(BA_ro)
        )



    return(o)
} 

NUM_BA_STATISTICS <- 39
NMS_BA_STATISTICS <- c(
        'Bias_1','LoA_L_1','LoA_U_1','Min_1','Max_1','SD_1',
        'Bias_2','LoA_L_2','LoA_U_2','DLoA_2',
        'Bias_b','Bias_a','LoA_L_b','LoA_L_a','LoA_U_b','LoA_U_a',
        'TDI','Max_4','TDI_Q_05_b','TDI_Q_05_a','TDI_Q_095_b','TDI_Q_095_a','TDI_Q_0975_b','TDI_Q_0975_a',
        'Bias_min','Bias_max','LoA_L_min','LoA_L_max','LoA_U_min','LoA_U_max','DLoA_min','DLoA_max','TDI_Q_05_min','TDI_Q_05_max','TDI_Q_095_min','TDI_Q_095_max','TDI_Q_0975_min','TDI_Q_0975_max',
        'BA_ro')

sv1quantile <- function(x, probs) {
  #Sfakianakis–Verginis
  n <- length(x)
  if (n <= 2)
    return(quantile(x, probs))
  x <- sort(x)
  sapply(probs, function(p) {
    B <- function(x) dbinom(x, n, p)
    B(0) * (x[1] + x[2] - x[3]) / 2 +
      sum(sapply(1:n, function(i) (B(i) + B(i - 1)) * x[i] / 2)) +
      B(n) * (-x[n-2] + x[n-1] + x[n]) / 2
  })
}

GetQuantiles <- function(x,p=c(0.025,0.5,0.975,1)){
    #For sample sizes exceeding 80 observations, more advanced quantile estimators, 
    #such as the Harrell–Davis and estimators of Sfakianakis–Verginis type
    #https://www.mdpi.com/2571-905X/3/3/22
    #v1
    #q <- quantile(x, p, type=3) #na.rm=T
    #v2 (Sfakianakis–Verginis)
    q <- sv1quantile(x, p)
    #v3 (Harrell–Davis)
    #q <- hdquantile(x, p) 
    return(q) 
}

#########################################################################################################
#########################################################################################################
#########################################################################################################

get_BA_data <- function(y1, y2, fl_trsfrm=F, fl_rerr=F, true_val_type="mean"){
    #true_val_type = "mean", "y1", "y2"

    #логарифмічна трансформація даних
    if (fl_trsfrm)
    {
        y1 <- log(y1)
        y2 <- log(y2)
        #потребує зворотньої трансформації (можлво не для всіх типів CI)
    }
    D <- y1-y2
    if (true_val_type=="mean")
        A <- (y1+y2)/2
    if (true_val_type=="y1")
        A <- y1
    if (true_val_type=="y2")
        A <- y2
    #по вісі ординат на діаграмі БлентаАльтмана буде відносна помилка у %
    if (fl_rerr)
    {
        i_na <- which(D==0 & A==0) #наприклад для pNN50 обидва значення можуть = 0 
        D <- D / A * 100 
        D[i_na] <- 0 #помилка = 0
    }    

    return( data.frame(A,D) )

}

#########################################################################################################
#########################################################################################################
#########################################################################################################

ba_boot_from_expr <- function(y1,y2, bs_nrep, fl_trsfrm=F, fl_rerr=F, ba_ci_conf=0.95, ci_conf=0.95, true_val_type="mean"){

# bs_nrep - кільк реплікац
# fl_rerr - БА діаграма, різниці (F) або відносна помилка (T)
# ba_ci_conf=0.95 - LoA
# ci_conf=0.95 - CI для LoA, Bias ....

    ba_stats_ <- BA_sfun(get_BA_data(y1,y2, fl_trsfrm, fl_rerr, true_val_type), 1:length(y1)) 
    NMS_BA_PR <-  names( ba_stats_ )

            DATA <- get_BA_data(y1, y2, fl_trsfrm, fl_rerr, true_val_type)
            bs <- boot(data=DATA, BA_sfun, R = bs_nrep,  parallel = "multicore", ncpus=4,  ,ba_ci_conf=ba_ci_conf) 

            ci1 <- getCI_boot(bs, ci_conf, 'basic')
            ci2 <- getCI_boot(bs, ci_conf, 'norm')
            ci3 <- getCI_boot(bs, ci_conf, 'perc')
            ci4 <- getCI_boot(bs, ci_conf, 'bca')
            #ci5 <- getCI_boot(bs, BA_ci_conf, 'stud')

            BA_PR_0 <- bs$t0 
            BA_PR_MN <- apply(bs$t, 2, mean) #2 by col
            BA_PR_BS_L <- ci1$CI_L; BA_PR_BS_R <- ci1$CI_R
            BA_PR_NR_L <- ci2$CI_L; BA_PR_NR_R <- ci2$CI_R
            BA_PR_PC_L <- ci3$CI_L; BA_PR_PC_R <- ci3$CI_R
            BA_PR_BCA_L <- ci4$CI_L; BA_PR_BCA_R <- ci4$CI_R
            BA_D <- DATA$D 
            BA_A <- DATA$A



    ba <- list(BA_D=BA_D, BA_A=BA_A, 
               BA_PR_0=BA_PR_0, BA_PR_MN=BA_PR_MN, 
               BA_PR_BS_L=BA_PR_BS_L, BA_PR_BS_R=BA_PR_BS_R, 
               BA_PR_NR_L=BA_PR_NR_L, BA_PR_NR_R=BA_PR_NR_R, 
               BA_PR_PC_L=BA_PR_PC_L, BA_PR_PC_R=BA_PR_PC_R, 
               BA_PR_BCA_L=BA_PR_BCA_L, BA_PR_BCA_R=BA_PR_BCA_R,
               BA_NMS_BA_PR=NMS_BA_PR)

    return(ba)
}

#########################################################################################################
#########################################################################################################
#########################################################################################################

getCI_boot <- function(bs,conf=0.95,citype){

#citype = "norm" | "basic" | "perc" | "stud" | "bca"
    
    n_stat <- length(bs$t0)
    CI_L<-c()
    CI_R<-c()
    for(k in seq(1,n_stat)){ #for all statistics
        #cat(k,' ')
        #CI_L[k] <- NA
        #CI_R[k] <- NA
        res = tryCatch(
            expr = {
                ci <- boot.ci(bs, conf, type = citype, index=k)
                if(citype == "norm"){
                    CI_L[k] <- ci$norm[2]
                    CI_R[k] <- ci$norm[3]}
                else if(citype == "basic"){
                    CI_L[k] <- ci$basic[4]
                    CI_R[k] <- ci$basic[5]}
                else if(citype == "perc"){
                    CI_L[k] <- ci$perc[4]
                    CI_R[k] <- ci$perc[5]}
                else if(citype == "stud"){
                    CI_L[k] <- ci$student[4]
                    CI_R[k] <- ci$student[5]}
                else if(citype == "bca"){
                    CI_L[k] <- ci$bca[4]
                    CI_R[k] <- ci$bca[5]}            
                else{ print('Error citype'); return() }
            },
            error = function(e) {
                CI_L[k] <<- NA
                CI_R[k] <<- NA
            }
        )#res
    }
    
    return( list(CI_L=CI_L,CI_R=CI_R) )
}

###########################################################################
