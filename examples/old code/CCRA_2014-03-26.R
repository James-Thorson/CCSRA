

##############################
#
# NOTES
#
# 1. Catch curves
#   - Requires that h=1 (so R is constant) and estimating a single F
#
# 2. Conventional SRA
#   - Doesn't appear to have information to update RecDevs (SE always equals SigmaR)
#
##############################

#######################
# Header
#######################

# File structure
PlatformFile = "C:/Users/James.Thorson/Desktop/UW Hideaway/Collaborations/"
  #PlatformFile = "C:/Users/thorsonja/Dropbox/"
SourceFile = paste(PlatformFile,"2013 -- catch curve reduction analysis/",sep="")
AdmbFile = paste(SourceFile,"ADMB/",sep="")

# Libraries
library(vioplot)

# Functions
source(paste(SourceFile,"Fn_CCRA_2013-11-21.R",sep=""))

# Date file
Date = Sys.Date()
  #Date = "2013-12-16_EstH_EstMweak_SigmaR=0.9_EstRecDev_Fast"
  DateFile = paste(SourceFile,Date,"/",sep="")
  dir.create(DateFile)
FigFile = paste(DateFile,"Figs/",sep="")
  dir.create(FigFile)
  
#######################
# Simulate data
#######################

# Settings
AdmbVersion = "CCRA_v4"   # v3: Added priors on h and M; v4: Added RecDevs
AgeMax = 20
Nyears = 20
Ncomp_per_year = 100
Nmethods = 3  # 1: Catch curve; 2: CC-SRA; 3:DB-SRA
MethodSet = c(-1, 1, 1) # -2: Catch-curve with instantanous F; -1: Catch-curve with continuous F; 1: CCRA with continuous F; 2=CCRA with instantaneous F using Pope's approximation
IncludeDepletionPrior = c(FALSE, FALSE, TRUE) # 0: Off; 1: On
IncludeAgeComp = c(TRUE, TRUE, FALSE) # 0: Off; 1: On
                                # Slow=Periodic (high-steepness, late-maturity, high survival) "red snapper" from fishbase
# Biological parameters         # Fast=Opportunistic (low-steepness, early maturity, low survival) "Pacific sardine" from fishbase
SpeciesType = 2 # 1=Slow; 2=Fast
K = c(0.1, 0.2)[SpeciesType]        # Slow: K=0.1; Fast: K=0.2
Linf = c(60, 30)[SpeciesType]       # Slow: Linf=60; Fast: Linf=30
LMASS = c(2, 1)[SpeciesType] # log-maximum annual spawners per spawner ; Slow: LMASS=2 (h=0.83); Fast: LMASS=1 (h=0.65)
L0 = 1
Amat = log( 3*(Linf-L0)/Linf ) / K
S50 = Amat
Sslope = 1
W_alpha = 0.01
W_beta = 3.04
R0 = 1e9
SigmaR = 0.9

# Effort dynamics  parameters
Fequil = 0.25
Frate = 0.2
F1 = 0.1
SigmaF = 0.2

# Derived
M = 1.84 * K
L_a = Linf - (Linf - L0) * exp(-K*0:AgeMax)
W_a = W_alpha * L_a^W_beta      # In grams
S_a = 1 / (1 + exp( -Sslope * (0:AgeMax - S50) ))
M_a = pnorm( (0:AgeMax-Amat)/(0.25*Amat) * 1.96 )
LMLSS = LMASS - log(M)
h = exp(LMLSS) / (4 + exp(LMLSS))
SB0 = sum( R0 * exp(-M * 0:AgeMax) * W_a * M_a )
SBPR0 = SB0 / R0

##### Simulate data
# Recruits (e.g., N_at[1,1]) arise from spawning biomass in that year (i.e., N_at[,1])
# Changes in effort (e.g., F[t+1]) arise from spawning biomass in the earlier year (e.g., SB[t])
# Catches (e.g., C[t]) arise from removals in that year (i.e., F[t])
# Cw_t -- Catch (weight) in year t
# Cn_at -- Catch (numbers) in year t and age a
# Dn_at -- Natural mortality (numbers) in year t and age a
# Zn_at -- total mortality (numbers) in year t and age a

# Simulation settings
Nrep = 200
RandomSeed = ceiling(runif(1, min=1,max=1e6))

# Saving objects
TimeseriesResults = array(NA, dim=c(1+Nmethods,Nrep,3,Nyears,2), dimnames=list(c("True",paste("Method=",1:Nmethods)),paste("Rep=",1:Nrep),c("ln_F_t","ln_SB_t","ln_D_t"),paste("Year=",1:Nyears),c("Est","SE")))
RecDevResults = array(NA, dim=c(1+Nmethods,Nrep,Nyears+AgeMax,2), dimnames=list(c("True",paste("Method=",1:Nmethods)),paste("Rep=",1:Nrep),c(paste("Age=",AgeMax:0),paste("Year=",2:Nyears)),c("Est","SE")))
ParamResults = array(NA, dim=c(1+Nmethods,Nrep,5,2), dimnames=list(c("True",paste("Method=",1:Nmethods)),paste("Rep=",1:Nrep),c("ln_R0","M","h","S50","Sslope"),c("Est","SE")))

# Simulation loop
  # RepI = MethodI = 1
for(RepI in 1:Nrep){
  for(MethodI in 1:Nmethods){
  
    # Simulation settings
    set.seed(RepI + RandomSeed)
    F_method = MethodSet[MethodI]   
  
    # Data objects
    Cw_t = SB_t = F_t = rep(NA, Nyears)
    Zn_at = Dn_at = Cn_at = N_at = matrix(NA, nrow=AgeMax+1, ncol=Nyears)
    RecDev = rnorm( Nyears + AgeMax, mean=-SigmaR^2/2, sd=SigmaR )
    
    # Initialization
    F_t[1] = F1
    N_at[,1] = R0 * exp(-M * 0:AgeMax) * exp( RecDev[1:(AgeMax+1)] )
    SB_t[1] = sum( N_at[,1] * W_a * M_a )
    if(F_method==-1 | F_method==1){
      Zn_at[,1] = N_at[,1] * (1 - exp( -M - F_t[1]*S_a ))
      Dn_at[,1] = Zn_at[,1] * (M) / (M + F_t[1]*S_a)
      Cn_at[,1] = Zn_at[,1] * (F_t[1]*S_a) / (M + F_t[1]*S_a)
    }
    if(F_method==-2 | F_method==2){
      Dn_at[,1] = N_at[,1] * (1 - exp(-M/2))
      Cn_at[,1] = N_at[,1] * exp(-M/2) * F_t[1]*S_a
      Dn_at[,1] = Dn_at[,1] + N_at[,1] * exp(-M/2) * (1 - F_t[1]*S_a) * (1 - exp(-M/2))
      Zn_at[,1] = Dn_at[,1] + Cn_at[,1]
    }
    
    # Projection
    for(YearI in 2:Nyears){
      # Fishing effort
      F_t[YearI] = F_t[YearI-1] * (SB_t[YearI-1] / (Fequil * SB0))^Frate * exp( rnorm(1, mean=-SigmaF^2/2, sd=SigmaF) )
        if( F_t[YearI] > 0.95 & (F_method==-2 | F_method==2) ) F_t[YearI] = 0.95
      # Survival
      N_at[-1,YearI] = N_at[-AgeMax,YearI-1] - Zn_at[-AgeMax,YearI-1]
      # Spawning biomass
      SB_t[YearI] = sum( (N_at[,YearI] * W_a * M_a)[-1] )
      # Recruitment
      N_at[1,YearI] = 4 * h * R0 * SB_t[YearI] / ( SB0*(1-h) + SB_t[YearI]*(5*h-1) ) * exp( RecDev[AgeMax+YearI] )
      # Removals
      if(F_method==-1 | F_method==1){
        Zn_at[,YearI] = N_at[,YearI] * (1 - exp( -M - F_t[YearI]*S_a ))
        Dn_at[,YearI] = Zn_at[,YearI] * (M) / (M + F_t[YearI]*S_a)
        Cn_at[,YearI] = Zn_at[,YearI] * (F_t[YearI]*S_a) / (M + F_t[YearI]*S_a)
      }
      if(F_method==-2 | F_method==2){
        Dn_at[,YearI] = N_at[,YearI] * (1 - exp(-M/2))
        Cn_at[,YearI] = N_at[,YearI] * exp(-M/2) * F_t[YearI]*S_a
        Dn_at[,YearI] = Dn_at[,YearI] + N_at[,YearI] * exp(-M/2) * (1 - F_t[YearI]*S_a) * (1 - exp(-M/2))
        Zn_at[,YearI] = Dn_at[,YearI] + Cn_at[,YearI]
      }
    }
    Cw_t = (W_a %*% Cn_at)[1,]
    
    # Plot projection
    matplot( cbind(SB_t/SB0,F_t,Cw_t/max(Cw_t)), type="l", col=c("black","red","blue"), lty="solid")
    
    #######################
    # Estimate model
    #######################
    
    # Generate data
    AgeComp_at = array(0, dim=dim(N_at))
    if(IncludeAgeComp[MethodI]==TRUE){
      #AgeComp_at[,Nyears-1] = rmultinom(n=1, size=Ncomp_per_year, prob=Cn_at[,Nyears-1])[,1]
      AgeComp_at[,Nyears] = rmultinom(n=1, size=Ncomp_per_year, prob=Cn_at[,Nyears])[,1]
    }
    if(F_method==-1 | F_method==-2) Cw_t_input = rep(0, length(Cw_t))
    if(F_method==1 | F_method==2) Cw_t_input = Cw_t
    
    # Transformations
    h_alpha = ((h-0.2)/0.8) * ( ((h-0.2)/0.8) * (1 - ((h-0.2)/0.8)) / 0.1^2 -1)
    h_beta = (1 - ((h-0.2)/0.8)) * (((h-0.2)/0.8) * (1-((h-0.2)/0.8)) / 0.1^2 - 1)
    if( F_method==-2 ){
      F_t_HI = 0.95
      F_t_length = 1
    }
    if( F_method==-1 ){
      F_t_HI = 3
      F_t_length = 1
    }
    if( F_method==1 ){
      F_t_HI = 3
      F_t_length = Nyears
    }
    if( F_method==2 ){
      F_t_HI = 0.95
      F_t_length = 0
    }
    F_t_phase = 2
    
    # Priors
    ln_R0_prior = c(10, 30, 20, 999, 999, 1)
    F_t_prior = c( 0, F_t_HI, 0.1, 999, 999, F_t_phase, F_t_length )
    h_prior = c(0.2, 1.0, ifelse(F_method<0,0.9999,h), h_alpha, h_beta, ifelse(F_method<0,-1,4))
    M_prior = c(0, 1, M, M, 0.5, 4)
    S50_prior = c(999, 999, ifelse(IncludeAgeComp[MethodI]==TRUE,5,S50), 999, 999, ifelse(IncludeAgeComp[MethodI]==TRUE,3,-1))
    Sslope_prior = c(999, 999, ifelse(IncludeAgeComp[MethodI]==TRUE,1,10), 999, 999, ifelse(IncludeAgeComp[MethodI]==TRUE,3,-1))
    D_prior = c(0.4, 0.2, ifelse(IncludeDepletionPrior[MethodI],1,0))
    SigmaR_prior = c( 0, 1, SigmaR, 0.6, 0.1, -1)
    RecDev_prior = c( -3, 3, 0, 999, 999, 5 )
    
    # Fit twice for bias adjustment if estimating recruitment deviations
    if(RecDev_prior[6]>0){ Nloop = 2 }else{ Nloop = 1 }
    for(LoopI in 1:Nloop){
      # Bias adjustment for each loop
      if(LoopI==1) RecDev_biasadj = rep(0, Nyears+AgeMax)
      if(LoopI==2 & paste(tolower(AdmbVersion),".std",sep="")%in%list.files(DateFile) & file.info(paste(DateFile,AdmbVersion,".std",sep=""))$size>0){
        STD = read.table( file=paste(DateFile,AdmbVersion,".std",sep=""), skip=1)        
        RecDev_biasadj = 1 - STD[which(STD[,'V2']=="RecDev"),'V4']^2 / SigmaR_prior[3]^2
      } 
    
      # Write to file
      write(c("# Nyears",Nyears), file=paste(DateFile,"CCRA.dat",sep=""), append=FALSE)
        write(c("# AgeMax",AgeMax), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# F_method", F_method), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# ln_R0_prior", ln_R0_prior), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# M_prior", M_prior), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# h_prior", h_prior), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# S50_prior", S50_prior), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# Sslope_prior", Sslope_prior), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# F_t_prior", F_t_prior), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# D_prior", D_prior), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# SigmaR_prior", SigmaR_prior), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# RecDev_prior", RecDev_prior), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# RecDev_biasadj", RecDev_biasadj), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# Cw_t", Cw_t_input), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# W_a", W_a), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write(c("# M_a", M_a), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write("# AgeComp_at", file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
        write.table(AgeComp_at, file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE, row.names=FALSE, col.names=FALSE)
        write(c("# Test value", 123456), file=paste(DateFile,"CCRA.dat",sep=""), append=TRUE)
      
      # Run model
      CleanAdmbFn(SaveFile=DateFile, PrefixVec=tolower(AdmbVersion), KeepVec=c("dat"))
        setwd(DateFile)
        file.copy(from=paste(AdmbFile,AdmbVersion,".exe",sep=""), to=paste(DateFile,AdmbVersion,".exe",sep=""), overwrite=TRUE)
        shell(paste(AdmbVersion,".exe",sep=""), intern=TRUE)
        Sys.sleep(0.5)
    
      # Read in timeseries output
      if( paste(tolower(AdmbVersion),".std",sep="")%in%list.files(DateFile) & file.info(paste(DateFile,AdmbVersion,".std",sep=""))$size>0 ){
        STD = read.table( file=paste(DateFile,AdmbVersion,".std",sep=""), skip=1)
        ln_F_t_hat = STD[which(STD[,'V2']=="ln_F_t"),c('V3','V4')]
        ln_SB_t_hat = STD[which(STD[,'V2']=="ln_SB_t"),c('V3','V4')]
        ln_D_t_hat = STD[which(STD[,'V2']=="ln_D_t"),c('V3','V4')]
        RecDev_hat = STD[which(STD[,'V2']=="RecDev_hat"),c('V3','V4')]
        Param_hat = STD[which(STD[,'V2']=="Param_hat"),c('V3','V4')]
      }else{
        ln_F_t_hat = cbind( scan(file=paste(DateFile,AdmbVersion,".rep",sep=""), skip=1, nlines=1, quiet=TRUE), NA )
        ln_SB_t_hat = cbind( scan(file=paste(DateFile,AdmbVersion,".rep",sep=""), skip=3, nlines=1, quiet=TRUE), NA )
        ln_D_t_hat = cbind( scan(file=paste(DateFile,AdmbVersion,".rep",sep=""), skip=5, nlines=1, quiet=TRUE), NA )
        RecDev_hat = cbind( scan(file=paste(DateFile,AdmbVersion,".rep",sep=""), skip=7, nlines=1, quiet=TRUE), NA )
        Param_hat = cbind( scan(file=paste(DateFile,AdmbVersion,".rep",sep=""), skip=9, nlines=1, quiet=TRUE), NA )
      }
      # Correct output for catch curves
      if(F_method==-1 | F_method==-2){
        ln_SB_t_hat[] = ln_D_t_hat[] = NA
        Param_hat[c(1,3),] = NA
      }
      # Correct output for conventional SRA
      if(IncludeAgeComp[MethodI]==FALSE){
        RecDev_hat[] = NA
      } 
    }  # End fitting loop

    # Possible correct non-converged models
    if(FALSE){
      if( all(exp(ln_F_t_hat[1:Nyears,1])<0.01) & (!is.na(ln_D_t_hat[Nyears,1]) & exp(ln_D_t_hat[Nyears,1])>1) ){
        stop("check convergence")
      }
    }
    
    # Record status results
    TimeseriesResults[paste("Method=",MethodI),RepI,'ln_F_t',,] = as.matrix(ln_F_t_hat)
    TimeseriesResults[paste("Method=",MethodI),RepI,'ln_SB_t',,] = as.matrix(ln_SB_t_hat)
    TimeseriesResults[paste("Method=",MethodI),RepI,'ln_D_t',,] = as.matrix(ln_D_t_hat)
    # Record recruitment deviations
    RecDevResults[paste("Method=",MethodI),RepI,,] = as.matrix(RecDev_hat)
    # Record parameter estimates
    ParamResults[paste("Method=",MethodI),RepI,c("ln_R0","M","h","S50","Sslope"),] = as.matrix(Param_hat[1:5,])
  
    # Debugging
    if(FALSE){
      par(mfrow=c(1,2))
      matplot( cbind(exp(ln_F_t_hat[,1]),F_t), type="l", col=c("red","black"), lty="solid")
      matplot( cbind(exp(ln_SB_t_hat[,1]),SB_t), type="l", col=c("red","black"), lty="solid")
    }
  }  # End Method loop

  # Record status results
  TimeseriesResults['True',RepI,'ln_F_t',,'Est'] = log(F_t)
  TimeseriesResults['True',RepI,'ln_SB_t',,'Est'] = log(SB_t)
  TimeseriesResults['True',RepI,'ln_D_t',,'Est'] = log(SB_t / SB0)
  #print( cbind(t(TimeseriesResults[,RepI,'ln_F_t',,'Est']), t(TimeseriesResults[,RepI,'ln_F_t',,'SE'])) )
  # Record recruitment deviations
  RecDevResults['True',RepI,,'Est'] = RecDev
  # Record parameter estimates
  ParamResults['True',RepI,c("ln_R0","M","h","S50","Sslope"),'Est'] = c( log(R0), M, h, S50, Sslope )
  
  # Plot results
  png( file=paste(FigFile,"Timeseries_RepI=",RepI,".png",sep=""), width=6, height=6, res=200, units="in")
    par(mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
    for(VarI in 1:4){
      if(VarI==4){ Mat = RecDevResults[,RepI,,] }else{ Mat = TimeseriesResults[,RepI,c('ln_F_t','ln_SB_t','ln_D_t')[VarI],,]}
      matplot( exp(t(Mat[,,'Est'])), type="l", lty="solid", col=c("black","red","green","blue"), ylim=range(na.rm=TRUE,exp(rbind(Mat[,,'Est'],Mat[,,'Est']+1.96*Mat[,,'SE']))), lwd=2, xlab="Year", ylab="", main=c("F","SB","D","RecDev")[VarI], log="y" )
      for(MethodI in 1:length(MethodSet)){
        polygon( x=c(1:dim(Mat)[2],dim(Mat)[2]:1), y=exp(c(Mat[paste("Method=",MethodI),,'Est']+1.96*Mat[paste("Method=",MethodI),,'SE'],rev(Mat[paste("Method=",MethodI),,'Est']-1.96*Mat[paste("Method=",MethodI),,'SE']))), col=list(rgb(1,0,0,alpha=0.2),rgb(0,1,0,alpha=0.2),rgb(0,0,1,alpha=0.2))[[MethodI]], border=NA )
      }
    }
  dev.off() 
}

## Save objects
if( !("Results.RData" %in% list.files(DateFile)) ){
  # Results
  Results = list('TimeseriesResults'=TimeseriesResults, 'RecDevResults'=RecDevResults, 'ParamResults'=ParamResults)
  save(Results , file=paste(DateFile,"Results.RData",sep=""))
  # Cleaning
  CleanAdmbFn(SaveFile=DateFile, PrefixVec=tolower(AdmbVersion), KeepVec=c("dat"))
  # Record
  RecordList = list(AdmbVersion=AdmbVersion, AgeMax=AgeMax, Nyears=Nyears, Ncomp_per_year=Ncomp_per_year, Nmethods=Nmethods, MethodSet=MethodSet, IncludeDepletionPrior=IncludeDepletionPrior, IncludeAgeComp=IncludeAgeComp, K=K, Linf=Linf, L0=L0, Amat=Amat, S50=S50, Sslope=Sslope, W_alpha=W_alpha, W_beta=W_beta, R0=R0, LMASS=LMASS, SigmaR=SigmaR, Fequil=Fequil, Frate=Frate, F1=F1, SigmaF=SigmaF, M=M, L_a=L_a, W_a=W_a, S_a=S_a, M_a=M_a, LMLSS=LMLSS, h=h, SB0=SB0, SBPR0=SBPR0, Nrep=Nrep, RandomSeed=RandomSeed)
  capture.output(RecordList, file=paste(DateFile,"RecordList.txt",sep=""))
  save(RecordList, file=paste(DateFile,"RecordList.RData",sep=""))
}

## Summary figure
TimeseriesEst = TimeseriesResults[-1,,,Nyears,'Est']
TimeseriesTrue = outer(rep(1,length(MethodSet)),TimeseriesResults['True',,,Nyears,'Est']) 
ParamEst = ParamResults[-1,,,'Est']
ParamTrue = outer(rep(1,length(MethodSet)),ParamResults['True',,,'Est'])

# Plot time series results
square = function(Obj) Obj^2
png( file=paste(FigFile,"Timeseries_Summary.png",sep=""), width=16, height=8, res=200, units="in")
  par(mfrow=c(2,4), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(0,2,0,0))
  for(VarI in 1:3){
    boxplot( t(TimeseriesEst[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]]-TimeseriesTrue[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]]), main=c("F","SB","D")[VarI], ylim=range(na.rm=TRUE,TimeseriesEst[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]]-TimeseriesTrue[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]]) + c(-0,0.2)*diff(range(na.rm=TRUE,TimeseriesEst[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]]-TimeseriesTrue[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]])) )
    RMSE = sqrt(colMeans(square(t(TimeseriesEst[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]]-TimeseriesTrue[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]])), na.rm=TRUE))
    Bias = colMeans(t(TimeseriesEst[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]]-TimeseriesTrue[,,c('ln_F_t','ln_SB_t','ln_D_t')[VarI]]), na.rm=TRUE)
    text( x=1:Nmethods, y=rep(par()$usr[4],2), labels=paste(format(RMSE,digits=2,nsmall=3),'\n',format(Bias,digits=2,nsmall=3),sep=""), pos=1, cex=2)
    abline(h=0, lwd=3, lty="dotted")
  }
  for(ParamI in 1:5){
    boxplot( t(ParamEst[,,ParamI]-ParamTrue[,,ParamI]), main=dimnames(ParamResults)[[3]][ParamI], ylim=range(na.rm=TRUE,ParamEst[,,ParamI]-ParamTrue[,,ParamI]) + c(-0,0.2)*diff(range(na.rm=TRUE,ParamEst[,,ParamI]-ParamTrue[,,ParamI])) )
    RMSE = sqrt(colMeans(square(t(ParamEst[,,ParamI]-ParamTrue[,,ParamI])), na.rm=TRUE))
    Bias = colMeans(t(ParamEst[,,ParamI]-ParamTrue[,,ParamI]), na.rm=TRUE)
    text( x=1:Nmethods, y=rep(par()$usr[4],2), labels=paste(format(RMSE,digits=2,nsmall=3),'\n',format(Bias,digits=2,nsmall=3),sep=""), pos=1, cex=2)
    abline(h=0, lwd=3, lty="dotted")
  }
  mtext(side=2, text="Error (log-space)", outer=TRUE) 
dev.off()
    