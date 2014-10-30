

##############################
#
# NOTES
#
# 1. Catch curves
#   - Requires that h=1 (so expected R is constant) and estimating a single F
#
# 2. Conventional SRA
#   - Doesn't appear to have information to update RecDevs (SE always equals SigmaR)
#
# 3. Simulate data
#   - Recruits (e.g., N_at[1,1]) arise from spawning biomass in that year (i.e., N_at[,1])
#   - Changes in effort (e.g., F[t+1]) arise from spawning biomass in the earlier year (e.g., SB[t])
#   - Catches (e.g., C[t]) arise from removals in that year (i.e., F[t])
#   - Cw_t -- Catch (weight) in year t
#   - Cn_at -- Catch (numbers) in year t and age a
#   - Dn_at -- Natural mortality (numbers) in year t and age a
#   - Zn_at -- total mortality (numbers) in year t and age a
#
##############################

#######################
# Header
#######################

# Install package
#install.packages("devtools")
library("devtools")
install_github("James-Thorson/CCSRA")

# Libraries
library(CCSRA)
library(TMB)

# File structure
TmbFile = paste0(system.file("executables", package="CCSRA"),"/")

# Date file
Date = Sys.Date()
  DateFile = paste0(getwd(),"/",Date,"/")
  dir.create(DateFile)
FigFile = paste0(DateFile,"Figs/")
  dir.create(FigFile)

# Compile model
Version = "CCSRA_v4"   # v3: Added priors on h and M; v4: Added RecDevs
setwd(TmbFile)
compile( dynlib(Version) )
  
#######################
# Settings
#######################

# General
AgeMax = 20
Nyears = 20
Ncomp_per_year = 100
MethodSet = c("CC","CCSRA","SRA") 
Nmethods = 3  # 1: Catch curve; 2: CC-SRA; 3:DB-SRA
                                # Slow=Periodic (high-steepness, late-maturity, high survival) "red snapper" from fishbase
# Biological parameters         # Fast=Opportunistic (low-steepness, early maturity, low survival) "Pacific sardine" from fishbase
SpeciesType = 1 # 1=Slow; 2=Fast
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
SigmaR = 0.4

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
Mat_a = pnorm( (0:AgeMax-Amat)/(0.25*Amat) * 1.96 )
LMLSS = LMASS - log(M)
h = exp(LMLSS) / (4 + exp(LMLSS))
SB0 = sum( R0 * exp(-M * 0:AgeMax) * W_a * Mat_a )
SBPR0 = SB0 / R0

# Simulation settings
Nrep = 10
RandomSeed = ceiling(runif(1, min=1,max=1e6))

#######################
# Loop through simulation experiment
#######################

# Saving objects
TimeseriesResults = array(NA, dim=c(1+Nmethods,Nrep,3,Nyears,2), dimnames=list(c("True",paste("Method=",1:Nmethods)),paste("Rep=",1:Nrep),c("ln_F_t","ln_SB_t","ln_D_t"),paste("Year=",1:Nyears),c("Est","SE")))
RecDevResults = array(NA, dim=c(1+Nmethods,Nrep,Nyears+AgeMax,2), dimnames=list(c("True",paste("Method=",1:Nmethods)),paste("Rep=",1:Nrep),c(paste("Age=",AgeMax:0),paste("Year=",2:Nyears)),c("Est","SE")))
ParamResults = array(NA, dim=c(1+Nmethods,Nrep,5,2), dimnames=list(c("True",paste("Method=",1:Nmethods)),paste("Rep=",1:Nrep),c("ln_R0","M","h","S50","Sslope"),c("Est","SE")))

# Simulation loop
  RepI=1;  MethodI=2
for(RepI in 1:Nrep){
  for(MethodI in 1:Nmethods){
  
    # Simulation settings
    set.seed(RepI + RandomSeed)
    F_method = switch(MethodSet[MethodI], "CC"=-1, "CCSRA"=1, "SRA"=1)
    
    # Simulate data
    DataList = SimData_Fn( Nyears=Nyears, AgeMax=AgeMax, SigmaR=SigmaR, M=M, F1=F1, W_a=W_a, S_a=S_a, Mat_a=Mat_a, h=h, SB0=SB0, Frate=Frate, Fequil=Fequil, SigmaF=SigmaF, Ncomp_per_year=Ncomp_per_year )
    DataList[['AgeComp_at']][,1:(Nyears-1)] = 0
    
    # Plot projection
    matplot( cbind(DataList[['SB_t']]/SB0,DataList[['F_t']],DataList[['Cw_t']]/max(DataList[['Cw_t']])), type="l", col=c("black","red","blue"), lty="solid")
    
    #######################
    # Estimate model
    #######################
    
    #if(F_method==-1 | F_method==-2) Cw_t_input = rep(0, length(Cw_t))
    #if(F_method==1 | F_method==2) Cw_t_input = Cw_t
    
    # Transformations
    Nloop = 2
    
    # Fit twice for bias adjustment if estimating recruitment deviations
    LoopI = 1
    for(LoopI in 1:Nloop){
      # Bias adjustment for each loop
      if(LoopI==1) RecDev_biasadj = rep(0, Nyears+AgeMax)
      if(LoopI==2 && !("condition" %in% names(attributes(Sdreport))) ){
        SD = summary(Sdreport)      
        RecDev_biasadj = 1 - SD[which(rownames(SD)=="RecDev_hat"),'Std. Error']^2 / Report$SigmaR^2
      } 
     
      # Format inputs
      Method = c("CC", "CCSRA", "SRA")[MethodI]
      InputList = FormatInput_Fn(Method=Method, M_prior=c(M,0.5), h_prior=c(h,0.1), D_prior=c(0.4,0.2), SigmaR_prior=c(0.6,0.2), AgeComp_at=DataList[['AgeComp_at']], Cw_t=DataList[['Cw_t']], W_a=W_a, Mat_a=Mat_a, RecDev_biasadj=RecDev_biasadj)
      
      # Compile 
      dyn.load( paste0(TmbFile,dynlib(Version)) )
      Obj <- MakeADFun(data=InputList[['Data']], parameters=InputList[['Parameters']], map=InputList[['Map']], random=InputList[['Random']], inner.control=list(maxit=100) )
      Obj$env$inner.control <- c(Obj$env$inner.control, "step.tol"=1e-12, "tol10"=1e-8, "grad.tol"=1e-12) 
      Obj$fn( Obj$par )
    
      # Run
      Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-14) )
      Opt[["final_gradient"]] = Obj$gr( Opt$par )
      Report = Obj$report()
      Sdreport = try( sdreport(Obj) )
    
      # Read in timeseries output
      if( !("condition" %in% names(attributes(Sdreport))) ){
        SD = summary(Sdreport)
        ln_F_t_hat = SD[which(rownames(SD)=="ln_F_t"),c('Estimate','Std. Error')]
        ln_SB_t_hat = SD[which(rownames(SD)=="ln_SB_t"),c('Estimate','Std. Error')]
        ln_D_t_hat = SD[which(rownames(SD)=="ln_D_t"),c('Estimate','Std. Error')]
        RecDev_hat = SD[which(rownames(SD)=="RecDev_hat"),c('Estimate','Std. Error')]
        Param_hat = SD[which(rownames(SD)=="Param_hat"),c('Estimate','Std. Error')]
      }else{
        ln_F_t_hat = cbind( Report$ln_F_t, NA )
        ln_SB_t_hat = cbind( Report$ln_SB_t, NA )
        ln_D_t_hat = cbind( Report$ln_D_t, NA )
        RecDev_hat = cbind( Report$RecDev_hat, NA )
        Param_hat = cbind( Report$Param_hat, NA )
      }
      # Correct output for catch curves
      if(Method=="CC"){
        ln_SB_t_hat[] = ln_D_t_hat[] = NA
        Param_hat[c(1,3),] = NA
      }
      # Correct output for conventional SRA
      if(Method=="SRA"){
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
  
  }  # End Method loop

  # Record status results
  TimeseriesResults['True',RepI,'ln_F_t',,'Est'] = log(DataList$F_t)
  TimeseriesResults['True',RepI,'ln_SB_t',,'Est'] = log(DataList$SB_t)
  TimeseriesResults['True',RepI,'ln_D_t',,'Est'] = log(DataList$SB_t / SB0)
  # Record recruitment deviations
  RecDevResults['True',RepI,,'Est'] = DataList$RecDev
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
  RecordList = list("R.version"=R.version, AdmbVersion=AdmbVersion, AgeMax=AgeMax, Nyears=Nyears, Ncomp_per_year=Ncomp_per_year, Nmethods=Nmethods, MethodSet=MethodSet, K=K, Linf=Linf, L0=L0, Amat=Amat, S50=S50, Sslope=Sslope, W_alpha=W_alpha, W_beta=W_beta, R0=R0, LMASS=LMASS, SigmaR=SigmaR, Fequil=Fequil, Frate=Frate, F1=F1, SigmaF=SigmaF, M=M, L_a=L_a, W_a=W_a, S_a=S_a, Mat_a=Mat_a, LMLSS=LMLSS, h=h, SB0=SB0, SBPR0=SBPR0, Nrep=Nrep, RandomSeed=RandomSeed)
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
    