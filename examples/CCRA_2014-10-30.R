

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
#   - Recruits (numbers at age 0, i.e., N_at[1,1]) arise from spawning biomass in that year (i.e., N_at[,1])
#   - Changes in effort (e.g., F[t+1]) arise from spawning biomass in the earlier year (e.g., SB[t])
#   - Catches (e.g., C[t]) arise from removals in that year (i.e., F[t])
#   - Cw_t -- Catch (weight) in year t
#   - Cn_at -- Catch (numbers) in year t and age a
#   - Dn_at -- Natural mortality (numbers) in year t and age a
#   - Zn_at -- total mortality (numbers) in year t and age a
#
# CONVERGENCE SUGGESTIONS
#  1. Put prior on Sslope so it doesn't go >4, because S50 becomes knife-edge and inestimable
#  2. Explore a prior on logF to keep it <3
#
##############################

setwd("C:/Users/James.Thorson/Desktop/")

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

# If running locally from J. Thorson's laptop
if(TRUE){
  TmbFile = paste( "C:/Users/James.Thorson/Desktop/Project_git/CCSRA/inst/executables/" )
  SourceFile = "C:/Users/James.Thorson/Desktop/Project_git/CCSRA/R/"
  File = list.files( SourceFile )
  for(i in 1:length(File)) source( paste0(SourceFile,File[i]))
}
  
# Date file
Date = paste0(Sys.Date(),"b")
  DateFile = paste0(getwd(),"/",Date,"/")
  dir.create(DateFile)
FigFile = paste0(DateFile,"Figs/")
  dir.create(FigFile)

# Compile model
Version = "CCSRA_v4"   # v3: Added priors on h and M; v4: Added RecDevs
setwd(TmbFile)
dyn.unload( paste0(Version,".dll") )
file.remove( paste(Version,c(".dll",".o"),sep="") )
compile( paste0(Version,".cpp") )

#######################
# Settings
#######################

# General
AgeMax = 20
Nyears = 20
Ncomp_per_year = 100
SurveyCV = 0.4
MethodSet = c("CC", "CCSRA", "SRA", "AS" ) # 1: Catch curve; 2: CC-SRA; 3:DB-SRA; 4: Age-structured 
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
TimeseriesResults = array(NA, dim=c(1+length(MethodSet),Nrep,3,Nyears,2), dimnames=list(c("True",MethodSet),paste("Rep=",1:Nrep),c("ln_F_t","ln_SB_t","ln_D_t"),paste("Year=",1:Nyears),c("Est","SE")))
RecDevResults = array(NA, dim=c(1+length(MethodSet),Nrep,Nyears+AgeMax,2), dimnames=list(c("True",MethodSet),paste("Rep=",1:Nrep),c(paste("Age=",AgeMax:0),paste("Year=",2:Nyears)),c("Est","SE")))
ParamResults = array(NA, dim=c(1+length(MethodSet),Nrep,5,2), dimnames=list(c("True",MethodSet),paste("Rep=",1:Nrep),c("ln_R0","M","h","S50","Sslope"),c("Est","SE")))
BiasCorr = array(NA, dim=c(length(MethodSet),Nrep,Nyears+AgeMax,3), dimnames=list(paste("Method=",MethodSet),paste("Rep=",1:Nrep),c(paste("Age=",AgeMax:0),paste("Year=",2:Nyears)),c("Orig","Epsilon","AdHoc")))

# Simulation loop
  RepI=1;  MethodI=2
for(RepI in 1:Nrep){
  for(MethodI in 1:length(MethodSet)){
  
    # Simulation settings
    set.seed( (RepI+RandomSeed) %% 1e6 )
    Method = MethodSet[MethodI]
    F_method = switch(Method, "CC"=-1, "CCSRA"=1, "SRA"=1, "AS"=1) # 1: Explicit F; 2: Hybrid (not implemented)
    
    # Simulate data
    DataList = SimData_Fn( Nyears=Nyears, AgeMax=AgeMax, SigmaR=SigmaR, M=M, F1=F1, W_a=W_a, S_a=S_a, Mat_a=Mat_a, h=h, SB0=SB0, Frate=Frate, Fequil=Fequil, SigmaF=SigmaF, Ncomp_per_year=Ncomp_per_year, SurveyCV=SurveyCV )
    # Exclude all age-comp except for final year, except for age-structured model
    if( Method%in%c("CC","CCSRA") ) DataList[['AgeComp_at']][,1:(Nyears-1)] = 0
    if( Method=="SRA" ) DataList[['AgeComp_at']][] = 0
    if( Method=="AS" ) DataList[['AgeComp_at']][,1:ceiling(Nyears/2)] = 0
    # Turn off index except for age-structured model
    if( Method%in%c("CC","CCSRA","SRA") ) DataList[['Index_t']][,1] = NA

    # Plot projection
    png( file=paste(FigFile,"Simulation_RepI=",RepI,".png",sep=""), width=3*1, height=3*1, res=200, units="in")
      par(mfrow=c(1,1), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
      matplot( cbind(DataList[['SB_t']]/SB0,DataList[['F_t']],DataList[['Cw_t']]/max(DataList[['Cw_t']])), type="l", col=c("black","red","blue"), lty="solid")
    dev.off()
    
    #######################
    # Estimate model
    #######################
    
    #if(F_method==-1 | F_method==-2) Cw_t_input = rep(0, length(Cw_t))
    #if(F_method==1 | F_method==2) Cw_t_input = Cw_t
    
    # xTransformations
    if(Method%in%c("CC","CCSRA","SRA","AS")) Nloop = 2
    
    # Fit twice for bias adjustment if estimating recruitment deviations
    LoopI = 1
    for(LoopI in 1:Nloop){
      # Bias adjustment for each loop
      if(LoopI==2 && !("condition" %in% names(attributes(Sdreport))) ){
        SD = summary(Sdreport)      
        RecDev_biasadj = 1 - SD[which(rownames(SD)=="RecDev_hat"),'Std. Error']^2 / Report$SigmaR^2
      }else{
        RecDev_biasadj = rep(1, Nyears+AgeMax)
      } 
     
      # Format inputs
      InputList = FormatInput_Fn(Method=Method, M_prior=c(M,0.5), h_prior=c(h,0.1), Sslope_prior=c(1,1,1), D_prior=c(0.4,0.2), SigmaR_prior=c(0.6,0.2), AgeComp_at=DataList[['AgeComp_at']], Index_t=DataList[['Index_t']], Cw_t=DataList[['Cw_t']], W_a=W_a, Mat_a=Mat_a, RecDev_biasadj=RecDev_biasadj)
      # Keep SigmaR for CC and CCSRA
      if(Method%in%c("CC","CCSRA","AS") & "ln_SigmaR"%in%names(InputList$Map)) InputList$Map = InputList$Map[-which("ln_SigmaR"==names(InputList$Map))]
      # Turn off RecDev for SRA (because it is often singular)
      #if(Method%in%c("SRA") & !("RecDev_hat"%in%names(InputList$Map))) InputList$Map = c(InputList$Map, list("RecDev_hat"=factor(rep(NA,Nyears+AgeMax))))
      # Increase starting value for ln_R0
      InputList$Parameters$ln_R0 = 22
      # Assume M and h are known
      if(!("input_h"%in%names(InputList$Map))) InputList$Map = c(InputList$Map, list("input_h"=factor(NA)))
      if(!("ln_M"%in%names(InputList$Map))) InputList$Map = c(InputList$Map, list("ln_M"=factor(NA)))
      # Turn off SigmaR in 2nd loop
      if(LoopI==2 & !("ln_SigmaR"%in%names(InputList$Map)) ) InputList$Map = c( InputList$Map, list("ln_SigmaR"=factor(NA)) )
      # Change priors on F_t
      InputList$Data$F_t_prior = c(0.0001, 2, 0.1, 999, 0.1, 1, Nyears)
      # Change priors S50
      InputList$Data$S50_prior = c(999, 999, 999, S50, 3, 999)
      # Change priors S50
      InputList$Data$Sslope_prior = c(Sslope, Sslope, 999, 999, 1, 999)
      
      # Compile 
      dyn.load( paste0(TmbFile,dynlib(Version)) )
      if(LoopI==1) Obj <- MakeADFun(data=InputList[['Data']], parameters=InputList[['Parameters']], map=InputList[['Map']], random=InputList[['Random']], inner.control=list(maxit=1e3) )
      if(LoopI==2) Obj <- MakeADFun(data=InputList[['Data']], parameters=ParList, map=InputList[['Map']], random=InputList[['Random']], inner.control=list(maxit=1e3) )
      Obj$env$inner.control <- c(Obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
      Obj$fn( Obj$par )
      
      # Set bounds
      Upr = rep(Inf, length(Obj$par))
        Upr[match("ln_SigmaR",names(Obj$par))] = log(2)
        Upr[match("ln_F_t_input",names(Obj$par))] = log(2)
        Upr[match("Sslope",names(Obj$par))] = log(5)
        Upr[match("S50",names(Obj$par))] = AgeMax*1.5
      Lwr = rep(-Inf, length(Obj$par))
      
      # Run
      Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr, control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) )
      Opt[["final_gradient"]] = Obj$gr( Opt$par )
      # Gr = grad( Obj$fn, Opt$par )
      
      # Re-fit 
      for(i in 1:10){
        if( abs(min(Opt[["final_gradient"]]))>0.1 | Opt$message=="false convergence (8)" ){
          # changes
          InputList$Map = c(InputList$Map, list("ln_SigmaR"=factor(NA)))
          InputList$Random = c( InputList$Random, c("S50","Sslope") )
          Obj <- MakeADFun(data=InputList[['Data']], parameters=InputList[['Parameters']], map=InputList[['Map']], random=InputList[['Random']], inner.control=list(maxit=1e3) )
          Opt = nlminb( start=Obj$env$last.par.best[-Obj$env$random]+rnorm(length(Obj$par)), objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr, control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) )
        }else{
          break()
        }
        Opt[["final_gradient"]] = Obj$gr( Opt$par )
      }
      
      # Standard errors
      Report = Obj$report()
      if( "bias.correct" %in% names(formals(sdreport)) ){
        Obj$env$MCcontrol <- list("doMC"=FALSE, "seed"=RandomSeed%%1e6, "n"=1e4)
        Sdreport = try( sdreport(Obj, bias.correct=c(TRUE,FALSE)[LoopI], importance.sample=FALSE) )
        if(LoopI==1){
          Opt0 = Opt
          ParList = Obj$env$parList( Obj$env$last.par.best )
          BiasCorr[MethodI,RepI,,c("Orig","Epsilon")] = cbind( Sdreport$value, Sdreport$unbiased$value )[grep("RecMult_t",names(Sdreport$value)),]
        }
        if(LoopI==2){
          BiasCorr[MethodI,RepI,,"AdHoc"] = Sdreport$value[grep("RecMult_t",names(Sdreport$value))]
        }
      }else{
        Sdreport = try( sdreport(Obj) )
      }
    
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

    # Record status results
    TimeseriesResults[MethodSet[MethodI],RepI,'ln_F_t',,] = as.matrix(ln_F_t_hat)
    TimeseriesResults[MethodSet[MethodI],RepI,'ln_SB_t',,] = as.matrix(ln_SB_t_hat)
    TimeseriesResults[MethodSet[MethodI],RepI,'ln_D_t',,] = as.matrix(ln_D_t_hat)
    # Record recruitment deviations
    RecDevResults[MethodSet[MethodI],RepI,,] = as.matrix(RecDev_hat)
    # Record parameter estimates
    ParamResults[MethodSet[MethodI],RepI,c("ln_R0","M","h","S50","Sslope"),] = as.matrix(Param_hat[1:5,])
  
  }  # End Method loop

  # Record status results
  TimeseriesResults['True',RepI,'ln_F_t',,'Est'] = log(DataList$F_t)
  TimeseriesResults['True',RepI,'ln_SB_t',,'Est'] = log(DataList$SB_t)
  TimeseriesResults['True',RepI,'ln_D_t',,'Est'] = log(DataList$SB_t / SB0)
  # Record recruitment deviations
  RecDevResults['True',RepI,,'Est'] = DataList$RecDev
  # Record parameter estimates
  ParamResults['True',RepI,c("ln_R0","M","h","S50","Sslope"),'Est'] = c( log(R0), M, h, S50, Sslope )
  
  # Plotting
  if( exists("BiasCorr") ){
    png( file=paste(FigFile,"BiasCorr_RepI=",RepI,".png",sep=""), width=3*2, height=3*2, res=200, units="in")
      par(mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
      for(MethodI in 1:length(MethodSet)){
        matplot( BiasCorr[MethodI,RepI,,], type=c("l","p","p"), main=MethodSet[MethodI])
        text( x=mean(par()$usr[1:2]), y=c(1,0.9,0.8)*par()$usr[4], pos=1, labels=formatC(colMeans(BiasCorr[MethodI,RepI,,]),digits=3,format="f"), col=c("black","red","green"))
      }
    dev.off()
  }
  # Plot results
  png( file=paste(FigFile,"Timeseries_RepI=",RepI,".png",sep=""), width=6, height=6, res=200, units="in")
    par(mfrow=c(2,2), mar=c(3,3,2,0), mgp=c(2,0.5,0), tck=-0.02)
    for(VarI in 1:4){
      if(VarI==4){ Mat = RecDevResults[,RepI,,] }else{ Mat = TimeseriesResults[,RepI,c('ln_F_t','ln_SB_t','ln_D_t')[VarI],,]}
      Ylim = range(na.rm=TRUE,exp(rbind(Mat[,,'Est'],Mat[,,'Est']+1.96*Mat[,,'SE'])))
      matplot( exp(t(Mat[,,'Est'])), type="l", lty="solid", col=c("black","red","green","blue","orange"), ylim=c(Ylim[1],Ylim[2]), lwd=2, xlab="Year", ylab="", main=c("F","SB","D","RecDev")[VarI], log="y" )
      for(MethodI in 1:length(MethodSet)){
        polygon( x=c(1:dim(Mat)[2],dim(Mat)[2]:1), y=exp(c(Mat[MethodSet[MethodI],,'Est']+1.96*Mat[MethodSet[MethodI],,'SE'],rev(Mat[MethodSet[MethodI],,'Est']-1.96*Mat[MethodSet[MethodI],,'SE']))), col=list(rgb(1,0,0,alpha=0.2),rgb(0,1,0,alpha=0.2),rgb(0,0,1,alpha=0.2),rgb(1,1,0,alpha=0.2))[[MethodI]], border=NA )
      }
    }
  dev.off() 
}

    