

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

# Date file
Date = paste0(Sys.Date(),"b")
  DateFile = paste0(getwd(),"/",Date,"/")
  dir.create(DateFile)
FigFile = paste0(DateFile,"Figs/")
  dir.create(FigFile)

# Compile model
Version = "CCSRA_v5"   # v3: Added priors on h and M; v4: Added RecDevs; v5: fixed bug in steepness bounding
setwd(TmbFile)
#dyn.unload( paste0(Version,".dll") )
#file.remove( paste(Version,c(".dll",".o"),sep="") )
compile( paste0(Version,".cpp") )

#######################
# Settings
#######################

# General
AgeMax = 20
Nyears = 20
Ncomp_per_year = 100
SurveyCV = 0.4
MethodSet = c("CC", "CCSRA", "SRA", "AS" )[4] # 1: Catch curve; 2: CC-SRA; 3:DB-SRA; 4: Age-structured
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

#######################
# Simulate data
#######################

# Which model?
MethodI = 1
Method = MethodSet[MethodI]

# Simulation settings
F_method = switch(Method, "CC"=-1, "CCSRA"=1, "SRA"=1, "AS"=1) # 1: Explicit F; 2: Hybrid (not implemented)

# Simulate data
DataList = SimData_Fn( F_method=F_method, Nyears=Nyears, AgeMax=AgeMax, SigmaR=SigmaR, M=M, F1=F1, W_a=W_a, S_a=S_a, Mat_a=Mat_a, h=h, SB0=SB0, Frate=Frate, Fequil=Fequil, SigmaF=SigmaF, Ncomp_per_year=Ncomp_per_year, SurveyCV=SurveyCV )
# Exclude all age-comp except for final year, except for age-structured model
if( Method%in%c("CC","CCSRA") ) DataList[['AgeComp_at']][,1:(Nyears-1)] = 0
if( Method=="SRA" ) DataList[['AgeComp_at']][] = 0
# Turn off index except for age-structured model
if( Method%in%c("CC","CCSRA","SRA") ) DataList[['Index_t']][,1] = NA

# Plot time series
matplot( cbind(DataList[['SB_t']]/SB0,DataList[['F_t']],DataList[['Cw_t']]/max(DataList[['Cw_t']])), type="l", col=c("black","red","blue"), lty="solid")

#######################
# Estimate model
#######################

# Fit twice for bias adjustment if estimating recruitment deviations
for(LoopI in 1:2){
  # Bias adjustment for each loop
  if(LoopI==2 && !("condition" %in% names(attributes(Sdreport))) ){
    SD = summary(Sdreport)      
    RecDev_biasadj = 1 - SD[which(rownames(SD)=="RecDev_hat"),'Std. Error']^2 / Report$SigmaR^2
  }else{
    RecDev_biasadj = rep(1, Nyears+AgeMax)
  } 
 
  # Format inputs
  InputList = FormatInput_Fn(Method=Method, M_prior=c(M,0.5), h_prior=c(h,0.1), Sslope_prior=c(1,1,1), D_prior=c(0.4,0.2), SigmaR_prior=c(0.6,0.2), AgeComp_at=DataList[['AgeComp_at']], Index_t=DataList[['Index_t']], Cw_t=DataList[['Cw_t']], W_a=W_a, Mat_a=Mat_a, RecDev_biasadj=RecDev_biasadj)
  
  # Compile 
  dyn.load( paste0(TmbFile,dynlib(Version)) )
  if(LoopI==1) Obj <- MakeADFun(data=InputList[['Data']], parameters=InputList[['Parameters']], map=InputList[['Map']], random=InputList[['Random']], inner.control=list(maxit=1e3) )
  if(LoopI==2) Obj <- MakeADFun(data=InputList[['Data']], parameters=ParList, map=InputList[['Map']], random=InputList[['Random']], inner.control=list(maxit=1e3) )
  InitVal = Obj$fn( Obj$par )
  
  # Check for bad start
  if( is.nan(InitVal) & LoopI==2 ){
    Obj <- MakeADFun(data=InputList[['Data']], parameters=InputList[['Parameters']], map=InputList[['Map']], random=InputList[['Random']], inner.control=list(maxit=1e3) )
    InitVal = Obj$fn( Obj$par )
  }
  
  # Set bounds
  Obj$env$inner.control <- c(Obj$env$inner.control, "step.tol"=c(1e-8,1e-12)[1], "tol10"=c(1e-6,1e-8)[1], "grad.tol"=c(1e-8,1e-12)[1]) 
  Upr = rep(Inf, length(Obj$par))
    Upr[match("ln_SigmaR",names(Obj$par))] = log(2)
    Upr[match("ln_F_t_input",names(Obj$par))] = log(2)
    Upr[match("Sslope",names(Obj$par))] = log(5)
    Upr[match("S50",names(Obj$par))] = AgeMax*1.5
  Lwr = rep(-Inf, length(Obj$par))
  
  # Run
  Opt = nlminb( start=Obj$par, objective=Obj$fn, gradient=Obj$gr, upper=Upr, lower=Lwr, control=list(trace=1, eval.max=1e4, iter.max=1e4, rel.tol=1e-10) )
  Opt[["final_gradient"]] = Obj$gr( Opt$par )
  
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
  Sdreport = try( sdreport(Obj) )
  ParList = Obj$env$parList( Opt$par )

  # Derived quantities
  Derived = Calc_derived_quants( Obj )
}  # End fitting loop

    