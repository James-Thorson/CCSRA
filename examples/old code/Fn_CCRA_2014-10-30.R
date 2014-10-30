
#Method = c("CC", "CCSRA", "SRA")    
FormatInput_Fn = function( Method, M_prior, h_prior, D_prior, SigmaR_prior, AgeComp_at, Cw_t, W_a, Mat_a, RecDev_biasadj ){
  
  # Calculate derived stuff
  Nyears = ncol(AgeComp_at)
  MaxAge = nrow(AgeComp_at)-1 # includes age 0
  
  # Calculate beta prior parameters for h
  h_alpha = ((h_prior[1]-0.2)/0.8) * ( ((h_prior[1]-0.2)/0.8) * (1 - ((h_prior[1]-0.2)/0.8)) / h_prior[2]^2 -1)
  h_beta = (1 - ((h_prior[1]-0.2)/0.8)) * (((h_prior[1]-0.2)/0.8) * (1-((h_prior[1]-0.2)/0.8)) / h_prior[2]^2 - 1)
  
  # Priors
  ln_R0_prior = c(10, 30, 20, 999, 999, 1)
  F_t_prior = c( 0, 3, 0.1, 999, 999, 1, Nyears )
  h_prior = c(0.2, 1.0, ifelse(Method=="CC",0.9999,0.8), h_alpha, h_beta, 1)
  M_prior = c(0, 1, 0.2, M_prior[1], M_prior[2], 4)
  S50_prior = c(999, 999, 5, 999, 999, 3)
  Sslope_prior = c(999, 999, ifelse(Method=="SRA",10,1), 999, 999, 1)
  D_prior = c(D_prior[1], D_prior[2], ifelse(Method=="SRA",1,0))
  SigmaR_prior = c( 0, 1, 0.6, SigmaR_prior[1], SigmaR_prior[2], -1)
  RecDev_prior = c( -3, 3, 0, 999, 999, 5 )

  # Compile TMB inputs -- Data
  CatchCV = 0.01
  Data = list("Nyears"=Nyears, "AgeMax"=AgeMax, "F_method"=F_method, "CatchCV"=CatchCV, "ln_R0_prior"=ln_R0_prior, "M_prior"=M_prior, "h_prior"=h_prior, 
            "S50_prior"=S50_prior, "Sslope_prior"=Sslope_prior, "F_t_prior"=F_t_prior, "D_prior"=D_prior, "SigmaR_prior"=SigmaR_prior, 
            "RecDev_prior"=RecDev_prior, "RecDev_biasadj"=RecDev_biasadj, "Cw_t"=Cw_t, "W_a"=W_a, "Mat_a"=Mat_a, "AgeComp_at"=AgeComp_at)
  if(Method=="SRA") Data$AgeComp_at[] = 0

  # Compile TMB inputs -- Parameters
  Parameters = list( "ln_R0"=ln_R0_prior[3], "M"=M_prior[3], "h"=h_prior[3], "S50"=S50_prior[3]+rnorm(1), "Sslope"=Sslope_prior[3]+rnorm(1), "ln_SigmaR"=log(SigmaR_prior[3]), "ln_F_t_input"=log(rep(0.1,Nyears)), "RecDev_hat"=rep(0,AgeMax+Nyears))
  
  # Compile TMB inputs -- Turn off parameters      
  Map = list()
  #if(M_prior[6]<=0) Map[["M"]] = factor(NA)
  #if(h_prior[6]<=0) Map[["h"]] = factor(NA)
  
  # Method-specific stuff
  if(Method=="CC"){
    Map[["ln_F_t_input"]] = factor( rep(1,length(Parameters$ln_F_t_input)) )
    Map[["ln_R0"]] = factor( NA )
  }
  if(Method=="SRA"){
    #Map[["M"]] = factor(NA)
    #Map[["h"]] = factor(NA)
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["S50"]] = factor(NA)
    Map[["Sslope"]] = factor(NA)
    Map[["ln_F_t_input"]] = factor( rep(NA,length(Parameters$ln_F_t_input)) )
  }
  if(Method=="CCSRA"){
    Map[["ln_F_t_input"]] = factor( rep(NA,length(Parameters$ln_F_t_input)) )
  }
  
  #Map[["SigmaR"]] = factor(NA)
  #Map[["RecDev"]] = factor( rep(NA,Nyears+AgeMax) )
  
  # declare random
  Random = NULL
  if( !("RecDev_hat"%in%names(Map)) ) Random = c("RecDev_hat")
  
  # Compile if necessary
  if(FALSE){
    setwd(AdmbFile)
    dyn.unload( paste0(AdmbFile,dynlib(Version)) )
    file.remove( paste0(Version,c(".dll",".o")) )
    compile( paste0(Version,".cpp") )
  }
  
  # Input
  InputList = list("Data"=Data, "Parameters"=Parameters, "Random"=Random, "Map"=Map)
  return(InputList)
}

SimData_Fn = function( Nyears, AgeMax, SigmaR, M, F1, W_a, S_a, Mat_a, h, SB0, Frate, Fequil, SigmaF, Ncomp_per_year ){
  # Data objects
  Cw_t = SB_t = F_t = rep(NA, Nyears)
  Zn_at = Dn_at = Cn_at = N_at = matrix(NA, nrow=AgeMax+1, ncol=Nyears)
  RecDev = rnorm( Nyears + AgeMax, mean=-SigmaR^2/2, sd=SigmaR )
  
  # Initialization
  F_t[1] = F1
  N_at[,1] = R0 * exp(-M * 0:AgeMax) * exp( RecDev[1:(AgeMax+1)] )
  SB_t[1] = sum( N_at[,1] * W_a * Mat_a )
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
    SB_t[YearI] = sum( (N_at[,YearI] * W_a * Mat_a)[-1] )
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

  # Generate data
  AgeComp_at = array(0, dim=dim(N_at))
  for(YearI in 2:Nyears){
    AgeComp_at[,YearI] = rmultinom(n=1, size=Ncomp_per_year, prob=Cn_at[,YearI])[,1]
  }

  # List
  DataList = list("Cw_t"=Cw_t, "SB_t"=SB_t, "F_t"=F_t, "N_at"=N_at, "AgeComp_at"=AgeComp_at, "RecDev"=RecDev)
  return(DataList)
}    
