FormatInput_Fn <-
function( Method, M_prior, h_prior, D_prior, SigmaR_prior, AgeComp_at, Cw_t, W_a, Mat_a, RecDev_biasadj ){
  
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
    Map[["ln_SigmaR"]] = factor(NA)
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
    Map[["ln_SigmaR"]] = factor(NA)
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
