
#' Build inputs for CCSRA
#'
#' \code{make_inputs} builds inputs necessary for running CCSRA
#'
#' @param composition_likelihood Integer specifying the likelihood for compositional data, from the available options below:
#' \describe{
#'   \item{\code{composition_likelihood=0}}{Multinomial}
#'   \item{\code{composition_likelihood=1}}{Dirichlet-multinomial using linear parameterization}
#'   \item{\code{composition_likelihood=2}}{Multivariate-Tweedie using fixed power parameter \code{p = 1.2}}
#'   \item{\code{composition_likelihood=3}}{Multivariate-Tweedie estimating both parameters phi and p}
#' }
#'
#' @return a list with inputs for building CCSRA TMB object
#'
#' @export
make_inputs <-
function( Version = "CCSRA_v9",
          Method,
          M_prior,
          h_prior,
          D_prior,
          SigmaR_prior,
          Sslope_prior = c(-999,999,1),
          AgeComp_at,
          Index_t,
          Cw_t,
          W_a,
          Mat_a,
          RecDev_biasadj,
          F_method,
          rec_method = "dev",
          estimate_recdevs = TRUE,
          composition_likelihood = 0 ){
  
  # Generate defaults
  if( missing(F_method) ){
    F_method = switch(Method, "CC"=-1, "CCSRA"=1, "SRA"=1, "AS"=1, "ASSP"=1 ) # 1: Explicit F; 2: Hybrid (not implemented)
  }

  # Calculate derived stuff
  Nyears = ncol(AgeComp_at)
  MaxAge = nrow(AgeComp_at)-1 # includes age 0
  
  # Calculate beta prior parameters for h
  h_alpha = ((h_prior[1]-0.2)/0.8) * ( ((h_prior[1]-0.2)/0.8) * (1 - ((h_prior[1]-0.2)/0.8)) / h_prior[2]^2 -1)
  h_beta = (1 - ((h_prior[1]-0.2)/0.8)) * (((h_prior[1]-0.2)/0.8) * (1-((h_prior[1]-0.2)/0.8)) / h_prior[2]^2 - 1)
  
  # Priors
  ln_R0_prior = c(10, 30, 20, 999, 999, 1)
  F_t_prior = c( 0, 3, 0.1, 999, 0.1, 1, Nyears )
  h_prior = c(0.2, 1.0, ifelse(Method=="CC",0.9999,h_prior[1]), h_alpha, h_beta, 1)
  M_prior = c(0, 1, M_prior[1], M_prior[1], M_prior[2], 4)
  S50_prior = c(999, 999, 5, 999, 999, 3)
  Sslope_prior = c(Sslope_prior[1], Sslope_prior[2], ifelse(Method=="SRA",10,1), 999, Sslope_prior[3], 1)
  D_prior = c(D_prior[1], D_prior[2], ifelse(Method=="SRA",1,0))
  SigmaR_prior = c( 0, 1, 0.6, SigmaR_prior[1], SigmaR_prior[2], -1)
  RecDev_prior = c( -3, 3, 0, NA, 999, 5 )
  if(rec_method=="dev") RecDev_prior[4] = 1
  if(rec_method=="ln_R0_ratio") RecDev_prior[4] = 2

  # Check for problems
  if( estimate_recdevs==FALSE & rec_method!="dev" ){
    stop("Must use rec-dev parameterization to turn them off")
  }

  # Override defaults
  if( Method %in% c("SRA","ASSP","CC") & composition_likelihood!=0 ){
    message("Over-riding inputs to set `composition_likelihood=0` given choice of `Method`")
    composition_likelihood = 0
  }

  # Compile TMB inputs -- Data
  CatchCV = 0.01
  if(Version %in% c("CCSRA_v9")){
    Data = list("Nyears"=Nyears, "AgeMax"=AgeMax, "F_method"=F_method, "composition_likelihood"=composition_likelihood,
            "CatchCV"=CatchCV, "ln_R0_prior"=ln_R0_prior, "M_prior"=M_prior, "h_prior"=h_prior,
            "S50_prior"=S50_prior, "Sslope_prior"=Sslope_prior, "F_t_prior"=F_t_prior, "D_prior"=D_prior, "SigmaR_prior"=SigmaR_prior,
            "RecDev_prior"=RecDev_prior, "RecDev_biasadj"=RecDev_biasadj, "Cw_t"=Cw_t, "W_a"=W_a, "Mat_a"=Mat_a, "AgeComp_at"=AgeComp_at, "Index_t"=Index_t)
  }
  if(Version %in% c("CCSRA_v8")){
    Data = list("Nyears"=Nyears, "AgeMax"=AgeMax, "F_method"=F_method, "use_dirmult"=composition_likelihood, "CatchCV"=CatchCV, "ln_R0_prior"=ln_R0_prior, "M_prior"=M_prior, "h_prior"=h_prior,
            "S50_prior"=S50_prior, "Sslope_prior"=Sslope_prior, "F_t_prior"=F_t_prior, "D_prior"=D_prior, "SigmaR_prior"=SigmaR_prior,
            "RecDev_prior"=RecDev_prior, "RecDev_biasadj"=RecDev_biasadj, "Cw_t"=Cw_t, "W_a"=W_a, "Mat_a"=Mat_a, "AgeComp_at"=AgeComp_at, "Index_t"=Index_t)
  }
  if(Version %in% c("CCSRA_v7","CCSRA_v6","CCSRA_v5","CCSRA_v4")){
    Data = list("Nyears"=Nyears, "AgeMax"=AgeMax, "F_method"=F_method, "CatchCV"=CatchCV, "ln_R0_prior"=ln_R0_prior, "M_prior"=M_prior, "h_prior"=h_prior,
            "S50_prior"=S50_prior, "Sslope_prior"=Sslope_prior, "F_t_prior"=F_t_prior, "D_prior"=D_prior, "SigmaR_prior"=SigmaR_prior, 
            "RecDev_prior"=RecDev_prior, "RecDev_biasadj"=RecDev_biasadj, "Cw_t"=Cw_t, "W_a"=W_a, "Mat_a"=Mat_a, "AgeComp_at"=AgeComp_at, "Index_t"=Index_t)
  }

  # Compile TMB inputs -- Parameters
  if(Version %in% c("CCSRA_v9")){
    Parameters = list( "ln_R0"=ln_R0_prior[3], "ln_M"=log(M_prior[3]), "input_h"=qlogis((h_prior[3]-0.2)/0.8), "S50"=S50_prior[3],
      "Sslope"=Sslope_prior[3], "ln_SigmaR"=log(SigmaR_prior[4]), "ln_theta"=log(3), "ln_phi"=0, "power_prime"=qlogis(1.2-1),
      "Survey_par"=c(0,log(0.0001)), "ln_F_t_input"=log(rep(0.1,Nyears)), "Rec_par"=rep(0,AgeMax+Nyears))
  }
  if(Version %in% c("CCSRA_v8")){
    Parameters = list( "ln_R0"=ln_R0_prior[3], "ln_M"=log(M_prior[3]), "input_h"=qlogis((h_prior[3]-0.2)/0.8), "S50"=S50_prior[3],
      "Sslope"=Sslope_prior[3], "ln_SigmaR"=log(SigmaR_prior[4]), "ln_theta"=log(3), "Survey_par"=c(0,log(0.0001)),
      "ln_F_t_input"=log(rep(0.1,Nyears)), "Rec_par"=rep(0,AgeMax+Nyears))
  }
  if(Version %in% c("CCSRA_v7")){
    Parameters = list( "ln_R0"=ln_R0_prior[3], "ln_M"=log(M_prior[3]), "input_h"=qlogis((h_prior[3]-0.2)/0.8), "S50"=S50_prior[3],
      "Sslope"=Sslope_prior[3], "ln_SigmaR"=log(SigmaR_prior[4]), "Survey_par"=c(0,log(0.0001)),
      "ln_F_t_input"=log(rep(0.1,Nyears)), "Rec_par"=rep(0,AgeMax+Nyears))
  }
  if(Version %in% c("CCSRA_v6","CCSRA_v5","CCSRA_v4")){
    Parameters = list( "ln_R0"=ln_R0_prior[3], "ln_M"=log(M_prior[3]), "input_h"=qlogis((h_prior[3]-0.2)/0.8), "S50"=S50_prior[3],
      "Sslope"=Sslope_prior[3], "ln_SigmaR"=log(SigmaR_prior[4]), "Survey_par"=c(0,log(0.0001)), "ln_F_t_input"=log(rep(0.1,Nyears)),
      "RecDev_hat"=rep(0,AgeMax+Nyears))
  }

  # Compile TMB inputs -- Turn off parameters      
  Map = list()

  # Method-specific stuff
  if(Method=="CC"){
    Map[["ln_F_t_input"]] = factor( rep(1,length(Parameters$ln_F_t_input)) )
    Map[["ln_R0"]] = factor( NA )
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["input_h"]] = factor(NA)
    Map[["Survey_par"]] = factor( rep(NA,2) )
  }
  if(Method=="SRA"){
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["S50"]] = factor(NA)
    Map[["Sslope"]] = factor(NA)
    Map[["Survey_par"]] = factor( rep(NA,2) )
  }
  if(Method=="CCSRA"){
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["Survey_par"]] = factor( rep(NA,2) )
  }
  if(Method=="AS"){
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["Survey_par"]] = factor( c(1,NA) )
  }
  if(Method=="ASSP"){
    Map[["ln_SigmaR"]] = factor(NA)
    Map[["S50"]] = factor(NA)
    Map[["Sslope"]] = factor(NA)
    Map[["Survey_par"]] = factor( c(1,NA) )
  }

  # Map off overdispersion parameters
  if( composition_likelihood == 0 ){
    if("ln_theta" %in% names(Parameters)) Map[["ln_theta"]] = factor(NA)
    if("ln_phi" %in% names(Parameters)) Map[["ln_phi"]] = factor(NA)
    if("power_prime" %in% names(Parameters)) Map[["power_prime"]] = factor(NA)
  }
  if( composition_likelihood == 1 ){
    if("ln_phi" %in% names(Parameters)) Map[["ln_phi"]] = factor(NA)
    if("power_prime" %in% names(Parameters)) Map[["power_prime"]] = factor(NA)
  }
  if( composition_likelihood == 2 ){
    if("ln_theta" %in% names(Parameters)) Map[["ln_theta"]] = factor(NA)
    if("power_prime" %in% names(Parameters)) Map[["power_prime"]] = factor(NA)
  }
  if( composition_likelihood == 3 ){
    if("ln_theta" %in% names(Parameters)) Map[["ln_theta"]] = factor(NA)
  }

  # declare random
  Random = NULL
  if(Method=="SRA"){
    Random = c(Random, "ln_F_t_input")
  }
  if(Method=="CCSRA"){
    Random = c(Random, "ln_F_t_input")
  }
  if(Method=="AS"){
    Random = c(Random, "ln_F_t_input")
  }
  if(Method=="ASSP"){
    Random = c(Random, "ln_F_t_input")
  }
  if( "RecDev_hat" %in% names(Parameters) ){
    Random = union(Random, "RecDev_hat")
  }
  if( "Rec_par" %in% names(Parameters) ){
    Random = union(Random, "Rec_par")
  }

  # Turn off rec-devs
  if( estimate_recdevs==FALSE ){
    Random = setdiff(Random, c("RecDev_hat","Rec_par") )
    if("Rec_par" %in% names(Parameters)){
      Map[["Rec_par"]] = factor( rep(NA,length(Parameters$Rec_par)) )
    }
    if("RecDev_hat" %in% names(Parameters)){
      Map[["RecDev_hat"]] = factor( rep(NA,length(Parameters$RecDev_hat)) )
    }
    Map[["ln_SigmaR"]] = factor(NA)
    Parameters[["ln_SigmaR"]] = log(0.001)
  }

  ################
  # Drop data for models that dont use it
  ################

  # Exclude all age-comp except for final year for catch curve and CCSRA
  if( Method %in% c("CC","CCSRA") ){
    Data[['AgeComp_at']][,1:(Nyears-1)] = 0
  }
  # Exclude all age-comps for SRA and age-structured production model
  if( Method %in% c("SRA","ASSP") ){
    Data[['AgeComp_at']][] = 0
  }
  # Turn off index except for age-structured model and age-structured production model
  if( Method %in% c("CC","CCSRA","SRA") ){
    Data[['Index_t']][,1] = NA
  }


  # Input
  InputList = list("Data"=Data, "Parameters"=Parameters, "Random"=Random, "Map"=Map)
  return(InputList)
}
