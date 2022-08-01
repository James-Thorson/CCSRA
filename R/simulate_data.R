#' Simulate data for fitting CCSRA
#'
#' \code{simulate_data} simulates data for use in testing CCSRA
#'
#' @export
simulate_data <-
function( Nyears,
          AgeMax,
          SigmaR,
          M,
          F1,
          W_a,
          S_a,
          Mat_a,
          h,
          SB0,
          R0,
          Frate,
          Fequil,
          SigmaF,
          F_method = 1,
          Fdynamics = "Endogenous",
          Fmax = NA,
          Ncomp_per_year,
          SurveyCV,
          Recruitment_dynamics = "Stationary",
          Regime_multiplier = NA,
          SigmaS = 0,
          rhoA = 0.8,
          rhoT = 0.8 ){

  # Data objects
  Cw_t = SB_t = F_t = Bexploit_t = rep(NA, Nyears)
  Zn_at = Dn_at = Cn_at = N_at = matrix(NA, nrow=AgeMax+1, ncol=Nyears)
  if( Fdynamics=="Ramp" ) Framp_t = c( "rampup"=seq(F1,Fmax,length=floor(Nyears/2)), "peak"=rep(Fmax,floor((Nyears-floor(Nyears/2))/2)), "managed"=rep(Fmax/3,Nyears-floor(Nyears/2)-floor((Nyears-floor(Nyears/2))/2)))

  # Recruitment
  if( Recruitment_dynamics=="Stationary" ){
    Regime_multiplier = rep(0, Nyears+AgeMax)
  }
  if( Recruitment_dynamics=="Regime" ){
    if( is.na(Regime_multiplier) || length(Regime_multiplier)!=(Nyears+AgeMax) ) stop("Please specify `RegimeRmultiplier` for regime shift dynamics")
  }
  RecDev = rnorm(Nyears + AgeMax, mean=0, sd=SigmaR) + Regime_multiplier
  RecMult = exp( RecDev - SigmaR^2/2 )

  # Selectivity
  S_at = Sdev_at = array( NA, dim=c(AgeMax+1,Nyears) )
  for(YearI in 1:Nyears){
    Cov_aa = SigmaS^2 * rhoA ^ as.matrix(dist(0:AgeMax,diag=TRUE,upper=TRUE))
    if( YearI==1 ){
      Sdev_at[,YearI] = mvtnorm::rmvnorm( n=1, mean=rep(0,AgeMax+1), sigma=Cov_aa )[1,]
    }else{
      Sdev_at[,YearI] = rhoT * Sdev_at[,YearI-1] + sqrt(1-rhoT^2) * mvtnorm::rmvnorm( n=1, mean=rep(0,AgeMax+1), sigma=Cov_aa )[1,]
    }
    S_at[,YearI] = S_a * exp(Sdev_at[,YearI])
  }

  # Initialization
  if(Fdynamics=="Endogenous") F_t[1] = F1
  if(Fdynamics=="Ramp") F_t[1] = Framp_t[1]
  N_at[,1] = R0 * exp(-M * 0:AgeMax) * RecMult[(AgeMax+1):1]
  SB_t[1] = sum( N_at[,1] * W_a * Mat_a )
  Bexploit_t[1] = sum( N_at[,1] * W_a * S_at[,1] )
  if(F_method==-1 | F_method==1){
    Zn_at[,1] = N_at[,1] * (1 - exp( -M - F_t[1]*S_at[,1] ))
    Dn_at[,1] = Zn_at[,1] * (M) / (M + F_t[1]*S_at[,1])
    Cn_at[,1] = Zn_at[,1] * (F_t[1]*S_at[,1]) / (M + F_t[1]*S_at[,1])
  }
  if(F_method==-2 | F_method==2){
    Dn_at[,1] = N_at[,1] * (1 - exp(-M/2))
    Cn_at[,1] = N_at[,1] * exp(-M/2) * F_t[1]*S_at[,1]
    Dn_at[,1] = Dn_at[,1] + N_at[,1] * exp(-M/2) * (1 - F_t[1]*S_at[,1]) * (1 - exp(-M/2))
    Zn_at[,1] = Dn_at[,1] + Cn_at[,1]
  }
  
  # Projection
  for(YearI in 2:Nyears){
    # Fishing effort
    if(Fdynamics=="Endogenous"){
      F_t[YearI] = F_t[YearI-1] * (SB_t[YearI-1] / (Fequil * SB0))^Frate * exp( rnorm(1, mean=-SigmaF^2/2, sd=SigmaF) )
      if( F_t[YearI] > 0.95 & (F_method==-2 | F_method==2) ) F_t[YearI] = 0.95
    }
    if(Fdynamics=="Ramp"){
      F_t[YearI] = Framp_t[YearI]
    }
    # Survival
    N_at[-1,YearI] = N_at[-(AgeMax+1),YearI-1] - Zn_at[-(AgeMax+1),YearI-1]
    # Spawning biomass
    SB_t[YearI] = sum( (N_at[,YearI] * W_a * Mat_a)[-1] )
    # Recruitment
    N_at[1,YearI] = 4 * h * R0 * SB_t[YearI] / ( SB0*(1-h) + SB_t[YearI]*(5*h-1) ) * RecMult[AgeMax+YearI]
    # Exploitable biomass
    Bexploit_t[YearI] = sum( N_at[,YearI] * W_a * S_at[,YearI] )
    # Removals
    if(F_method==-1 | F_method==1){
      Zn_at[,YearI] = N_at[,YearI] * (1 - exp( -M - F_t[YearI]*S_at[,YearI] ))
      Dn_at[,YearI] = Zn_at[,YearI] * (M) / (M + F_t[YearI]*S_at[,YearI])
      Cn_at[,YearI] = Zn_at[,YearI] * (F_t[YearI]*S_at[,YearI]) / (M + F_t[YearI]*S_at[,YearI])
    }
    if(F_method==-2 | F_method==2){
      Dn_at[,YearI] = N_at[,YearI] * (1 - exp(-M/2))
      Cn_at[,YearI] = N_at[,YearI] * exp(-M/2) * F_t[YearI]*S_at[,YearI]
      Dn_at[,YearI] = Dn_at[,YearI] + N_at[,YearI] * exp(-M/2) * (1 - F_t[YearI]*S_at[,YearI]) * (1 - exp(-M/2))
      Zn_at[,YearI] = Dn_at[,YearI] + Cn_at[,YearI]
    }
  }
  Cw_t = (W_a %*% Cn_at)[1,]

  # Generate compositional data
  AgeComp_at = array(0, dim=dim(N_at))
  for(YearI in 1:Nyears){
    AgeComp_at[,YearI] = rmultinom(n=1, size=Ncomp_per_year, prob=Cn_at[,YearI])[,1]
  }

  # Generate index
  Index_t = cbind( Bexploit_t * exp(rnorm(Nyears, mean=0, sd=SurveyCV)), SurveyCV )

  # List
  DataList = list("Cw_t"=Cw_t, "SB_t"=SB_t, "F_t"=F_t, "N_at"=N_at, "AgeComp_at"=AgeComp_at, "Index_t"=Index_t,
    "RecDev"=RecDev, "RecMult"=RecMult, "Bexploit_t"=Bexploit_t, "survey_q"=1, "Cn_at"=Cn_at, "Regime_multiplier"=Regime_multiplier,
    "S_at"=S_at, "Sdev_at"=Sdev_at)
  return(DataList)
}
