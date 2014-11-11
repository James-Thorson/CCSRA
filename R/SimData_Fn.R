SimData_Fn <-
function( Nyears, AgeMax, SigmaR, M, F1, W_a, S_a, Mat_a, h, SB0, Frate, Fequil, SigmaF, Ncomp_per_year, SurveyCV ){
  # Data objects
  Cw_t = SB_t = F_t = Bexploit_t = rep(NA, Nyears)
  Zn_at = Dn_at = Cn_at = N_at = matrix(NA, nrow=AgeMax+1, ncol=Nyears)
  RecDev = rnorm( Nyears + AgeMax, mean=-SigmaR^2/2, sd=SigmaR )
  
  # Initialization
  F_t[1] = F1
  N_at[,1] = R0 * exp(-M * 0:AgeMax) * exp( RecDev[(AgeMax+1):1] )
  SB_t[1] = sum( N_at[,1] * W_a * Mat_a )
  Bexploit_t[1] = sum( N_at[,1] * W_a * S_a )
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
    # Exploitable biomass
    Bexploit_t[YearI] = sum( N_at[,YearI] * W_a * S_a )
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

  # Generate compositional data
  AgeComp_at = array(0, dim=dim(N_at))
  for(YearI in 2:Nyears){
    AgeComp_at[,YearI] = rmultinom(n=1, size=Ncomp_per_year, prob=Cn_at[,YearI])[,1]
  }

  # Generate index
  Index_t = cbind( Bexploit_t/mean(Bexploit_t) * exp(rnorm(Nyears, mean=0, sd=SurveyCV)), SurveyCV)

  # List
  DataList = list("Cw_t"=Cw_t, "SB_t"=SB_t, "F_t"=F_t, "N_at"=N_at, "AgeComp_at"=AgeComp_at, "Index_t"=Index_t, "RecDev"=RecDev, "Bexploit_t"=Bexploit_t)
  return(DataList)
}
