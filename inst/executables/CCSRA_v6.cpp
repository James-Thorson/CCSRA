#include <TMB.hpp>
template<class Type>
Type objective_function<Type>::operator() ()
{
  
  DATA_INTEGER(Nyears);          
  DATA_INTEGER(AgeMax);
  DATA_INTEGER(F_method); // F_method=1: Explicit F; F_method=2: Hybrid method (not implemented)
  DATA_SCALAR( CatchCV );
  DATA_VECTOR( ln_R0_prior );   // 
  DATA_VECTOR( M_prior);
  DATA_VECTOR( h_prior ); 
  DATA_VECTOR( S50_prior );
  DATA_VECTOR( Sslope_prior );
  DATA_VECTOR( F_t_prior );
  DATA_VECTOR( D_prior );
  DATA_VECTOR( SigmaR_prior );
  DATA_VECTOR( RecDev_prior );
  DATA_VECTOR( RecDev_biasadj );
  DATA_VECTOR( Cw_t );
  DATA_VECTOR( W_a );
  DATA_VECTOR( Mat_a );
  DATA_MATRIX( AgeComp_at );
  DATA_MATRIX( Index_t );
   
  // =====================================================================

  PARAMETER( ln_R0 );
  PARAMETER( ln_M );
  PARAMETER( input_h );
  PARAMETER( S50 );
  PARAMETER( Sslope );
  PARAMETER( ln_SigmaR );
  PARAMETER_VECTOR( Survey_par );
  PARAMETER_VECTOR( ln_F_t_input );
  PARAMETER_VECTOR( RecDev_hat );

  // =====================================================================

  Type tmp_sum;
  Type Joint; 

  vector<Type> Param_hat(6);
  vector<Type> F_t(Nyears);
  vector<Type> ln_F_t(Nyears);
  vector<Type> SB_t(Nyears);
  vector<Type> Bexploit_t(Nyears);
  vector<Type> ln_SB_t(Nyears);
  vector<Type> D_t(Nyears);
  vector<Type> ln_D_t(Nyears);
  vector<Type> Cw_t_hat(Nyears);
  vector<Type> S_a(AgeMax+1);
  vector<Type> Rprop_t(Nyears);
  
  matrix<Type> N_at(AgeMax+1,Nyears);
  matrix<Type> Zn_at(AgeMax+1,Nyears);
  matrix<Type> Dn_at(AgeMax+1,Nyears);
  matrix<Type> Cn_at(AgeMax+1,Nyears);
  
  // =====================================================================

  // other global variables
  Type jnll = 0;
  vector<Type> jnll_comp(6);
  // Slot 0:  Catch
  // Slot 1:  Index
  // Slot 2:  Comps
  // Slot 3:  RecDev
  // Slot 4:  F
  // Slot 5:  Priors
  jnll_comp.setZero();
  vector<Type> jnll0_t(Nyears);
  Type pi = 3.141592;
  Type infinity = 1.0 / 0.0;
  
  // Transform parameters
  Type survey_q = exp( Survey_par(0) );
  Type survey_extrasd = exp( Survey_par(1) );
  Type M = exp( ln_M );
  Type h = 0.2001 + 0.7998*1/(1+exp(-input_h));
  Type SigmaR = exp( ln_SigmaR );
  Type R0 = exp( ln_R0 );
  Type SBPR0 = 0;
  for(int AgeI=0; AgeI<=AgeMax; AgeI++){ SBPR0 += exp(-M*Type(AgeI)) * W_a(AgeI) * Mat_a(AgeI); }
  Type SB0 = R0 * SBPR0; 
  for(int YearI=0; YearI<Nyears; YearI++) F_t(YearI) = exp( ln_F_t_input(YearI) );
  
  // Initialization
  // Abundance
  for(int AgeI=0; AgeI<=AgeMax; AgeI++){
    S_a(AgeI) = 1 / (1 + exp( -Sslope * (Type(AgeI) - S50) )); 
    if(RecDev_prior[3]==1) N_at(AgeI,0) = R0 * exp(-M * Type(AgeI)) * exp(RecDev_hat(AgeMax-AgeI) - RecDev_biasadj(AgeMax-AgeI)*pow(SigmaR,2)/2);
  }
  // Calculate F for Pope's approximation
  if(F_method==2){ 
    tmp_sum = 0;
    for(int AgeI=0; AgeI<=AgeMax; AgeI++) tmp_sum += S_a(AgeI) * W_a(AgeI) * N_at(AgeI,0); 
    F_t(0) += Cw_t(0) / (exp(-M/2) * tmp_sum);
    Joint = 1 / (1 + exp( Type(30) * (F_t(0)-Type(0.95)) ));
    F_t(0) =  Joint*F_t(0) + (1-Joint)*0.95;
  }
  // Deaths and removals
  for(int AgeI=0; AgeI<=AgeMax; AgeI++){
    if(F_method==-1 | F_method==1){   // Continuous F
      Zn_at(AgeI,0) = N_at(AgeI,0) * (1 - exp( -M - F_t(0)*S_a(AgeI) )); // Zn = total deaths
      Dn_at(AgeI,0) = Zn_at(AgeI,0) * (M) / (M + F_t(0)*S_a(AgeI));      // Dn = deaths due to natural mortality
      Cn_at(AgeI,0) = Zn_at(AgeI,0) * (F_t(0)*S_a(AgeI)) / (M + F_t(0)*S_a(AgeI));  // Cn = deaths due to fishing
    }
    if(F_method==-2 | F_method==2){   // Instantaneous F at mid-point of year 
      Dn_at(AgeI,0) = N_at(AgeI,0) * (1 - exp(-M/2));
      Cn_at(AgeI,0) = N_at(AgeI,0) * exp(-M/2) * S_a(AgeI)*F_t(0);
      Dn_at(AgeI,0) += N_at(AgeI,0) * exp(-M/2) * (1 - S_a(AgeI)*F_t(0)) * (1 - exp(-M/2));
      Zn_at(AgeI,0) = Cn_at(AgeI,0) + Dn_at(AgeI,0);
    }
  }
  // Summaries
  SB_t(0) = 0;
  Cw_t_hat(0) = 0;
  Bexploit_t(0) = 0; 
  for(int AgeI=0; AgeI<=AgeMax; AgeI++){
    SB_t(0) += W_a(AgeI) * Mat_a(AgeI) * N_at(AgeI,0);
    Cw_t_hat(0) += Cn_at(AgeI,0) * W_a(AgeI);
    Bexploit_t(0) += W_a(AgeI) * S_a(AgeI) * N_at(AgeI,0);
  }
  
  // Projection
  for(int YearI=1; YearI<Nyears; YearI++){
    // Survival
    for(int AgeI=1; AgeI<=AgeMax; AgeI++){
      N_at(AgeI,YearI) = N_at(AgeI-1,YearI-1) - Zn_at(AgeI-1,YearI-1);
    }
    // Spawning and exploitable biomass
    SB_t(YearI) = 0;
    for(int AgeI=1; AgeI<=AgeMax; AgeI++) SB_t(YearI) += W_a(AgeI) * Mat_a(AgeI) * N_at(AgeI,YearI);  // Start at AgeI=1, because recruitment hasn't been calculated yet
    // Recruitment
    if(RecDev_prior[3]==1){
      N_at(0,YearI) = Type(4) * h * R0 * SB_t(YearI) / ( SB0*(Type(1)-h) + SB_t(YearI)*(Type(5)*h-Type(1)) )  * exp(RecDev_hat(AgeMax+YearI) - RecDev_biasadj(AgeMax+YearI)*pow(SigmaR,2)/Type(2));
    }
    // Exploitable biomass
    Bexploit_t(YearI) = 0;
    for(int AgeI=0; AgeI<=AgeMax; AgeI++) Bexploit_t(YearI) += W_a(AgeI) * S_a(AgeI) * N_at(AgeI,YearI);
    // Calculate F for Pope's approximation
    if(F_method==2){ 
      tmp_sum = 0;
      for(int AgeI=0; AgeI<=AgeMax; AgeI++) tmp_sum += S_a(AgeI) * W_a(AgeI) * N_at(AgeI,YearI); 
      F_t(YearI) = Cw_t(YearI) / (exp(-M/2) * tmp_sum);
      Joint = 1 / (1 + exp( Type(30) * (F_t(YearI)-Type(0.95)) ));
      F_t(YearI) =  Joint*F_t(YearI) + (Type(1)-Joint)*Type(0.95);
    }
    // Deaths and removals
    for(int AgeI=0; AgeI<=AgeMax; AgeI++){
      // Removals
      if(F_method==-1 | F_method==1){   // Continuous F
        Zn_at(AgeI,YearI) = N_at(AgeI,YearI) * (Type(1) - exp( -M - F_t(YearI)*S_a(AgeI) ));
        Dn_at(AgeI,YearI) = Zn_at(AgeI,YearI) * (M) / (M + F_t(YearI)*S_a(AgeI));
        Cn_at(AgeI,YearI) = Zn_at(AgeI,YearI) * (F_t(YearI)*S_a(AgeI)) / (M + F_t(YearI)*S_a(AgeI));
      }
      if(F_method==-2 | F_method==2){    // Instantaneous F at mid-point of year
        Dn_at(AgeI,YearI) = N_at(AgeI,YearI) * (Type(1) - exp(-M/2));
        Cn_at(AgeI,YearI) = N_at(AgeI,YearI) * exp(-M/2) * S_a(AgeI)*F_t(YearI);
        Dn_at(AgeI,YearI) += N_at(AgeI,YearI) * exp(-M/2) * (1 - S_a(AgeI)*F_t(YearI)) * (1 - exp(-M/2));
        Zn_at(AgeI,YearI) = Cn_at(AgeI,YearI) + Dn_at(AgeI,YearI);
      }
    }
    // Catch
    Cw_t_hat(YearI) = 0; 
    for(int AgeI=0; AgeI<=AgeMax; AgeI++) Cw_t_hat(YearI) += Cn_at(AgeI,YearI) * W_a(AgeI);
  }
    
  // Reporting 
  ln_F_t = log(F_t);
  ln_SB_t = log(SB_t);
  D_t = SB_t / (SBPR0 * R0);
  ln_D_t = log(D_t);
  Param_hat(0) = ln_R0;
  Param_hat(1) = M;
  Param_hat(2) = h;
  Param_hat(3) = S50;
  Param_hat(4) = Sslope;
  Param_hat(5) = SigmaR;
  vector<Type> Rec_t( Nyears );
  Rec_t = N_at.row(0);

  // Objective function -- catches
  if(F_method==1 | F_method==2){ 
    for(int YearI=0; YearI<Nyears; YearI++){  
      if(Cw_t(YearI)>0) jnll0_t(YearI) = -1 * dnorm( log(Cw_t_hat(YearI)), log(Cw_t(YearI)), CatchCV, true ) - log(Cw_t(YearI));
    }
    jnll_comp(0) = jnll0_t.sum();
  }

  // Objective function -- fishery index
  for(int YearI=0; YearI<Nyears; YearI++){
    //if( Index_t(YearI,0)!=infinity ){
    if( !R_IsNA(asDouble(Index_t(YearI,0))) ){
      jnll_comp(1) -= dnorm( log(Index_t(YearI,0)), log(survey_q)+log(Bexploit_t(YearI)), Index_t(YearI,1)+survey_extrasd, true);
    }
  }
  
  // Objective function -- compositional data
  for(int AgeI=0; AgeI<=AgeMax; AgeI++){  
  for(int YearI=0; YearI<Nyears; YearI++){
    tmp_sum = 0;
    for(int AgeII=0; AgeII<=AgeMax; AgeII++) tmp_sum += Cn_at(AgeII,YearI);
    jnll_comp(2) -= AgeComp_at(AgeI,YearI) * log( Cn_at(AgeI,YearI)/tmp_sum*0.9999 + 0.0001/(AgeMax+1) );
  }}
  
  // RecDevs
  for(int Index=0; Index<(Nyears+AgeMax); Index++){
    if( SigmaR>0 ) jnll_comp(3) -= dnorm( RecDev_hat(Index), Type(0), SigmaR, true);
  }

  // F_t
  vector<Type> pen_F(Nyears);
  for(int YearI=0; YearI<Nyears; YearI++){
    pen_F(YearI) = 0;
    if( true ){
      if( log(F_t(YearI)) < log(F_t_prior(0)) ) pen_F(YearI) = ( dnorm( log(F_t(YearI)), log(F_t_prior(0)), F_t_prior(4), true ) - dnorm( log(F_t_prior(0)), log(F_t_prior(0)), F_t_prior(4), true ) );
      if( log(F_t(YearI)) > log(F_t_prior(1)) ) pen_F(YearI) = ( dnorm( log(F_t(YearI)), log(F_t_prior(1)), F_t_prior(4), true ) - dnorm( log(F_t_prior(1)), log(F_t_prior(1)), F_t_prior(4), true ) );
    }
  }
  jnll_comp(4) -= sum(pen_F);

  // Objective function -- priors
  // M: lognormal
  jnll_comp(5) -= ( dnorm( log(M), log(M_prior(3)), M_prior(4), true) - log(M_prior(3)) );
  // h: some kinda beta, I forget
  jnll_comp(5) -= ( (h_prior(3)-1)*log((h-h_prior(0))/(h_prior(1)-h_prior(0))) + (h_prior(4)-1)*log(1-(h-h_prior(0))/(h_prior(1)-h_prior(0))) ); // h
  // Depletion
  if(D_prior(2)==1) jnll_comp(5) -= ( dnorm( log(D_t(Nyears-1)), log(D_prior(0)), D_prior(1), true) - log(D_prior(0)) );
  // SigmaR
  if(SigmaR_prior(3)>0) jnll_comp(5) -= ( dnorm( log(SigmaR), log(SigmaR_prior(3)), SigmaR_prior(4), true) - log(SigmaR_prior(3)) );
  // Selex-S50
  jnll_comp(5) -= dnorm( S50, S50_prior(3), S50_prior(4), true);
  // Selex-slope
  // if bounds aren't equal, Sslope is normal past bounds (slots 1 and 2) with input SD (slot 5)
  if( Sslope_prior(0) < Sslope_prior(1) ){
    if( Sslope < Sslope_prior(0) ) jnll_comp(5) -= ( dnorm( Sslope, Sslope_prior(0), Sslope_prior(4), true ) - dnorm( Sslope_prior(0), Sslope_prior(0), Sslope_prior(4), true ) );
    if( Sslope > Sslope_prior(1) ) jnll_comp(5) -= ( dnorm( Sslope, Sslope_prior(1), Sslope_prior(4), true ) - dnorm( Sslope_prior(1), Sslope_prior(1), Sslope_prior(4), true ) );
  }
  // if bounds are equal, Sslope is normal with input SD (slot 5) and input mean (slot 1)
  if( Sslope_prior(0) >= Sslope_prior(1) ){
    jnll_comp(5) -= dnorm( Sslope, (Sslope_prior(0)+Sslope_prior(1))/2, Sslope_prior(4), true );
  }

  vector<Type> RecMult_t(Nyears+AgeMax);
  for(int t=0; t<(Nyears+AgeMax); t++){
    RecMult_t(t) = exp(RecDev_hat(t) - RecDev_biasadj(t)*pow(SigmaR,2)/Type(2));
  }
  
  // Total objective function
  jnll = jnll_comp.sum();

  // ===========================================================================
  REPORT(RecMult_t);
  REPORT(S_a);
  REPORT(D_t);
  REPORT(F_t);
  REPORT(Param_hat);
  REPORT(SB_t);
  REPORT(N_at);
  REPORT(ln_F_t);
  REPORT(Cw_t_hat);
  REPORT(ln_SB_t);
  REPORT(ln_D_t);
  REPORT(RecDev_hat);
  REPORT(SigmaR);
  REPORT(jnll);
  REPORT(pen_F);
  REPORT( survey_q );
  REPORT( survey_extrasd );
  REPORT( Bexploit_t );
  REPORT( jnll_comp );
  REPORT( jnll0_t );
  REPORT( CatchCV );
  REPORT( Cn_at );
  REPORT( Rec_t );

  //ADREPORT(S_a);
  ADREPORT(Param_hat);
  //ADREPORT(F_t);
  ADREPORT(SB_t);
  ADREPORT(D_t);
  ADREPORT(ln_F_t);
  ADREPORT(ln_SB_t);
  ADREPORT(ln_D_t);  
  //ADREPORT(SigmaR);
  ADREPORT(RecMult_t);
  ADREPORT( Rec_t );
  // ===========================================================================
  return jnll;
}

