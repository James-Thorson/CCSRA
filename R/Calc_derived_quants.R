
#' Calculate biological reference points
#'
#' \code{Calc_derived_quants} calculates derived quanties for status or productivity
#'
#' @param Obj, the fitted TMB object

#' @return List, a tagged list of potentially useful benchmarks



Calc_derived_quants = function( Obj ){
  # Extract elements
  Data = Obj$env$data
  ParHat = Obj$env$parList()
  Report = Obj$report()

  # Total biomass
  TB_t = as.vector( Data$W_a %*% Report$N_at )

  # MSY calculations
  Yield_Fn = function( Fmean, Return_type="Yield" ){
    # Modify data
    Data_new = Data
    Data_new[["Nyears"]] = 1000
    Data_new[["RecDev_biasadj"]] = rep(0, Data_new[["Nyears"]]+Data_new[["AgeMax"]])
    Data_new[["Cw_t"]] = rep(1, Data_new[["Nyears"]])
    Data_new[["AgeComp_at"]] = matrix(0, nrow=Data_new[["AgeMax"]]+1, ncol=Data_new[["Nyears"]])
    Data_new[["Index_t"]] = cbind( NA, rep(1,Data_new[["Nyears"]]))
    # Modify parameters
    ParHat_new = ParHat
    ParHat_new[["ln_F_t_input"]] = rep( log(Fmean+1e-10), Data_new[["Nyears"]])
    ParHat_new[["RecDev_hat"]] = rep(0, Data_new[["Nyears"]]+Data_new[["AgeMax"]])
    Obj_new = MakeADFun(data=Data_new, parameters=ParHat_new, inner.control=list(maxit=1e3) )
    # Extract and return stuff
    Report_new = Obj_new$report()
    if(Return_type=="Yield") Return = rev(Report_new$Cw_t_hat)[1]
    if(Return_type=="Report") Return = Report_new
    return(Return)
  }

  # Calculate Fmsy
  Fmsy = optimize( f=Yield_Fn, interval=c(0,3), maximum=TRUE)$maximum
  Report_msy = Yield_Fn( Fmean=Fmsy, Return_type="Report" )
  Report_0 = Yield_Fn( Fmean=0, Return_type="Report" )
  TBmsy = rev(as.vector(Data$W_a %*% Report_msy$N_at))[1]
  SBmsy = rev(Report_msy$SB_t)[1]
  MSY = rev(Report_msy$Cw_t_hat)[1]
  TB0 = rev(as.vector(Data$W_a %*% Report_0$N_at))[1]
  SB0 = rev(Report_0$SB_t)[1]

  # Return
  Return = list("Fmsy"=Fmsy, "SB0"=SB0, "TB0"=TB0, "TB_t"=TB_t, "SB_t"=Report$SB_t, "MSY"=MSY, "TBmsy"=TBmsy, "SBmsy"=SBmsy)
  return( Return )
}
