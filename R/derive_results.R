
#' Calculate biological reference points
#'
#' \code{derive_results} calculates derived quanties for status or productivity
#'
#' @param Obj, the fitted TMB object

#' @return List, a tagged list of potentially useful benchmarks
derive_results = function( Obj, SD=NULL ){
  # Extract elements
  if( length(Obj$env$random)>0 ){
    par = Obj$env$last.par[-Obj$env$random]
  }else{
    par = Obj$env$last.par
  }
  Data = Obj$env$data
  Report = Obj$report()
  ParHat = Obj$env$parList(par)

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
    if("RecDev_hat"%in%names(ParHat)) ParHat_new[["RecDev_hat"]] = rep(0, Data_new[["Nyears"]]+Data_new[["AgeMax"]])
    if("Rec_par"%in%names(ParHat)) ParHat_new[["Rec_par"]] = rep(0, Data_new[["Nyears"]]+Data_new[["AgeMax"]])
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

  # Bundle
  Return = list("Fmsy"=Fmsy, "SB0"=SB0, "TB0"=TB0, "TB_t"=TB_t, "D_t"=Report$D_t, "SB_t"=Report$SB_t, "MSY"=MSY, "TBmsy"=TBmsy, "SBmsy"=SBmsy)

  # Add standard errors if available
  if( !is.null(SD) ){
    Summ = summary(SD)
    if( "Est. (bias.correct)" %in% colnames(Summ)){
      Summ[,'Estimate'] = ifelse( is.na(Summ[,'Est. (bias.correct)']), Summ[,'Estimate'], Summ[,'Est. (bias.correct)'] )
    }
    Return$ln_D_t = Summ[ which(rownames(Summ)=="ln_D_t"), ]
    Return$ln_SB_t = Summ[ which(rownames(Summ)=="ln_SB_t"), ]
    Return$D_t = Summ[ which(rownames(Summ)=="D_t"), ]
    Return$SB_t = Summ[ which(rownames(Summ)=="SB_t"), ]
    Return$Rec_t = Summ[ which(rownames(Summ)=="Rec_t"), ]
  }

  # Return
  return( Return )
}
