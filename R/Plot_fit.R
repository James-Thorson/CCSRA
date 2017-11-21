plot_fit = function( Obj, plotdir=getwd() ){
  Report = Obj$report()
  Data = Obj$env$data

  # Time series fits
  ThorsonUtilities::save_fig( file=paste0(plotdir,"/Timeseries"), width=4, height=4 )
    par( mfrow=c(2,1), mar=c(2,2,1,0), mgp=c(2,0.5,0), tck=-0.02 )
    # Catches
    plot( x=1:Data$Nyears, y=Data$Cw_t, main="Catch (kg)" )
    lines( x=1:Data$Nyears, y=Report$Cw_t_hat, col="red", lwd=1.5 )
    # Abundance
    plot( x=1:Data$Nyears, y=Data$Index_t[,1], ylim=c(0,max(c(1,Data$Index_t[,1]),na.rm=TRUE)), main="Index (kg)" )
    lines( x=1:Data$Nyears, y=Report$survey_q*Report$Bexploit_t, col="red", lwd=1.5 )
  dev.off()

  # Abundance at age
  Dim = c( ceiling(sqrt(Data$Nyears)), ceiling(Data$Nyears/ceiling(sqrt(Data$Nyears))) )
  ThorsonUtilities::save_fig( file=paste0(plotdir,"/Abundance-at-age"), width=Dim[2]*2, height=Dim[1]*2 )
    par( mfrow=Dim, mar=c(0,2,1,0), mgp=c(2,0.5,0), tck=-0.02, oma=c(2,0,0,0) )
    # Catches
    for( tI in 1:Data$Nyears ){
      plot( x=0:Data$AgeMax, y=Data$AgeComp_at[,tI]/max(1,sum(Data$AgeComp_at[,tI])), xaxt="n", cex=2 )
      lines( x=0:Data$AgeMax, y=Report$Cn_at[,tI]/sum(Report$Cn_at[,tI]), col="red", lwd=1.5)
      if( tI > (Nyears-Dim[2]) ) axis(1)
    }
  dev.off()
}
