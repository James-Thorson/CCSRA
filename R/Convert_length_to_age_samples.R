
#' Convert length-composition to age-composition data
#'
#' \code{Convert_length_to_age_samples} converts a matrix of length-composition data to age-composition data for use in CCSRA
#'
#' @param K, Brody growth coefficient
#' @param Linf, Asymptotic maximum length
#' @param L0, Expected length at age 0 from fitting a von Bertalanffy growth curve
#' @param Lcv, coefficient of variation for size given expected length at age
#' @param Lbin_mat, a matrix with two columns, where each row specifies the lower and upper length for a given length bin (the lowest should be -Inf, and highest Inf)
#' @param LengthComp_lt, a matrix of length-composition data, where columns are samples in year t, and cells are the count of samples with a given length and year
#' @param checkforbugs, a boolean stating whether to check for bugs or not

#' @return AgeComp_at, a matrix of expected age-composition data, where columns are samples in year t, and cells are the count of samples with a given age and year

Convert_length_to_age_samples = function(K, Linf, L0, Lcv, Lbin_mat, LengthComp_lt, checkforbugs=TRUE){
  # warnings
  if( checkforbugs==TRUE ){
    # Bounds
    if( Lbin_mat[1,1]>-Inf ) stop("I recommend using 0 as the lowest length of the first length bin")
    if( Lbin_mat[nrow(Lbin_mat),2]<Inf ) stop("I recommend using Inf as the highest length of the last length bin")
    # Mis-match in inputs
    if( nrow(AgeComp_lt)!=nrow(Lbin_mat) ) stop("Must have as many length bins as rows of data AgeComp_lt")
  }

  # Calculate expected length at age
  L_a = Linf - (Linf-L0)*exp(-K*0:nrow(AgeComp_at))

  # Make key to convert age to length
  K_a2l = matrix(NA, ncol=nrow(AgeComp_at), nrow=nrow(Lbin_mat))
  for( aI in 1:ncol(K_a2l)){
  for( lI in 1:nrow(K_a2l)){
    K_a2l[lI,aI] = pnorm(Lbin_mat[lI,2], mean=L_a[aI], sd=Lcv*L_a[aI]) - pnorm(Lbin_mat[lI,1], mean=L_a[aI], sd=Lcv*L_a[aI])
  }}

  # Make key to convert length to age
  K_l2a = t(K_a2l)
  K_l2a = K_l2a / outer(rep(1,nrow(K_l2a)),colSums(K_l2a))

  # Calculate
  AgeComp_at = K_l2a %*% LengthComp_lt

  # Return
  Return = list("AgeComp_lt"=AgeComp_at)
  return( Return )
}

