\name{SimData_Fn}
\alias{SimData_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for simulating data to test CCSRA
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SimData_Fn(Nyears, AgeMax, SigmaR, M, F1, W_a, S_a, Mat_a, h, SB0, Frate, Fequil, SigmaF, Ncomp_per_year)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Nyears}{
%%     ~~Describe \code{Nyears} here~~
}
  \item{AgeMax}{
%%     ~~Describe \code{AgeMax} here~~
}
  \item{SigmaR}{
%%     ~~Describe \code{SigmaR} here~~
}
  \item{M}{
%%     ~~Describe \code{M} here~~
}
  \item{F1}{
%%     ~~Describe \code{F1} here~~
}
  \item{W_a}{
%%     ~~Describe \code{W_a} here~~
}
  \item{S_a}{
%%     ~~Describe \code{S_a} here~~
}
  \item{Mat_a}{
%%     ~~Describe \code{Mat_a} here~~
}
  \item{h}{
%%     ~~Describe \code{h} here~~
}
  \item{SB0}{
%%     ~~Describe \code{SB0} here~~
}
  \item{Frate}{
%%     ~~Describe \code{Frate} here~~
}
  \item{Fequil}{
%%     ~~Describe \code{Fequil} here~~
}
  \item{SigmaF}{
%%     ~~Describe \code{SigmaF} here~~
}
  \item{Ncomp_per_year}{
%%     ~~Describe \code{Ncomp_per_year} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Nyears, AgeMax, SigmaR, M, F1, W_a, S_a, Mat_a, h, 
    SB0, Frate, Fequil, SigmaF, Ncomp_per_year) 
{
    Cw_t = SB_t = F_t = rep(NA, Nyears)
    Zn_at = Dn_at = Cn_at = N_at = matrix(NA, nrow = AgeMax + 
        1, ncol = Nyears)
    RecDev = rnorm(Nyears + AgeMax, mean = -SigmaR^2/2, sd = SigmaR)
    F_t[1] = F1
    N_at[, 1] = R0 * exp(-M * 0:AgeMax) * exp(RecDev[1:(AgeMax + 
        1)])
    SB_t[1] = sum(N_at[, 1] * W_a * Mat_a)
    if (F_method == -1 | F_method == 1) {
        Zn_at[, 1] = N_at[, 1] * (1 - exp(-M - F_t[1] * S_a))
        Dn_at[, 1] = Zn_at[, 1] * (M)/(M + F_t[1] * S_a)
        Cn_at[, 1] = Zn_at[, 1] * (F_t[1] * S_a)/(M + F_t[1] * 
            S_a)
    }
    if (F_method == -2 | F_method == 2) {
        Dn_at[, 1] = N_at[, 1] * (1 - exp(-M/2))
        Cn_at[, 1] = N_at[, 1] * exp(-M/2) * F_t[1] * S_a
        Dn_at[, 1] = Dn_at[, 1] + N_at[, 1] * exp(-M/2) * (1 - 
            F_t[1] * S_a) * (1 - exp(-M/2))
        Zn_at[, 1] = Dn_at[, 1] + Cn_at[, 1]
    }
    for (YearI in 2:Nyears) {
        F_t[YearI] = F_t[YearI - 1] * (SB_t[YearI - 1]/(Fequil * 
            SB0))^Frate * exp(rnorm(1, mean = -SigmaF^2/2, sd = SigmaF))
        if (F_t[YearI] > 0.95 & (F_method == -2 | F_method == 
            2)) 
            F_t[YearI] = 0.95
        N_at[-1, YearI] = N_at[-AgeMax, YearI - 1] - Zn_at[-AgeMax, 
            YearI - 1]
        SB_t[YearI] = sum((N_at[, YearI] * W_a * Mat_a)[-1])
        N_at[1, YearI] = 4 * h * R0 * SB_t[YearI]/(SB0 * (1 - 
            h) + SB_t[YearI] * (5 * h - 1)) * exp(RecDev[AgeMax + 
            YearI])
        if (F_method == -1 | F_method == 1) {
            Zn_at[, YearI] = N_at[, YearI] * (1 - exp(-M - F_t[YearI] * 
                S_a))
            Dn_at[, YearI] = Zn_at[, YearI] * (M)/(M + F_t[YearI] * 
                S_a)
            Cn_at[, YearI] = Zn_at[, YearI] * (F_t[YearI] * S_a)/(M + 
                F_t[YearI] * S_a)
        }
        if (F_method == -2 | F_method == 2) {
            Dn_at[, YearI] = N_at[, YearI] * (1 - exp(-M/2))
            Cn_at[, YearI] = N_at[, YearI] * exp(-M/2) * F_t[YearI] * 
                S_a
            Dn_at[, YearI] = Dn_at[, YearI] + N_at[, YearI] * 
                exp(-M/2) * (1 - F_t[YearI] * S_a) * (1 - exp(-M/2))
            Zn_at[, YearI] = Dn_at[, YearI] + Cn_at[, YearI]
        }
    }
    Cw_t = (W_a \%*\% Cn_at)[1, ]
    AgeComp_at = array(0, dim = dim(N_at))
    for (YearI in 2:Nyears) {
        AgeComp_at[, YearI] = rmultinom(n = 1, size = Ncomp_per_year, 
            prob = Cn_at[, YearI])[, 1]
    }
    DataList = list(Cw_t = Cw_t, SB_t = SB_t, F_t = F_t, N_at = N_at, 
        AgeComp_at = AgeComp_at, RecDev = RecDev)
    return(DataList)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
