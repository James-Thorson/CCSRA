\name{FormatInput_Fn}
\alias{FormatInput_Fn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Function for generating inputs for parameter estimation using Template Model Builder
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
FormatInput_Fn(Method, M_prior, h_prior, D_prior, SigmaR_prior, AgeComp_at, Cw_t, W_a, Mat_a, RecDev_biasadj)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Method}{
%%     ~~Describe \code{Method} here~~
}
  \item{M_prior}{
%%     ~~Describe \code{M_prior} here~~
}
  \item{h_prior}{
%%     ~~Describe \code{h_prior} here~~
}
  \item{D_prior}{
%%     ~~Describe \code{D_prior} here~~
}
  \item{SigmaR_prior}{
%%     ~~Describe \code{SigmaR_prior} here~~
}
  \item{AgeComp_at}{
%%     ~~Describe \code{AgeComp_at} here~~
}
  \item{Cw_t}{
%%     ~~Describe \code{Cw_t} here~~
}
  \item{W_a}{
%%     ~~Describe \code{W_a} here~~
}
  \item{Mat_a}{
%%     ~~Describe \code{Mat_a} here~~
}
  \item{RecDev_biasadj}{
%%     ~~Describe \code{RecDev_biasadj} here~~
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
function (Method, M_prior, h_prior, D_prior, SigmaR_prior, AgeComp_at, 
    Cw_t, W_a, Mat_a, RecDev_biasadj) 
{
    Nyears = ncol(AgeComp_at)
    MaxAge = nrow(AgeComp_at) - 1
    h_alpha = ((h_prior[1] - 0.2)/0.8) * (((h_prior[1] - 0.2)/0.8) * 
        (1 - ((h_prior[1] - 0.2)/0.8))/h_prior[2]^2 - 1)
    h_beta = (1 - ((h_prior[1] - 0.2)/0.8)) * (((h_prior[1] - 
        0.2)/0.8) * (1 - ((h_prior[1] - 0.2)/0.8))/h_prior[2]^2 - 
        1)
    ln_R0_prior = c(10, 30, 20, 999, 999, 1)
    F_t_prior = c(0, 3, 0.1, 999, 999, 1, Nyears)
    h_prior = c(0.2, 1, ifelse(Method == "CC", 0.9999, 0.8), 
        h_alpha, h_beta, 1)
    M_prior = c(0, 1, 0.2, M_prior[1], M_prior[2], 4)
    S50_prior = c(999, 999, 5, 999, 999, 3)
    Sslope_prior = c(999, 999, ifelse(Method == "SRA", 10, 1), 
        999, 999, 1)
    D_prior = c(D_prior[1], D_prior[2], ifelse(Method == "SRA", 
        1, 0))
    SigmaR_prior = c(0, 1, 0.6, SigmaR_prior[1], SigmaR_prior[2], 
        -1)
    RecDev_prior = c(-3, 3, 0, 999, 999, 5)
    CatchCV = 0.01
    Data = list(Nyears = Nyears, AgeMax = AgeMax, F_method = F_method, 
        CatchCV = CatchCV, ln_R0_prior = ln_R0_prior, M_prior = M_prior, 
        h_prior = h_prior, S50_prior = S50_prior, Sslope_prior = Sslope_prior, 
        F_t_prior = F_t_prior, D_prior = D_prior, SigmaR_prior = SigmaR_prior, 
        RecDev_prior = RecDev_prior, RecDev_biasadj = RecDev_biasadj, 
        Cw_t = Cw_t, W_a = W_a, Mat_a = Mat_a, AgeComp_at = AgeComp_at)
    if (Method == "SRA") 
        Data$AgeComp_at[] = 0
    Parameters = list(ln_R0 = ln_R0_prior[3], M = M_prior[3], 
        h = h_prior[3], S50 = S50_prior[3] + rnorm(1), Sslope = Sslope_prior[3] + 
            rnorm(1), ln_SigmaR = log(SigmaR_prior[3]), ln_F_t_input = log(rep(0.1, 
            Nyears)), RecDev_hat = rep(0, AgeMax + Nyears))
    Map = list()
    if (Method == "CC") {
        Map[["ln_F_t_input"]] = factor(rep(1, length(Parameters$ln_F_t_input)))
        Map[["ln_R0"]] = factor(NA)
    }
    if (Method == "SRA") {
        Map[["ln_SigmaR"]] = factor(NA)
        Map[["S50"]] = factor(NA)
        Map[["Sslope"]] = factor(NA)
        Map[["ln_F_t_input"]] = factor(rep(NA, length(Parameters$ln_F_t_input)))
    }
    if (Method == "CCSRA") {
        Map[["ln_F_t_input"]] = factor(rep(NA, length(Parameters$ln_F_t_input)))
    }
    Random = NULL
    if (!("RecDev_hat" \%in\% names(Map))) 
        Random = c("RecDev_hat")
    if (FALSE) {
        setwd(AdmbFile)
        dyn.unload(paste0(AdmbFile, dynlib(Version)))
        file.remove(paste0(Version, c(".dll", ".o")))
        compile(paste0(Version, ".cpp"))
    }
    InputList = list(Data = Data, Parameters = Parameters, Random = Random, 
        Map = Map)
    return(InputList)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
