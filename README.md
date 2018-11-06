Description
=============

CCSRA
* Software for implementing catch curve stock reduction analysis (CCSRA; Thorson and Cope 2014 Fish Res) implemented in Template Model Builder


Instructions
=============
First, please install TMB (Template Model Builder) here: 
https://github.com/kaskr/adcomp

Next, please use R version >=3.1.1 and install the package:


    # Install package
    install.packages("devtools")
    library("devtools")
    install_github("James-Thorson/CCSRA")
    # Load package
    library(CCSRA)

Please see examples folder for an example of how to run the model:
https://github.com/James-Thorson/CCSRA/blob/master/examples/CCRA_2014-10-30.R

Known installation/usage issues
=============
none

Further reading
=============

For more details regarding development and testing of this delta-GLMM software please see:
* Thorson, J.T., Rudd, M.B., Winker, H., 2018. The case for estimating recruitment variation in data-moderate and data-poor age-structured models. Fish. Res. https://doi.org/10.1016/j.fishres.2018.07.007
* Thorson, J.T., Kristensen, K., 2016. Implementing a generic method for bias correction in statistical models using random effects, with spatial and population dynamics examples. Fish. Res. 175, 66–74. https://doi.org/10.1016/j.fishres.2015.11.016
* Thorson, J.T., and Cope, J.M. 2015. Catch curve stock-reduction analysis: an alternative solution to the catch equation. Fish. Res. 171: 33–41. http://www.sciencedirect.com/science/article/pii/S0165783614001507
* Thorson, J.T., and Prager, M.H. 2011. Better catch curves: Incorporating age-specific natural mortality and logistic selectivity. Trans. Am. Fish. Soc. 140: 356–366. 

