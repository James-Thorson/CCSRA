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
* Thorson, J.T., and Cope, J.M. 2015. Catch curve stock-reduction analysis: an alternative solution to the catch equation. Fish. Res. 171: 33–41. http://www.sciencedirect.com/science/article/pii/S0165783614001507
* Thorson, J.T., and Prager, M.H. 2011. Better catch curves: Incorporating age-specific natural mortality and logistic selectivity. Trans. Am. Fish. Soc. 140: 356–366. 

