Demographic study of Centaurea corymbosa using Integral Projection Models by Loïc Pages.
Project supervised by Ophélie Ronce, Eric Imbert and François Rousset.

See my InternshipReport2025_LPages.pdf for details on our work.

1. ModelSelection.Rmd gives a script for model selection of the different vital rates of C. corymbosa using generalized linear mixed models with package spaMM.

2. IPM folder :
   a. IPM_AIC.Rmd gives the discretized mean kernel over all populations and years. (same for IPM_BIC with vital rates selected by BIC);
   b. IPM_AIC_Pop_Year.Rmd gives the discretized kernals for all pair of year and population;
   c. You can plot all the results with the two Rmd files ResultsIPM_*;
   d. Predict.R contains the fitted vital rate function used for the IPM discretization.
