Code coming from my internship work at ISEM-CNRS on : Demographic study of Centaurea corymbosa using Integral Projection Models. 
By Loïc Pages in 2025.
Project supervised by Ophélie Ronce, Eric Imbert and François Rousset.

See my InternshipReport2025_LPages.pdf for details on our work.

1. Model selection folder : scripts for model selection of the different vital rates of C. corymbosa using generalized linear mixed models with package spaMM.

2. IPM folder :
  a. Predict.R contains the fitted vital rate function used for the IPM discretization.
  b. Kernals implementation : scripts to implement the integral projection model
  c. Results : 
    - Kernal analysis : global results of the IPM (growth rate, size distribution, ...) 
    - Climate : links between life history traits and optimal flowering size with climate
    - Optimal flowering strategy : determine if the flowering strategy is an ESS
    - Validation Flowering data : predictions of the model comparing with the number of observed flowering plants each year
    - Perturbation analysis : elasticity analysis on the vital rates
