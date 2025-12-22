Code coming from my internship work at ISEM-CNRS on : Demographic study of Centaurea corymbosa using Integral Projection Models. 
By Loïc Pages in 2025.
Project supervised by Ophélie Ronce, Eric Imbert and François Rousset.

See my InternshipReport2025_LPages.pdf for details on our work.

1. Model selection folder : scripts for model selection of the different vital rates of C. corymbosa using generalized linear mixed models with package spaMM.

2. IPM folder :
  a. Predict.R contains the fitted vital rate function used for the IPM discretization.
  
  b. Kernals implementation : scripts to implement the integral projection model and return kernal matrices
  
  c. Results : 
<<<<<<< HEAD
    - Kernal analysis : global results of the IPM (growth rate, size and age distribution, distribution of sizes at flowering, elasticity of the kernal matrix, comparison with matrix life history traits) 
    
    - Optimal flowering strategy : determine if the flowering strategy is an ESS, calculate the optimal growth rate and the corresponding optimal size at flowering according to variations of flowering intercept.
    
    - Climate : correlations between growth rate, size at flowering and life history traits with mean temperature and number of wet days
    
    - Validation Flowering data : predictions of the IPM comparing with the number of observed flowering plants each year
    
    - Perturbation analysis : elasticity analysis on the vital rates
=======
    - Kernal analysis : global results of the IPM (growth rate, size distribution, ...) 
    - Climate : links between life history traits and optimal flowering size with climate
    - Optimal flowering strategy : determine if the flowering strategy is an ESS
    - Validation Flowering data : predictions of the model comparing with the number of observed flowering plants each year
    - Perturbation analysis : elasticity analysis on the vital rates
>>>>>>> 86363ae7a10db5f790a0cb1506ea8e2d92a9eca9
