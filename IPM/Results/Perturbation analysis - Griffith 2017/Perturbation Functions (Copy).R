#
# -The functions below are for Approaches 1-4. Perturbations for Approach 5 are described
#  in 'Approach 5 templates.R' and are implemented within 'MAIN.R'.
#
# -Most of the functions below require an input argument that is a list containing
#  all of objects used to construct the overall IPM kernel. For the IPM described in
#  the text, some/all of the folling objects are required within the list, depending
#  on the perturbation approach:
#
# n         -Number of size classes in the discretized kernel
# z         -Midpoint size value for each size class
# z.extend  -Extended size classes for use in eviction correction
# h         -size class width
# K         -IPM projection kernel (n by n array)
# P         -Survival-growth additive sub-kernel (n by n array)
# F         -Reproduction additive sub-kernel (n by n array)
# S         -Size-dependent survival (length n numeric)
# G         -Size transition probabilities for surviving individuals (n by n array)
# G.mu      -Size-dependent future size mean (length n numeric)
# G.sig     -Size-independent future size SD (length 1 numeric)
# R         -Size-dependent fecundity (length n numeric)
# OS        -Size-independent offspring survival (length 1 numeric)
# OG        -Size-independent offspring size probabilities (length n numeric)
# OG.mu     -Size-independent offspring size mean (length 1 numeric)
# OG.sig    -Size-independent offspring size SD (length 1 numeric)
#
#
# -The list objects contained within 'Kernels.RData' adhere to this format and may be
#  passed directly into the functions below.


#Function to extract (first-level) objects from a list
extract <- function(ListObj){
   
    vars <- names(ListObj)
    
    for (i in 1:length(vars)){
        
        assign(vars[i], ListObj[[i]], pos=parent.frame(1))
        
    }
}


#Function to calculate lambda (dominant eigenvalue) from a matrix / discretized kernel
get.Lam <- function(K) Re(eigen(K)$values[1])


###                                                       ###
### - - - - - Approach 1 Perturbation Functions - - - - - ###
###                                                       ###


# - - Projection kernel (Eq. 4 & 5)
do.Pert1.K <- function(K){
    # -Requires K as a square matrix
    # -Note that this is not the most robust script for the analytical determination
    #  of sensitivity values. It's possible that it will fail depending on the values
    #  in the K matrix. The 'sensitivity' function in the 'popbio' package is more robust.
    
    
    #Stable stage distribution - dominant left eigenvector
    W <- eigen(K)$vectors
    w <- Re(W[,1]/sum(W[,1]))  
    w <- as.matrix(w)
    
    #Reproductive value - dominant right eigenvector
    V <- Conj(solve(W))
    v <- Re(V[1,]/V[1,1])
    v <- as.matrix(v)
    
    #Sensivity matrix
    SV <- (v %*% t(w)) / as.numeric((t(v) %*% w))
    EV <- SV * (K / get.Lam(K))
    return(list(SV=SV, EV=EV))
}


# - - Projection kernel - Numerical
do.Pert1.K.num <- function(K,p){
    #Requires K as a square matrix
    
    Lambda <- get.Lam(K)
    n <- dim(K)[1]
    
    #Preallocating
    SV <- K*0
    
    for (j in 1:n){
        
        for (i in 1:n){
            
            #Perturbation
            K.p <- K
            K.p[i,j] <- K[i,j] + p
            
            Lambda.p <- get.Lam(K.p) 
            
            SV[i,j] <- (Lambda.p - Lambda) / p
            
        }
    }
    
    EV <- SV * (K / Lambda)
    return(list(SV=SV, EV=EV))
    
}


###                                                       ###
### - - - - - Approach 2 Perturbation Functions - - - - - ###
###                                                       ###


# - - Any additive Sub-Kerrnel (Eq. 6)
do.Pert2.Sub <- function(K,Sub){
    #Requires K and Sub (e.g. P or F) as square matrices
    
    SV <- do.Pert1.K(K)$SV
    EV <- SV * (Sub / get.Lam(K))
    return(list(SV=SV, EV=EV))
    
}


# - - Survival/Growth Sub-Kernel - Numerical
do.Pert2.P.num <- function(Kernel,p){
    #Kernel must be a list object as described above
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    
    #Preallocating
    SV <- P*0
    
    for (j in 1:n){
        
        for (i in 1:n){
            
            #Perturbation
            P.p <- P
            P.p[i,j] <- P.p[i,j] + p
            
            K.p <- P.p + F
            
            Lambda.p <- get.Lam(K.p)
            
            SV[i,j] <- (Lambda.p - Lambda) / p
            
        }
    }
    
    EV <- SV * (P / Lambda)
    return(list(SV=SV, EV=EV))
    
}


# - - Reproduction Sub-Kernel - Numerical
do.Pert2.F.num <- function(Kernel,p){
    #Kernel must be a list object as described above
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    
    #Preallocating
    SV <- F*0
    
    for (j in 1:n){
        
        for (i in 1:n){
            
            #Perturbation
            F.p <- F
            F.p[i,j] <- F.p[i,j] + p
            
            K.p <- F.p + P
            
            Lambda.p <- get.Lam(K.p)
            
            SV[i,j] <- (Lambda.p - Lambda) / p
            
        }
    }
    
    EV <- SV * (F / Lambda)
    return(list(SV=SV, EV=EV))
    
}


###                                                       ###
### - - - - - Approach 3 Perturbation Functions - - - - - ###
###                                                       ###
library(popbio)
library(IPMpack)

# - - Survival (Eq. 9 & A3)
do.Pert3.S <- function(K){
    #K must be a list object as described above
    Lambda <- lambda(K)
    SV.K <- IPMpack::sens(K)
    
    #Preallocating
    SV <- S*0
    #column loop
    for (j in 1:n){
        SUM <- 0
        #row loop
        for (i in 1:n){
            SUM <- SUM + SV.K[i,j] * G[i,j] * (1-Fl[j])
            }
        SV[j] <- SUM
    }
    EV <- SV * (S / Lambda)
    return(list(SV=SV, EV=EV))
}


# - - Size transition probabilities (Eq. 12 & A4) - Proportional compensation - Integrated Sensitivities
do.Pert3.G <- function(K){
  Lambda <- lambda(K)
  SV.K <- IPMpack::sens(K)
    #Preallocating
    SV <- G
    for (j in 1:n){
        for (i in 1:n){
            SUM <- 0
            for (m in 1:n){
                if (m != i){
                    PD <- -(G[m,j] / sum(G[-i,j]))
                    SUM <- SUM + PD * SV.K[m,j] * S[j] * (1-Fl[j])
                }
            }
            SV[i,j] <- SV.K[i,j] * S[j] * (1-Fl[j]) + SUM
        }
    }
    EV <- SV * (G / Lambda )
    
    return(list(SV=SV, EV=EV))
}

# - - Flowering
do.Pert3.Fl <- function(K) {
  Lambda <- lambda(K)
  SV.K <- IPMpack::sens(K)
  
  #Preallocating
  SV <- Fl * 0
  
  for (j in 1:n) {
    SUM <- 0
    for (i in 1:n) {
      SUM <- SUM + SV.K[i, j] * R[j] * Er * OG[i]
    }
    SV[j] <- SUM
  }
  EV <- SV * (Fl / Lambda)
  
  return(list(SV = SV, EV = EV))
}


# - - Fecundity (Eq. 10 & A6)
do.Pert3.R <- function(K) {
  Lambda <- lambda(K)
  SV.K <- IPMpack::sens(K)
  
  #Preallocating
  SV <- R * 0
  
  for (j in 1:n) {
    SUM <- 0
    for (i in 1:n) {
      SUM <- SUM + SV.K[i, j] * Fl[j] * Er * OG[i]
    }
    SV[j] <- SUM
  }
  EV <- SV * (R / Lambda)
  
  return(list(SV = SV, EV = EV))
}


# - - Offspring size (Eq. A8) - Proportional compensation - Integrated Sensitivities
do.Pert3.OG <- function(K) {
  Lambda <- lambda(K)
  SV.K <- IPMpack::sens(K)
  
  #Preallocating
  SV <- array(0, c(n, n))
  EV <- SV
  
  for (j in 1:n) {
    for (i in 1:n) {
      SUM <- 0
      
      for (m in 1:n) {
        if (m != i) {
          PD <- -(OG[m] / sum(OG[-i]))
          SUM <- SUM + PD * SV.K[m, j] * Fl[j] * R[j] * Er
        }
      }
      SV[i, j] <- SV.K[i, j] * Fl[j] * R[j] * Er + SUM
      EV[i, j] <- SV[i, j] * (OG[i] / Lambda)
    }
  }
  return(list(SV = SV, EV = EV))
}


# - - Establishment rate
do.Pert3.Er <- function(K) {
  Lambda <- lambda(K)
  SV.K <- IPMpack::sens(K)
  
  #Preallocating
  SV <- Er * 0
  
  for (j in 1:n) {
    SUM <- 0
    for (i in 1:n) {
      SUM <- SUM + SV.K[i, j] * R[j] * Fl[j] * OG[i]
    }
    SV[j] <- SUM
  }
  EV <- SV * (Er / Lambda)
  
  return(list(SV = SV, EV = EV))
}

# - - Size transition probabilities (Eq. A9 & A10) - Proportional compensation - Numerical (SLOW)
do.Pert3.G.num <- function(Kernel,p){
    #Kernel must be a list object as described above
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    #Preallocating
    SV <- G*0
    P.p <- P*0
    G.sums <- colSums(G)
    
    for (j in 1:n){
        
        for (i in 1:n){
            
            G.p <- G
            
            #Rescales non-target elements
            B <- G[,j] * ( (G.sums[j] - G[i,j] - p) / (G.sums[j] - G[i,j]) )
            
            G.p[,j] <- B
            
            #Perturbation
            G.p[i,j] <- G[i,j] + p
            
            
            #make P kernel
            for(k in 1:n) P.p[,k] <- G.p[,k] * S[k]
            
            K.p <- P.p + F
            
            Lambda.p <- get.Lam(K.p) 
            
            SV[i,j] <- (Lambda.p - Lambda) / p
            
        }
    }
    
    EV <- SV * (G / Lambda)
    return(list(SV=SV, EV=EV))
    
}


###                                                       ###
### - - - - - Approach 4 Perturbation Functions - - - - - ###
###                                                       ###


# - - Survival - Numerical
do.Pert4.S <- function(Kernel,p){
    #Kernel must be a list object as described above
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    #Preallocating
    SV <- S*0
    P.p <- P*0
    
    for (j in 1:n){
        
        #Perturbation
        S.p <- S
        S.p[j] <- S.p[j] + p
        
        #make P kernel
        for(k in 1:n) P.p[,k] <- G[,k] * S.p[k]
        
        K.p <- P.p + F
        
        Lambda.p <- get.Lam(K.p)
        
        SV[j] <- (Lambda.p - Lambda) / p
        
    }
    
    EV <- SV * (S / Lambda)
    return(list(SV=SV, EV=EV))
    
}


# - - Mean future size - Numerical
do.Pert4.G.mu <- function(Kernel,p){
    #Kernel must be a list object as described above    
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    i.in <- which(is.element(z.extend, z)) #index of future sizes inside kernel
    i.above <- which(z.extend > max(z)) #index of too large
    i.below <- which(z.extend < min(z)) #index of too small
    
    #Preallocating
    SV <- G.mu*0
    P.p <- P*0
    
    for (j in 1:n){
        
        #Perturbation
        G.mu.p <- G.mu
        G.mu.p[j] <- G.mu.p[j] + p
        
        #G matrix
        #only need to correct for eviction in perturbed column j
        G.p <- G #start with unperturbed G matrix
        G.extend <- dnorm(z.extend, G.mu.p[j], G.sig) * h #calculated perturbed extended distribution (sum = 1)
        G.p[,j] <- G.extend[i.in] #insert perturbed column
        G.p[n,j] <- G.p[n,j] + sum(G.extend[i.above]) #insert large eviction prob 
        G.p[1,j] <- G.p[1,j] + sum(G.extend[i.below]) #insert small eviction prob
        
        #make P kernel
        for(k in 1:n) P.p[,k] <- G.p[,k] * S[k]
        
        K.p <- P.p + F
        
        Lambda.p <- get.Lam(K.p) 
        
        SV[j] <- (Lambda.p - Lambda) / p
        
    }
    
    EV <- SV * (G.mu / Lambda)
    return(list(SV=SV, EV=EV))
    
}


# - - Future size SD - Numerical
do.Pert4.G.sig <- function(Kernel,p){
    #Kernel must be a list object as described above    
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    i.in <- which(is.element(z.extend, z)) #index of future sizes inside kernel
    i.above <- which(z.extend > max(z)) #index of too large
    i.below <- which(z.extend < min(z)) #index of too small
    
    #Preallocating
    P.p <- P*0
    G.p <- G*0
    
    #Perturbation
    G.sig.p <- G.sig + p
    
    #G matrix
    for (j in 1:n){
        
        G.extend <- dnorm(z.extend, G.mu[j], G.sig.p) * h #calculated perturbed extended distribution (sum = 1)
        G.p[,j] <- G.extend[i.in] #insert into G matrix with eviction
        G.p[n,j] <- G.p[n,j] + sum(G.extend[i.above]) #insert large eviction prob 
        G.p[1,j] <- G.p[1,j] + sum(G.extend[i.below]) #insert small eviction prob
        
    }
    
    #make P kernel
    for(j in 1:n) P.p[,j] <- G.p[,j] * S[j]
    
    K.p <- P.p + F
    
    Lambda.p <- get.Lam(K.p) 
    
    SV <- (Lambda.p - Lambda) / p
    
    EV <- SV * (G.sig / Lambda)
    
    return(list(SV=SV, EV=EV))
    
}


# - - Fecundity - Numerical
do.Pert4.R <- function(Kernel,p){
    #Kernel must be a list object as described above
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    #Preallocating
    SV <- R*0
    F.p <- P*0
    
    for (j in 1:n){
        
        #Perturbation
        R.p <- R
        R.p[j] <- R.p[j] + p
        
        for(k in 1:n) F.p[,k] <- R.p[k] * OS * OG
        
        K.p <- P + F.p
        
        Lambda.p <- get.Lam(K.p)
        
        SV[j] <- (Lambda.p - Lambda) / p
        
    }
    
    EV <- SV * (R / Lambda)
    return(list(SV=SV, EV=EV))
    
}


# - - Offspring survival - Numerical
do.Pert4.OS <- function(Kernel,p){
    #Kernel must be a list object as described above
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    #Preallocating
    F.p <- P*0
    
    #Perturbation
    OS.p <- OS + p
    
    #make F kernel
    for(k in 1:n) F.p[,k] <- R[k] * OS.p * OG
    
    K.p <- P + F.p
    
    Lambda.p <- get.Lam(K.p)
    
    SV <- (Lambda.p - Lambda) / p
    
    EV <- SV * (OS / Lambda)
    
    return(list(SV=SV, EV=EV))
    
}


# - - Mean offspring size - Numerical
do.Pert4.OG.mu <- function(Kernel,p){
    #Kernel must be a list object as described above
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    #Preallocating
    SV <- OG.mu*0
    F.p <- F*0
    
    #Perturbation
    OG.mu.p <- OG.mu + p 
    
    OG.p.orig <- dnorm(z, OG.mu.p, OG.sig)*h
    OG.p <- OG.p.orig / sum(OG.p.orig)
    
    #make F kernel
    for(k in 1:n) F.p[,k] <- R[k] * OS * OG.p
    
    K.p <- P + F.p
    
    Lambda.p <- get.Lam(K.p) 
    
    SV <- (Lambda.p - Lambda) / p
    
    EV <- SV * (OG.mu / Lambda)
    return(list(SV=SV, EV=EV))
    
}


# - - Offspring size SD - Numerical
do.Pert4.OG.sig <- function(Kernel,p){
    #Kernel must be a list object as described above
    
    extract(Kernel)
    Lambda <- get.Lam(K)
    
    #Preallocating
    SV <- OG.sig*0
    F.p <- F*0
    
    #Perturbation
    OG.sig.p <- OG.sig + p 
    
    OG.p.orig <- dnorm(z, OG.mu, OG.sig.p)*h
    OG.p <- OG.p.orig / sum(OG.p.orig)
    
    #make F kernel
    for(k in 1:n) F.p[,k] <- R[k] * OS * OG.p
    
    K.p <- P + F.p
    
    Lambda.p <- get.Lam(K.p) 
    
    SV <- (Lambda.p - Lambda) / p
    
    EV <- SV * (OG.sig / Lambda)
    
    return(list(SV=SV, EV=EV))
    
}

