library(SplinesUtils)

#Pour un exemple :
# RegModel <- fitme(Survie ~ 1 + bs(Size0Mars, df = 4, degree = 2) +(bs(Age, degree = 3, knots = 6.5)) + (Age|year) + (1|Pop),
#                   family=binomial,
#                   data=survdata2,
#                   method="PQL/L")
# SplineTerm <- "bs(Size0Mars, df = 4, degree = 2)"

# La fonction d'origine est RegSplineAsPiecePoly

SplinesAsPolySpaMM <- function (RegModel, SplineTerm, shift = TRUE) {
  
  #Extrait les données du spline (position dans les effets fixes, degrés des polynomes, knots)
  RegSpline <- SplinesUtils:::ExtractSplineTerm(terms(RegModel), SplineTerm)
  pos <- RegSpline$pos
  
  #Permet de sortir les coefficients de chaque bout de splines
  
  ind0 <- attr(RegModel$X.pv, "namesOri") #Noms des coefficients
  ind <- integer(0)
  pattern <- paste0("^", gsub("\\(", "\\\\(", gsub("\\)", "\\\\)", SplineTerm)), "\\d+$") 
  #C'est pas très propre mais pas le choix car les noms des bouts de splines finissent avec un indice 1,2,...
  ind <- grep(pattern, ind0) #Détecte si le nom du coefficient correspond au spline voulu
  RegSplineCoef <- RegModel$fixef[ind] #Sort l'estimation des coeffs
  
  RegSplineCoef <- unname(RegSplineCoef)
  na <- is.na(RegSplineCoef)
  if (any(na)) {
    warning("NA coefficients found for SplineTerm; Replacing NA by 0")
    RegSplineCoef[na] <- 0
  }
  
  #Recalcule le polynome avec la méthode d'origine, je n'y ai pas touché
  PiecePolyCoef <- SplinesUtils:::PiecePolyRepara(RegSpline, RegSplineCoef, shift)
  structure(list(
    PiecePoly = list(coef = PiecePolyCoef, shift = shift),
    knots = RegSpline$knots),
    class = c("PiecePoly", "RegSpline"))
}