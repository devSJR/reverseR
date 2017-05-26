lmInfl <- function(model, alpha = 0.05, verbose = TRUE, ...) {
  
  if (class(model) != "lm") stop("Not an 'lm' model!")
  
  ## original parameters
  DATA <- eval(model$model)
  X <- DATA[, 2]
  Y <- DATA[, 1]
  PVAL.old <- summary(model)$coefficients[2, 4]
  SLOPE.old <- coef(model)[2]
  INTER.old <- coef(model)[1]
  SE.old <- summary(model)$coefficients[2, 2]
  wts <- model$weights
  
  ## print standard parameters
  if (verbose) {
    if (PVAL.old <= alpha) cat("\nModel is significant at p = ", PVAL.old, ".\n", sep = "")
    else cat("\nModel is insignificant at p = ", PVAL.old, ".\n", sep = "")
  }
  
  ## define change direction
  if (PVAL.old > alpha) change <- "down" else change <- "up"
  
  ## classical influence measures
  INFL <- influence.measures(model)[[1]]
  
  ## create Leave-One-Out models
  modList <- vector("list", length = length(X))
    
  for (i in 1:length(X)) {
    modList[[i]] <- lm(Y[-i] ~ X[-i], weights = wts, ...)
  }
    
  ## get coefficients, p-values, slopes and s.e.(slope)
  Coef <- sapply(modList, function(x) coef(x))
  PVAL.loo <- sapply(modList, function(x)  summary(x)$coefficients[2, 4])
  deltaP <- abs(PVAL.old - PVAL.loo)
  Slope <- sapply(modList, function(x)  summary(x)$coefficients[2, 1])
  Inter <- sapply(modList, function(x)  summary(x)$coefficients[1, 1])
  SE <- sapply(modList, function(x)  summary(x)$coefficients[2, 2])
  deltaSlope <- abs(Slope - SLOPE.old)
  deltaInter <- abs(Inter - INTER.old)
  deltaSE <- abs(SE - SE.old)
  
  ## return extended influence matrix 
  inflMat <- cbind(slope = Coef[2, ], Inter = Coef[1, ], 
                   dSlope = deltaSlope, dInter = deltaInter, INFL, 
                   leverage = INFL[, 6]/(1 - INFL[, 6]), looP = PVAL.loo, 
                   dP = deltaP, SE = SE, dSE = deltaSE)
  rownames(inflMat) <- 1:length(X)
  
  ## select all influencers 
  if (change == "up") SEL <- which(PVAL.loo > alpha)
  if (change == "down") SEL <- which(PVAL.loo <= alpha)
  
  if (length(SEL) > 0) {
    if (verbose) cat("Found the following significance reversers, in order of strength:\n")
    selP <- PVAL.loo[SEL]
    
    ## order influencers by decreasing strength
    ordP <- order(selP, decreasing = TRUE)
    ordSEL <- SEL[ordP]
    
    ## print ordered influencers
    if (verbose) {
      for (i in ordSEL) cat("  Point #", i, " => p = ", inflMat[i, "looP"], ".\n", sep = "")
    }  
  }  
  else {
    selMat <- NULL
    if (verbose) cat("No significance reversers found!\n")
    ordSEL <- NULL
  }
  
  OUT <- list(origModel = model, finalModels = modList, infl = inflMat, 
              sel = ordSEL, alpha = alpha, origP = PVAL.old)
  class(OUT) <- "influencers"
  return(OUT)
}

lmPlot <- function(infl, ...) {
  
  ## check for correct class
  if (class(infl) != "influencers") stop("object is not a result of 'lmInfl' !")
  
  ## plot full data
  DATA <- eval(infl$origModel$model)
  X <- DATA[, 2]
  Y <- DATA[, 1]
  plot(X, Y, pch = 16, xlab = colnames(DATA)[2], ylab = colnames(DATA)[1], ...)
  grid()
  
  ## define colors
  COL <- rainbow(length(infl$sel))
  CEX <- seq(from = 2, to = 1, length.out = length(infl$sel))
  
  ## add color/size for influencers
  points(X[infl$sel], Y[infl$sel], pch = 16, cex = 1.2, col = "red3", ...)
  points(X[infl$sel], Y[infl$sel], pch = 16, cex = CEX, col = COL, ...)
  
  ## original model regression line
  abline(infl$origModel, lwd = 2, col = "black", ...)
  
  ## LOO-models regression line
  if (length(infl$finalModels) > 0) {
    for (i in 1:length(infl$finalModels))
    abline(infl$finalModels[[i]], lwd = 2, col = COL[i], ...)
  }
}

pvalPlot <- function(infl, ...) {
  
  ## check for correct class
  if (class(infl) != "influencers") stop("object is not a result of 'lmInfl' !")
  
  ## extract leave-one-out p-values
  looP <- infl$infl[, "looP"]
  
  ## plot p-values
  par(mar = c(5, 5, 4, 2))
  plot(infl$origModel$model[, 2],  looP, pch = 16, type = "o", las = 0, 
       xlab = colnames(infl$origModel$model)[2], ylab = "P-value", 
       log = "y", ylim = c(0.75 * min(looP, na.rm = TRUE), 1.25 * max(looP, na.rm = TRUE)), 
       las = 1, ...)
  grid()
  
  ## add full model p-value and alpha border
  abline(h = infl$origP, col = "blue", lwd = 2, ...)
  abline(h = infl$alpha, col = "red3", lwd = 2, ...)
  legend("top", c("Original model p-value", "Selected alpha"), 
         fill = c("blue", "red3"), horiz = TRUE, inset = -0.25, xpd = TRUE, ...)
}

inflPlot <- function(infl, ...) {
  
  ## check for correct class
  if (class(infl) != "influencers") stop("object is not a result of 'lmInfl' !")
  
  ## extract influence values
  df_a <- infl$infl[, "dfb.a"]
  df_fit <- infl$infl[, "dffit"]
  cov_r <- infl$infl[, "cov.r"]
  cook_D <- infl$infl[, "cook.d"]
  looP <- infl$infl[, "looP"]
  
  ## define minimum y-axis plot value
  if (min(looP, na.rm = TRUE) <= infl$alpha) MIN <- min(looP, na.rm = TRUE) 
    else MIN <- infl$alpha
  
  ## define sample size
  N <- nrow(infl$origModel$model)
  
  ## plot influence values vs. p-value
  par(mfrow = c(2, 2))
  par(mar = (c(5, 4, 1, 1)))
  plot(df_a, looP, pch = 16, xlab = "dfbeta slope", ylab = "p-value", 
       las = 1, log = "y", ylim = c(MIN, 1), ...)
  points(df_a[infl$sel], looP[infl$sel], pch = 16, col = "red", ...)
  abline(h = infl$alpha, col = "blue", lwd = 2, ...)
  abline(v = 2/sqrt(N), col = "darkgreen", lwd = 2, ...) ## 2/sqrt(n)
  grid()
  
  plot(df_fit, looP, pch = 16, xlab = "dffits", ylab = "p-value", 
       las = 1, log = "y", ylim = c(MIN, 1), ...)
  points(df_fit[infl$sel], looP[infl$sel], pch = 16, col = "red", ...)
  abline(h = infl$alpha, col = "blue", lwd = 2, ...)
  abline(v = 2 * sqrt(2/N), col = "darkgreen", lwd = 2, ...)  ## 2*sqrt(2/n)
  grid()
  
  plot(abs(cov_r - 1), looP, pch = 16, xlab = "covratio", ylab = "p-value", 
       las = 1, log = "y", ylim = c(MIN, 1), ...)
  points(abs(cov_r - 1)[infl$sel], looP[infl$sel], pch = 16, col = "red", ...)
  abline(h = infl$alpha, col = "blue", lwd = 2, ...)
  abline(v = 3 * (2/N), col = "darkgreen", lwd = 2, ...)  ## 3*sqrt(2/n)
  grid()
  
  plot(cook_D, looP, pch = 16, xlab = "CookD", ylab = "p-value", 
       las = 1, log = "y", ylim = c(MIN, 1), ...)
  points(cook_D[infl$sel], looP[infl$sel], pch = 16, col = "red", ...)
  abline(h = infl$alpha, col = "blue", lwd = 2, ...)
  abline(v = 4/N, col = "darkgreen", lwd = 2, ...)
  abline(v = 3 * mean(cook_D, na.rm = TRUE), col = "darkgreen", lwd = 2, ...)
  grid()
}

slsePlot <- function(infl, ...) {
  
  ## check for correct class
  if (class(infl) != "influencers") stop("object is not a result of 'lmInfl' !")
  
  ## slope and SE of original model
  origSlope <- summary(infl$origModel)$coefficients[2, 1]
  origSE <- summary(infl$origModel)$coefficients[2, 2]
  
  ## extract LOO-slopes and LOO-SEs
  Slope <- infl$infl[, "slope"]
  SE <- infl$infl[, "SE"]
  
  ## plot both 
  plot(Slope, SE, type = "p", pch = 16, xlab = "Slope",
       ylab = "SE", las = 1)
  points(Slope[infl$sel], SE[infl$sel], type = "p", 
         cex = 1.2, pch = 16, col = "red3")
  grid()
  
  ## insert t-border for p = alpha
  Seq <- seq(0.5 * min(Slope), 2 * max(Slope), length.out = 2)
  Quant <- qt(1 - 0.5 * infl$alpha, df.residual(infl$origModel))
  lines(Seq, Seq/Quant, type = "l", col = "black", lwd = 2)
  
  ## insert original model Slope and SE
  abline(v = origSlope, lwd = 2, col = "blue", lty = 2)
  abline(h = origSE, lwd = 2, col = "darkgreen", lty = 2)
  legend("top", c("Original slope", "Original s.e.", "t-border"), 
         col = c("blue", "darkgreen", "black"), lty = c(3, 3, 1), lwd = 3, 
         xpd = TRUE, inset = -0.25, horiz = TRUE, ...)
}



