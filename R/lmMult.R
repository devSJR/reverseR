lmMult <- function(model, max = 5, n = 10000, alpha = 0.05, verbose = TRUE) {
  
  ## extract model data
  DATA <- model.frame(model)
  NR <- nrow(DATA)
  
  ## get original p-value
  origP <- summary(model)$coefficients[2, 4]
  
  ## preallocate 
  sampleMat <- matrix(NA, nrow = max * n, ncol = NR)
  pList <- gList <- vector("list", length = max)
  nRev <- rep(NA, max)
  iter <- 1
  
  ## create combination list
  for (i in 1:max) {
    if (verbose) cat("Sampling ", i, " out of ", NR, " => ", sep = "") 
    
    ## preallocate p-value vector and group vector
    pVec <- rep(NA, n)
    gVec <- rep(i, n)
    
    ## for all n, sample n of N, delete samples and calculate p
    for (j in 1:n) {
     SEL <- sample(1:NR, i, replace = FALSE)
     sampleMat[iter, SEL] <- 1
     X <- DATA[-SEL, 2]
     Y <- DATA[-SEL, 1]
     pVec[j] <- corr.test(X, Y)$p
     iter <- iter + 1
    }
    
    ### calculate number of significance reversers
    if (origP <= alpha) nRev[i] <- round(sum(pVec > alpha)/length(pVec) * 100, 2)
    else nRev[i] <- round(sum(pVec <= alpha)/length(pVec) * 100, 2)
    if (verbose) cat(nRev[i], "% significance reversers.\n", sep = "")
    
    ## store p-values and groups in list
    pList[[i]] <- pVec
    gList[[i]] <- gVec
  }
  
  ## concatenate sample matrix and p-values
  pAllVec <- unlist(pList)
  gAllVec <- unlist(gList)
  sampleMat <- cbind(sampleMat, pval = pAllVec, group = gAllVec)
  sampleMat <- as.data.frame(sampleMat)
  
  ## remove redundant samples
  sampleMat <- unique(sampleMat)
  if (verbose) cat(max * n, "samplings gave", nrow(sampleMat), "unique samples.\n")
  
  names(nRev) <- as.character(1:max)
  class(sampleMat) <- c("data.frame", "multinfl")
  OUT <- list(sample = sampleMat, stat = nRev)
  attr(OUT, "stats") <- c(origP, alpha)
  return(OUT)
}

multPlot <- function(mult, ...) {
  x <- mult
  if (class(x$sample)[2] != "multinfl") stop("object is not a result of 'searchInfl' !")
  
  ## get original model p-value and alpha border
  origP <- attr(x, "stats")[1]; alpha <- attr(x, "stats")[2]
  X <- x$sample[, "group"]; Y <- x$sample[, "pval"]
  
  ## plot in 'stripchart' style
  X <- X + rnorm(length(X), 0, 0.05)
  plot(X, Y, pch = 16, cex = 0.2, col = ifelse(Y <= alpha, "red3", "black"),
       las = 1, xlab = "left-out samples", ylab = "p-value", cex.lab = 1.5, ...)
}

