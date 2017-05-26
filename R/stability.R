stability <- function(x, plot = TRUE, ...) {
  
  ## for results of 'lmInfl'
  if (class(x) == "influencers") {
    stab <- 1 - (length(x$sel)/nrow(model.frame(x$origModel)))
    cat("Stability Index:", stab, "\n")
  }
  
  ## for results of 'lmMult'
  else if (!is.na(class(x$sample)[2])) {
    stab <- 100 - x$stat
    cat("Stability Index:\n")
    print(stab)
  }
  
  ## for results of 'lmThresh'
  else if (class(x) == "threshsearch") {
    ## get original data
    MODEL <- x$model
    DATA <- model.frame(MODEL); X <- DATA[, 2]; Y <- DATA[, 1]
    
    ## get prediction interval
    PRED <- predict(MODEL, interval = "prediction", level = 1-x$alpha)
    
    ## preallocat
    sVec <- pDiffVec <- as.numeric(rep(NA, nrow(DATA)))
    
    for (i in 1:nrow(DATA)) {
      ## for all fitted values + prediction intervals, calculate
      ## the corresponding scaled/shifted t-distribution by
      ## calculating the scaling parameters s from qt.scaled.
      ## Q = qt(P, df) * s + mu => s = (Q - mu)/qt(P, df)
      ## For example, for x[1]: Qupper = PRED[1, 3], mu = PRED[1 ,1] = y_hat
      s <- (PRED[i, 3] - PRED[i, 1])/qt(1 - (x$alpha/2), df.residual(MODEL))
      sVec[i] <- s
      
      ## calculate integral of threshold border to alpha border
      THRESH <- x$thresh[[i]]
      if (is.null(THRESH) | length(THRESH) == 1) {
        pDiffVec[i] <- NA
        next
      }
      if (THRESH[1] < PRED[i, 2] & THRESH[2] < PRED[i, 3])
        pDiff <- integrate(dt.scaled, lower = THRESH[2], upper = PRED[i, 3], 
                           df = df.residual(MODEL), mu = PRED[i, 1], s = s)$value
      if (THRESH[1] > PRED[i, 2] & THRESH[2] > PRED[i, 3]) 
        pDiff <- integrate(dt.scaled, lower = PRED[i, 2], upper = THRESH[1], 
                           df = df.residual(MODEL), mu = PRED[i, 1], s = s)$value
      if (THRESH[1] < PRED[i, 2] & THRESH[2] > PRED[i, 3]) 
        pDiff <- 0
      pDiffVec[i] <- 1 - pDiff
    }
    
    stab <- pDiffVec
    cat("Stability Index:\n")
    print(stab)
    if (plot) barplot(stab, col = "dodgerblue3", ylim = c(0, 1), ylab = "P(reversal)",
            xlab = "Predictor #", names.arg = as.character(1:nrow(DATA)), ...)
  }
  
    invisible(stab)
}

