simInfl <- function(
  x = 1:10, 
  slope = 0.02, 
  error = 0.05, 
  nrev = 1000, 
  ...)
{
  ## create combinations and result grid
  GRID <- expand.grid(1:nrev, slope, error)
  colnames(GRID) <- c("nrev", "slope", "error")
  MAT <- matrix(NA, nrow = nrow(GRID), ncol = 17)
  GRID <- cbind(GRID, MAT)
  
  ## initialize loop counter & reversal counter
  seeds <- 1
  i <- 1
  isRev <- 1
  
  ## loop until isRev = nrev
  while (isRev <= nrev) {
    ## create exact model 
    LME <- lmExact(x = x, slope = GRID[i, 2], error = GRID[i, 3], 
                   plot = FALSE, seed = seeds, verbose = FALSE, ...)
    ## get reversal results
    RES <- lmInfl(LME$lm, verbose = FALSE)
    
    ## if reversal, populate result grid and increase counters
    if (!is.null(RES$sel)) {
      counter(isRev)
      GRID[i, 4] <- seeds
      GRID[i, 5] <- RES$origP
      GRID[i, 6:20] <- RES$infl[RES$sel[1], ]
      isRev <- isRev + 1
      i <- i + 1
    }
    seeds <- seeds + 1
  }
 
  ## create column names and return
  colnames(GRID)[4:20] <- c("seed", "origP", colnames(RES$infl))
  return(GRID)
}
