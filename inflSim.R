counter <- function (i) 
{
  if (i%%10 == 0) 
    cat(i)
  else cat(".")
  if (i%%50 == 0) 
    cat("\n")
  flush.console()
}

n <- 10
error <- 0.05 #seq(1, 50, by = 5)/100
slope <- 0.02 #c(0.01, 0.02, 0.05, 0.1)
reps <-  1:100000

GRID <- expand.grid(reps, n, slope, error)
colnames(GRID) <- c("reps", "n", "slope", "error")
MAT <- matrix(NA, nrow = nrow(GRID), ncol = 16)
GRID <- cbind(GRID, MAT)

isRev <- 0

for (i in 1:nrow(GRID)) {
  LME <- lmExact(x = 1:GRID[i, 2], slope = GRID[i, 3], error = GRID[i, 4], 
                 plot = FALSE, seed = i, verbose = FALSE)
  RES <- lmInfl(LME$lm, verbose = FALSE)
  GRID[i, 5] <- RES$origP
  if (!is.null(RES$sel)) {
    counter(isRev)
    GRID[i, 6:20] <- RES$infl[RES$sel[1], ]
    isRev <- isRev + 1
  }
  if (isRev == 10000) break
}
colnames(GRID)[5:20] <- c("origP", colnames(RES$infl))
GRID <- GRID[!is.na(GRID$cook.d), ]

plot(GRID$cook.d, GRID$dP, pch = 16, cex = 0.5)
abline(v = 4/n, col = "darkred", lwd = 2)
abline(v = 3 * mean(GRID$cook.d, na.rm = TRUE), col = "darkred", lwd = 2)
plot(GRID$leverage, GRID$dP, pch = 16, cex = 0.5)
