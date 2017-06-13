## counter function for 'inflSim'
counter <- function (i) 
{
  if (i%%10 == 0) 
    cat(i)
  else cat(".")
  if (i%%50 == 0) 
    cat("\n")
  flush.console()
}

## define optimizing function
optFct <- function(par, x, y, iter, alpha) {
  y[iter] <- par
  tempMod <- lm(y ~ x)
  PVAL <- corr.test(x, y)$p 
  return(abs(alpha - PVAL))
}

## define interval function
is.between <- function(x, int) {
  int1 <- int[-length(int)]
  int2 <- int[-1]
  (x > int1 & x <= int2) | (x <= int1 & x > int2)
}

## define correlation function for fast p-value calculation
## taken from the 'corr.test' function of the 'psych' package
corr.test <- function (x, y) 
{
  r <- cor(x, y, use = "complete.obs", method = "pearson")
  n <- t(!is.na(x)) %*% (!is.na(y)); n <- min(n)
  t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
  p <- 2 * (1 - pt(abs(t), (n - 2)))
  se <- sqrt((1 - r * r)/(n - 2))
  return(list(p = p, t = t, se = se, r = r))
}

dt.scaled <- function(x, df, mu, s) 1/s * dt((x - mu)/s, df)
pt.scaled <- function(x, df, mu, s) pt((x - mu)/s, df)
qt.scaled <- function(prob, df, mu, s) qt(prob, df) * s + mu
rt.scaled <- function(n, df, mu, s) rt(n,df) * s + mu
