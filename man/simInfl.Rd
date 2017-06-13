\name{simInfl}
\alias{simInfl}
\encoding{latin1}

\title{Simulates significance reversals and calculates their influence parameters}

\description{
This function simulates linear regressions and stores the parameters of all simulations that resulted in a LOO significance reversal.
}

\usage{
simInfl(x = 1:10, slope = 0.02, error = 0.05, nrev = 1000, ...) 
}

\arguments{
\item{x}{the \eqn{x} values to be supplied to \code{\link{lmExact}}.}
\item{slope}{the slope \eqn{\beta_1} to be supplied to \code{\link{lmExact}}.}
\item{error}{the \eqn{\varepsilon} value to be supplied to \code{\link{lmExact}}.}
\item{nrev}{the number of desired significance reversals.}
\item{...}{other parameters to \code{\link{lmExact}}.}
}

\details{
Loops over an undefined number of EXACT regressions (\code{\link{lmExact}}) with incrementing random seeds, tests for significance reversal, stores the parameters (original \emph{p}-values, the output of \code{\link{lmInfl}}) when this is the case and terminates when \code{nrev} reversals are counted.
}

\value{
A matrix with the resulting parameters for all \code{nrev} reversals.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
\dontrun{
## Example with slight slope and 1000 reversals.
res <- inflSim(x = 1:10, slope = 0.02, error = 0.05, nrev = 1000)

## Plot Cook's D versus delta-P values
## and insert common cut-offs.
plot(res$cook.d, res$dP, pch = 16, cex = 0.5)
abline(v = 4/10, col = "darkred", lwd = 2)
abline(v = 3 * mean(res$cook.d, na.rm = TRUE), col = "darkred", lwd = 2)
}
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

