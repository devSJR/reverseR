\name{lmInfl}
\alias{lmInfl}
\encoding{latin1}

\title{Checks and analyzes leave-one-out (LOO) p-values in linear regression}

\description{
This function calculates leave-one-out (LOO) \emph{p}-values for all datapoints and identifies those resulting in "significance reversal", i.e. in traversing the user-defined \eqn{\alpha}-level.
}

\usage{
lmInfl(model, alpha = 0.05, verbose = TRUE, ...) 
}

\arguments{
\item{model}{the linear model of class \code{\link{lm}}.}
\item{alpha}{the \eqn{\alpha}-level to use as the threshold border.}
\item{verbose}{logical. If \code{TRUE}, results are dosplayed on the console.}
\item{...}{other arguments to \code{\link{lm}}.}
}

\details{
The algorithm\cr
1) calculates the \emph{p}-value of the full model (all points),\cr
2) calculates a LOO-\emph{p}-value for each point removed,\cr
3) checks for "significance reversal" in all data points and\cr
4) returns all models as well as classical \code{\link{influence.measures}} with LOO-\emph{p}-values, \eqn{\Delta}\emph{p}-values, slopes and standard errors attached.\cr
}

\value{
A list with the following items:\cr
\item{origModel}{the original model with all data points.} 
\item{finalModels}{a list of final models with the influencer(s) removed.}
\item{infl}{a matrix with the original data, the classical \code{\link{influence.measures}}, the leverages, LOO-\emph{p}-values, LOO-slopes/intercepts and their \eqn{\Delta}'s, LOO-standard errors and \eqn{R^2}s.}
\item{sel}{a vector with the influencers' indices.}
\item{alpha}{the selected \eqn{\alpha}-level.}
\item{origP}{the original model's \emph{p}-value.}
\item{stab}{the stability measure, as calculated above.}
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## Example #1 with single influencers and insignificant model (p = 0.115).
## Removal of #18 results in p = 0.0227!
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM1 <- lm(b ~ a)
res1 <- lmInfl(LM1) 
lmPlot(res1)
pvalPlot(res1)
inflPlot(res1)
slsePlot(res1)
stability(res1)

## Example #2 with multiple influencers and significant model (p = 0.0269).
## Removal of #2, #17, #18 or #20 result in crossing p = 0.05!
set.seed(125)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM2 <- lm(b ~ a)
res2 <- lmInfl(LM2) 
lmPlot(res2)
pvalPlot(res2)
inflPlot(res2)
slsePlot(res2)
stability(res2)

## Large Example #3 with top 10 influencers and significant model (p = 6.72E-8).
## Not possible to achieve a crossing of alpha with any point despite of
## strong noise.
set.seed(123)
a <- 1:100
b <- 5 + 0.08 * a + rnorm(100, 0, 5)
LM3 <- lm(b ~ a)
res3 <- lmInfl(LM3) 
lmPlot(res3)
stability(res3)

## Example #4 with replicates and single influencer (p = 0.114).
## Removal of #58 results in p = 0.039.
set.seed(123)
a <- rep(1:20, each = 3)
b <- 5 + 0.08 * a + rnorm(20, 0, 2)
LM4 <- lm(b ~ a)
res4 <- lmInfl(LM4) 
lmPlot(res4)
pvalPlot(res4)
inflPlot(res4)
slsePlot(res4)
stability(res4)

## As Example #1, but with weights.
## Removal of #18 results in p = 0.04747!
set.seed(123)
a <- 1:20
b <- 5 + 0.08 * a + rnorm(20, 0, 1)
LM5 <- lm(b ~ a, weights = 1:20)
res5 <- lmInfl(LM5) 
lmPlot(res5)
stability(res5)
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

