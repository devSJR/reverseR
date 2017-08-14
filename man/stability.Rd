\name{stability}
\alias{stability}
\encoding{latin1}

\title{Calculates stability values for results of 'lmInfl' and 'lmMult'}

\description{
This function calculates stability values for either LOO (\code{\link{lmInfl}}) or LMO  (\code{\link{lmMult}}) models.
}

\usage{
stability(x, plot = TRUE, ...) 
}

\arguments{
\item{x}{a result of either \code{\link{lmInfl}} or \code{\link{lmMult}}.}
\item{plot}{logical. If \code{TRUE}, the stability is plotted for each response value.}
\item{...}{other parameters to \code{\link{barplot}}.}
}

\details{
For results of \code{\link{lmInfl}}:\cr
A [0, 1]-bounded stability measure \eqn{S = 1-\frac{n}{N}}, with \eqn{n} = number of influencers (significance reversers) and \eqn{N} = total number of response values.\cr\cr
For results of \code{\link{lmMult}}:\cr
For each 1...\code{max}, the percentage of all resamples that did *NOT* result in significance reversal.\cr\cr
}

\value{
The stability value.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## See examples in 'lmInfl' and 'lmThresh'.
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

