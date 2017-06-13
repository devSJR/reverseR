\name{Influence plots}
\alias{lmPlot}
\alias{pvalPlot}
\alias{inflPlot}
\alias{slsePlot}
\alias{multPlot}
\encoding{latin1}

\title{Several diagnostic plots for checking p-value influencers}

\description{
Five different plot types that visualize \emph{p}-value influencers.\cr

1. \code{lmPlot}: plots the linear regression, marks the influencer(s) in red and displays trend lines for the full and LOO data set (black and red, respectively).\cr
2. \code{pvalPlot}: plots the \emph{p}-values for each leave-one-out (LOO) data point and displays the values as a full model/LOO model plot, together with the \code{alpha} border as defined in \code{\link{lmInfl}}.\cr
3. \code{inflPlot}: plots \code{\link{dfbeta}} for slope, \code{\link{dffits}}, \code{\link{covratio}}, \code{\link{cooks.distance}} and leverage (\code{\link{hatvalues}}) against the \eqn{\Delta}\emph{p}-value. Herewith, changes in these four parameters can be compared to the effect on the corresponding shift in \emph{p}-value.\cr
4. \code{slsePlot}: plots all LOO-slopes and their standard errors together with the corresponding original model values and a t-value border as calculated by \eqn{\mathit{Q_t}(1 - \frac{\alpha}{2}, n-2)}. LOO of points on the right of this border result in a significant model, and \emph{vice versa}.\cr
5. \code{multPlot}: for output of \code{\link{lmMult}}, plots all \code{n} LMO \emph{p}-values for all 1...\code{max} left-out samples in \code{\link{stripchart}} style.    
}

\usage{
lmPlot(infl, ...) 
pvalPlot(infl, ...) 
inflPlot(infl, ...)
slsePlot(infl, ...)
multPlot(mult, ...)
}

\arguments{
\item{infl}{an object obtained from \code{\link{lmInfl}}.}
\item{mult}{an object obtained from \code{\link{lmMult}}.}
\item{...}{other plotting parameters.}
}

\value{
The corresponding plot.
}

\note{
Cut-off values for the different influence measures are the following:\cr
\bold{dfbeta slope}: \eqn{> 2/\sqrt{n}}\cr
\bold{dffits}: \eqn{> 2\sqrt{2/n}}\cr
\bold{covratio}: \eqn{|\mathrm{covr} - 1| > 6/n}\cr
\bold{Cook's D}: \eqn{> 6/n} or \eqn{> 3\cdot\overline{\mathrm{CookD}}}\cr
\bold{leverage}: \eqn{> 2\cdot\overline{\mathrm{lev}}}
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## See Examples in 'lmInfl' and 'lmMult'.
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

