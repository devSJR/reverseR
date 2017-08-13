\name{shinyInfl}
\alias{shinyInfl}
\encoding{latin1}

\title{Initializes a Shiny App for significance reversal testing in linear regression}

\description{
A comprehensive Shiny App that facilitates the import of user-supplied data and subsequent detailed testing of Leave-One-Out significance reversal in linear regression. Analysis results can be exported as plots and text. 
}

\usage{
shinyInfl() 
}

\arguments{
None.
}

\value{
The analysis including plots and analysis tables.
}

\author{
Andrej-Nikolai Spiess
}   

\examples{
## shinyInfl() # <= to initialize.

## To use the example data.
cat("Example files are in", system.file("shiny", package = "reverseR")) 
}   

\keyword{optimize}
\keyword{models}
\keyword{linear}

